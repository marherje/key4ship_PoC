"""SQLite-backed trial store — the single source of truth of an experiment.

The store owns trial ids. That is not bookkeeping pedantry: the old controller
derived the next id by globbing the output directory, which is a race as soon as
two things run at once and silently reuses an id when a directory is cleaned up.
Here an id is allocated inside a `BEGIN IMMEDIATE` transaction, so concurrent
writers (threads, several loops, a resumed run) can never collide.

Resuming is `TrialStore(path)` on an existing file plus `completed()` and
`in_flight()`: nothing about the run lives only in the driver's memory.
"""

from __future__ import annotations

import json
import sqlite3
import threading
import time
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from .types import Result, Trial, TrialStatus

SCHEMA = """
CREATE TABLE IF NOT EXISTS trials (
    trial_id     INTEGER PRIMARY KEY,
    params       TEXT    NOT NULL,
    unit_x       TEXT    NOT NULL,
    seed         INTEGER NOT NULL,
    status       TEXT    NOT NULL,
    metrics      TEXT,
    error        TEXT,
    workdir      TEXT,
    job_id       TEXT,
    attempt      INTEGER NOT NULL DEFAULT 0,
    tag          TEXT,
    created_at   REAL,
    submitted_at REAL,
    finished_at  REAL
);
CREATE INDEX IF NOT EXISTS trials_status ON trials(status);

CREATE TABLE IF NOT EXISTS meta (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
"""

_COLUMNS = (
    "trial_id, params, unit_x, seed, status, metrics, error, workdir, job_id, "
    "attempt, tag, created_at, submitted_at, finished_at"
)


class TrialStore:
    """Persistent trial table. Safe to share between threads and processes."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._local = threading.local()
        with self._conn() as conn:
            conn.executescript(SCHEMA)

    # ── connection plumbing ──────────────────────────────────────────────────

    @contextmanager
    def _conn(self) -> Iterator[sqlite3.Connection]:
        """One connection per thread; sqlite3 objects are not thread-safe."""
        conn = getattr(self._local, "conn", None)
        if conn is None:
            conn = sqlite3.connect(str(self.path), timeout=60.0, isolation_level=None)
            conn.row_factory = sqlite3.Row
            # WAL lets readers work while a writer holds the lock, which is what
            # `mobo-status` on a live run needs.
            conn.execute("PRAGMA journal_mode=WAL")
            conn.execute("PRAGMA busy_timeout=60000")
            conn.execute("PRAGMA synchronous=FULL")
            self._local.conn = conn
        yield conn

    @contextmanager
    def _tx(self) -> Iterator[sqlite3.Connection]:
        """A write transaction that takes the lock up front (BEGIN IMMEDIATE).

        Without IMMEDIATE, two writers both read the current max id under a
        shared lock and then fight over the upgrade — exactly the duplicate-id
        bug this store exists to prevent.
        """
        with self._conn() as conn:
            conn.execute("BEGIN IMMEDIATE")
            try:
                yield conn
            except BaseException:
                conn.execute("ROLLBACK")
                raise
            conn.execute("COMMIT")

    def close(self) -> None:
        conn = getattr(self._local, "conn", None)
        if conn is not None:
            conn.close()
            self._local.conn = None

    # ── rows <-> trials ──────────────────────────────────────────────────────

    @staticmethod
    def _to_trial(row: sqlite3.Row) -> Trial:
        return Trial(
            trial_id=int(row["trial_id"]),
            params=json.loads(row["params"]),
            unit_x=json.loads(row["unit_x"]),
            seed=int(row["seed"]),
            status=TrialStatus(row["status"]),
            metrics=json.loads(row["metrics"]) if row["metrics"] else None,
            error=row["error"],
            workdir=row["workdir"],
            job_id=row["job_id"],
            attempt=int(row["attempt"]),
            tag=row["tag"],
            created_at=row["created_at"],
            submitted_at=row["submitted_at"],
            finished_at=row["finished_at"],
        )

    @staticmethod
    def _values(trial: Trial) -> tuple:
        return (
            trial.trial_id,
            json.dumps(trial.params, sort_keys=True),
            json.dumps(list(trial.unit_x)),
            int(trial.seed),
            trial.status.value,
            json.dumps(trial.metrics, sort_keys=True) if trial.metrics else None,
            trial.error,
            trial.workdir,
            trial.job_id,
            int(trial.attempt),
            trial.tag,
            trial.created_at,
            trial.submitted_at,
            trial.finished_at,
        )

    # ── writing ──────────────────────────────────────────────────────────────

    def create(
        self,
        params: dict[str, Any],
        unit_x: Sequence[float],
        seed: int,
        tag: str | None = None,
    ) -> Trial:
        """Allocate the next id and insert a PROPOSED trial. Atomic."""
        trial = Trial(
            trial_id=-1,
            params=dict(params),
            unit_x=[float(u) for u in unit_x],
            seed=int(seed),
            tag=tag,
            created_at=time.time(),
        )
        with self._tx() as conn:
            row = conn.execute(
                "SELECT COALESCE(MAX(trial_id) + 1, 0) FROM trials"
            ).fetchone()
            trial.trial_id = int(row[0])
            conn.execute(
                f"INSERT INTO trials ({_COLUMNS}) VALUES ({','.join('?' * 14)})",
                self._values(trial),
            )
        return trial

    def add(self, trial: Trial) -> Trial:
        """Insert a trial that already knows its id (resume, tests)."""
        with self._tx() as conn:
            conn.execute(
                f"INSERT INTO trials ({_COLUMNS}) VALUES ({','.join('?' * 14)})",
                self._values(trial),
            )
        return trial

    def update(self, trial: Trial) -> None:
        """Write a whole trial back."""
        with self._tx() as conn:
            conn.execute(
                "UPDATE trials SET params=?, unit_x=?, seed=?, status=?, metrics=?, "
                "error=?, workdir=?, job_id=?, attempt=?, tag=?, created_at=?, "
                "submitted_at=?, finished_at=? WHERE trial_id=?",
                self._values(trial)[1:] + (trial.trial_id,),
            )

    def mark(self, trial: Trial, status: TrialStatus, **fields: Any) -> Trial:
        """Set the status (and optionally workdir/job_id/attempt) of a trial.

        Timestamps follow from the status, so callers never set them by hand.
        """
        trial.status = status
        for key, value in fields.items():
            setattr(trial, key, value)
        if status is TrialStatus.SUBMITTED and trial.submitted_at is None:
            trial.submitted_at = time.time()
        if status.is_terminal:
            trial.finished_at = time.time()
        self.update(trial)
        return trial

    def complete(self, trial: Trial, result: Result) -> Trial:
        """Fold an evaluation result in and persist it."""
        trial.apply(result)
        trial.finished_at = time.time()
        self.update(trial)
        return trial

    # ── reading ──────────────────────────────────────────────────────────────

    def get(self, trial_id: int) -> Trial:
        with self._conn() as conn:
            row = conn.execute(
                f"SELECT {_COLUMNS} FROM trials WHERE trial_id=?", (trial_id,)
            ).fetchone()
        if row is None:
            raise KeyError(f"no trial {trial_id} in {self.path}")
        return self._to_trial(row)

    def _select(self, where: str = "", args: Sequence[Any] = ()) -> list[Trial]:
        with self._conn() as conn:
            rows = conn.execute(
                f"SELECT {_COLUMNS} FROM trials {where} ORDER BY trial_id", tuple(args)
            ).fetchall()
        return [self._to_trial(r) for r in rows]

    def all(self) -> list[Trial]:
        return self._select()

    def by_status(self, *statuses: TrialStatus) -> list[Trial]:
        if not statuses:
            return self.all()
        marks = ",".join("?" * len(statuses))
        return self._select(f"WHERE status IN ({marks})", [s.value for s in statuses])

    def completed(self) -> list[Trial]:
        return self.by_status(TrialStatus.COMPLETED)

    def in_flight(self) -> list[Trial]:
        """Submitted or running: what a resume has to re-adopt or re-queue."""
        return self.by_status(TrialStatus.SUBMITTED, TrialStatus.RUNNING)

    def pending(self) -> list[Trial]:
        """Everything not finished, in-flight ones included."""
        return self.by_status(
            TrialStatus.PROPOSED, TrialStatus.SUBMITTED, TrialStatus.RUNNING
        )

    def counts(self) -> dict[str, int]:
        with self._conn() as conn:
            rows = conn.execute(
                "SELECT status, COUNT(*) FROM trials GROUP BY status"
            ).fetchall()
        return {r[0]: int(r[1]) for r in rows}

    def __len__(self) -> int:
        with self._conn() as conn:
            return int(conn.execute("SELECT COUNT(*) FROM trials").fetchone()[0])

    def next_id(self) -> int:
        """The id `create()` would hand out next. Informational only."""
        with self._conn() as conn:
            return int(
                conn.execute(
                    "SELECT COALESCE(MAX(trial_id)+1, 0) FROM trials"
                ).fetchone()[0]
            )

    # ── metadata ─────────────────────────────────────────────────────────────

    def set_meta(self, key: str, value: Any) -> None:
        with self._tx() as conn:
            conn.execute(
                "INSERT INTO meta(key, value) VALUES(?, ?) "
                "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
                (key, json.dumps(value)),
            )

    def get_meta(self, key: str, default: Any = None) -> Any:
        with self._conn() as conn:
            row = conn.execute("SELECT value FROM meta WHERE key=?", (key,)).fetchone()
        return default if row is None else json.loads(row[0])

    def all_meta(self) -> dict[str, Any]:
        with self._conn() as conn:
            rows = conn.execute("SELECT key, value FROM meta").fetchall()
        return {r[0]: json.loads(r[1]) for r in rows}

    # ── export ───────────────────────────────────────────────────────────────

    def to_dataframe(self) -> Any:
        """One row per trial, with a column set that does not depend on the rows.

        The old CSV writer derived its header from whichever batch it was
        appending, so a run that changed its parameter set corrupted the file.
        Here the columns are the union over the whole table, sorted, and a
        trial that lacks a key gets NaN.
        """
        import pandas as pd

        trials = self.all()
        param_keys = sorted({k for t in trials for k in t.params})
        metric_keys = sorted({k for t in trials for k in (t.metrics or {})})

        records = []
        for t in trials:
            rec: dict[str, Any] = {
                "trial_id": t.trial_id,
                "status": t.status.value,
                "tag": t.tag,
                "seed": t.seed,
                "attempt": t.attempt,
                "error": t.error,
                "workdir": t.workdir,
                "job_id": t.job_id,
                "created_at": t.created_at,
                "submitted_at": t.submitted_at,
                "finished_at": t.finished_at,
            }
            for k in param_keys:
                rec[k] = t.params.get(k)
            for k in metric_keys:
                # A metric named like a parameter would otherwise overwrite it.
                rec[k if k not in param_keys else f"{k}_metric"] = (t.metrics or {}).get(
                    k
                )
            records.append(rec)

        columns = [
            "trial_id",
            "status",
            "tag",
            "seed",
            "attempt",
            *param_keys,
            *[k if k not in param_keys else f"{k}_metric" for k in metric_keys],
            "error",
            "workdir",
            "job_id",
            "created_at",
            "submitted_at",
            "finished_at",
        ]
        return pd.DataFrame.from_records(records, columns=columns)
