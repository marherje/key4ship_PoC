"""The executor contract, and the sentinel-file protocol both executors share.

Completion is detected by a file the payload writes as its very last act —
`DONE` or `FAILED` in the trial's workdir — and not by asking the batch system
whether the job is still queued. That choice is the difference between a loop
that works and one that lies to itself:

* `condor_q` forgets a job seconds after it ends, so a poll that misses the
  window cannot tell "finished fine" from "never existed";
* the job's exit code is not visible at all until the log is parsed, and on AFS
  that log may be seconds stale;
* the same mechanism then works identically for local subprocesses, so the
  local executor is a real rehearsal of the Condor path rather than a different
  code path that happens to produce the same files.

The batch system is still consulted, but only as a watchdog: a job that has
vanished without leaving a sentinel is lost, and gets retried.
"""

from __future__ import annotations

import contextlib
import os
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path

from ..core.types import Trial

DONE = "DONE"
FAILED = "FAILED"


@dataclass
class Handle:
    """What an executor needs to keep track of one running trial."""

    trial_id: int
    workdir: Path
    job_id: str | None = None
    attempt: int = 0
    submitted_at: float = field(default_factory=time.time)

    @property
    def name(self) -> str:
        return f"t{self.trial_id:05d}"


@dataclass
class ExitInfo:
    """How a trial ended, from the executor's point of view.

    `ok` only means "the payload said it finished"; whether the *metrics* are
    any good is the evaluator's call in `collect`.
    """

    ok: bool
    reason: str  # "done" | "failed" | "lost" | "timeout" | "cancelled"
    detail: str | None = None


def refresh(workdir: Path) -> None:
    """Force a directory-cache refresh before looking for a sentinel.

    AFS caches directory listings per client, so a file written by a worker node
    can take seconds to become visible here. Listing the directory invalidates
    that entry; without it the loop sees a job as still running long after it
    finished, and eventually declares it lost.
    """
    with contextlib.suppress(OSError):
        os.listdir(workdir)


def read_sentinel(workdir: Path) -> ExitInfo | None:
    """DONE/FAILED in the workdir, or None if the trial is still running."""
    refresh(workdir)
    done = workdir / DONE
    failed = workdir / FAILED
    if done.is_file():
        return ExitInfo(ok=True, reason="done", detail=_first_line(done))
    if failed.is_file():
        return ExitInfo(ok=False, reason="failed", detail=_first_line(failed))
    return None


def clear_sentinels(workdir: Path) -> None:
    """Remove stale sentinels before a (re)submission."""
    for name in (DONE, FAILED):
        (workdir / name).unlink(missing_ok=True)


def _first_line(path: Path) -> str | None:
    try:
        text = path.read_text().strip()
    except OSError:
        return None
    return text.splitlines()[0][:500] if text else None


class Executor(ABC):
    """Runs a command in a workdir. Knows nothing about physics or about GPs."""

    @abstractmethod
    def submit(self, trial: Trial, cmd: list[str], workdir: Path) -> Handle:
        """Start the trial. Returns immediately — this is never blocking."""

    @abstractmethod
    def poll(self) -> list[tuple[Handle, ExitInfo]]:
        """Trials that finished since the last call. Non-blocking."""

    @abstractmethod
    def cancel(self, handle: Handle) -> None:
        """Stop a running trial (best effort)."""

    def adopt(self, trial: Trial) -> Handle | None:
        """Re-attach to a trial that was in flight when the driver died.

        Returns None if this executor cannot re-attach (the local executor
        cannot: its children died with it), in which case the loop re-queues.
        """
        del trial
        return None

    def shutdown(self, cancel_running: bool = False) -> None:
        """Release resources. Running work survives unless asked otherwise."""
        del cancel_running

    @property
    @abstractmethod
    def in_flight(self) -> list[Handle]:
        """Everything submitted and not yet reported by `poll`."""

    @property
    def capacity(self) -> int | None:
        """Max concurrent trials, or None for "as many as you like"."""
        return None
