"""Run trials as local subprocesses. Debugging tool and integration-test rig.

Exactly the same payload command as the Condor executor, run in the same
workdir, detected by the same sentinel files — so "it works locally" is
evidence about the Condor path and not just about itself.
"""

from __future__ import annotations

import contextlib
import logging
import os
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import IO

from ..core.types import Trial
from .base import Executor, ExitInfo, Handle, clear_sentinels, read_sentinel

log = logging.getLogger(__name__)


@dataclass
class _Run:
    """One live subprocess: its handle, the process, and the log file we own."""

    handle: Handle
    proc: subprocess.Popen
    stream: IO[str]


class LocalExecutor(Executor):
    """A bounded pool of subprocesses.

    The bound is the loop's `max_in_flight`, not this class's business, but it
    is enforced here too so that a misconfigured loop cannot fork-bomb lxplus.
    """

    def __init__(
        self,
        max_parallel: int = 2,
        timeout_hours: float | None = None,
        env: dict[str, str] | None = None,
    ) -> None:
        self.max_parallel = int(max_parallel)
        self.timeout_s = timeout_hours * 3600.0 if timeout_hours else None
        self.env = dict(env) if env else None
        self._running: dict[int, _Run] = {}

    # ── contract ─────────────────────────────────────────────────────────────

    def submit(self, trial: Trial, cmd: list[str], workdir: Path) -> Handle:
        if len(self._running) >= self.max_parallel:
            raise RuntimeError(
                f"local executor is full ({self.max_parallel} slots); "
                "the loop should not have asked"
            )
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        clear_sentinels(workdir)

        env = os.environ.copy()
        if self.env:
            env.update(self.env)

        log_path = workdir / "payload.log"
        handle = Handle(trial_id=trial.trial_id, workdir=workdir, attempt=trial.attempt)
        # Closed in _reap; the subprocess writes to it until it exits.
        stream = open(log_path, "w")  # noqa: SIM115
        proc = subprocess.Popen(
            [str(c) for c in cmd],
            cwd=str(workdir),
            env=env,
            stdout=stream,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        handle.job_id = f"pid:{proc.pid}"
        self._running[trial.trial_id] = _Run(handle, proc, stream)
        log.info("submitted %s locally (pid %d) in %s", handle.name, proc.pid, workdir)
        return handle

    def poll(self) -> list[tuple[Handle, ExitInfo]]:
        finished: list[tuple[Handle, ExitInfo]] = []
        for trial_id, run in list(self._running.items()):
            info = self._status(run.handle, run.proc)
            if info is None:
                continue
            self._reap(trial_id)
            finished.append((run.handle, info))
        return finished

    def _status(self, handle: Handle, proc: subprocess.Popen) -> ExitInfo | None:
        rc = proc.poll()
        if rc is None:
            if self.timeout_s and time.time() - handle.submitted_at > self.timeout_s:
                self._kill(proc)
                return ExitInfo(False, "timeout", f"exceeded {self.timeout_s:.0f} s")
            return None

        # The process is gone: the sentinel is authoritative if it exists, since
        # the payload writes it after everything it cares about is on disk.
        sentinel = read_sentinel(handle.workdir)
        if sentinel is not None:
            return sentinel
        return ExitInfo(
            False,
            "lost",
            f"exited with code {rc} without writing a sentinel; see payload.log",
        )

    def cancel(self, handle: Handle) -> None:
        run = self._running.get(handle.trial_id)
        if run is None:
            return
        self._kill(run.proc)
        self._reap(handle.trial_id)

    def shutdown(self, cancel_running: bool = False) -> None:
        if not cancel_running:
            return
        for handle in list(self.in_flight):
            self.cancel(handle)

    @property
    def in_flight(self) -> list[Handle]:
        return [run.handle for run in self._running.values()]

    @property
    def capacity(self) -> int | None:
        return self.max_parallel

    # ── plumbing ─────────────────────────────────────────────────────────────

    def _reap(self, trial_id: int) -> None:
        """Drop a finished trial and close its log stream (we own the fd)."""
        run = self._running.pop(trial_id, None)
        if run is not None:
            with contextlib.suppress(OSError):
                run.stream.close()

    @staticmethod
    def _kill(proc: subprocess.Popen) -> None:
        """Kill the whole process group: ddsim and k4run spawn children."""
        try:
            os.killpg(os.getpgid(proc.pid), 15)
        except (ProcessLookupError, PermissionError):
            proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(os.getpgid(proc.pid), 9)
            except (ProcessLookupError, PermissionError):
                proc.kill()
