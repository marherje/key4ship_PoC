"""LocalExecutor: the sentinel protocol, capacity, timeouts and cancellation."""

from __future__ import annotations

import sys
import time

import pytest

from mobo.core.types import Trial
from mobo.exec.base import read_sentinel
from mobo.exec.local import LocalExecutor


def trial(n: int = 0) -> Trial:
    return Trial(trial_id=n, params={}, unit_x=[0.5], seed=1)


def wait_for_poll(executor, timeout=30.0):
    """Poll until something finishes, or give up."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        finished = executor.poll()
        if finished:
            return finished
        time.sleep(0.02)
    raise AssertionError("nothing finished within the timeout")


def sentinel_cmd(kind: str, sleep: float = 0.0) -> list[str]:
    return [
        sys.executable,
        "-c",
        f"import time, pathlib; time.sleep({sleep}); "
        f"pathlib.Path('{kind}').write_text('from the payload\\n')",
    ]


def test_done_sentinel_is_success(tmp_path):
    ex = LocalExecutor(max_parallel=1)
    handle = ex.submit(trial(), sentinel_cmd("DONE"), tmp_path / "t0")
    ((got, info),) = wait_for_poll(ex)
    assert got.trial_id == 0
    assert info.ok and info.reason == "done"
    assert info.detail == "from the payload"
    assert handle.job_id.startswith("pid:")


def test_failed_sentinel_is_a_classified_failure(tmp_path):
    ex = LocalExecutor(max_parallel=1)
    ex.submit(trial(), sentinel_cmd("FAILED"), tmp_path / "t0")
    ((_, info),) = wait_for_poll(ex)
    assert not info.ok and info.reason == "failed"
    assert info.detail == "from the payload"


def test_a_process_that_dies_without_a_sentinel_is_lost(tmp_path):
    """Not the same as a failure: the payload never got to say anything."""
    ex = LocalExecutor(max_parallel=1)
    ex.submit(trial(), [sys.executable, "-c", "import sys; sys.exit(7)"], tmp_path / "t0")
    ((_, info),) = wait_for_poll(ex)
    assert not info.ok and info.reason == "lost"
    assert "code 7" in info.detail


def test_the_sentinel_wins_over_the_exit_code(tmp_path):
    """ddsim is known to crash in teardown after writing a good output file.

    The payload decides whether its work was good; a non-zero exit code after a
    DONE is the runtime's problem, not the trial's.
    """
    ex = LocalExecutor(max_parallel=1)
    cmd = [
        sys.executable,
        "-c",
        "import pathlib, sys; pathlib.Path('DONE').write_text('ok\\n'); sys.exit(1)",
    ]
    ex.submit(trial(), cmd, tmp_path / "t0")
    ((_, info),) = wait_for_poll(ex)
    assert info.ok and info.reason == "done"


def test_stale_sentinels_are_cleared_on_submit(tmp_path):
    """A retry must not read the previous attempt's verdict."""
    workdir = tmp_path / "t0"
    workdir.mkdir()
    (workdir / "FAILED").write_text("from a previous attempt\n")

    ex = LocalExecutor(max_parallel=1)
    ex.submit(trial(), sentinel_cmd("DONE", sleep=0.2), workdir)
    assert read_sentinel(workdir) is None, "the stale sentinel survived the submit"
    ((_, info),) = wait_for_poll(ex)
    assert info.ok


def test_capacity_is_enforced(tmp_path):
    ex = LocalExecutor(max_parallel=2)
    for i in range(2):
        ex.submit(trial(i), sentinel_cmd("DONE", sleep=2.0), tmp_path / f"t{i}")
    assert len(ex.in_flight) == 2
    with pytest.raises(RuntimeError, match="full"):
        ex.submit(trial(2), sentinel_cmd("DONE"), tmp_path / "t2")
    ex.shutdown(cancel_running=True)


def test_cancel_stops_a_running_trial(tmp_path):
    ex = LocalExecutor(max_parallel=1)
    handle = ex.submit(trial(), sentinel_cmd("DONE", sleep=30.0), tmp_path / "t0")
    assert len(ex.in_flight) == 1
    ex.cancel(handle)
    assert ex.in_flight == []
    assert not (tmp_path / "t0" / "DONE").exists()


def test_timeout_kills_and_reports(tmp_path):
    ex = LocalExecutor(max_parallel=1, timeout_hours=0.5 / 3600.0)  # 0.5 s
    ex.submit(trial(), sentinel_cmd("DONE", sleep=30.0), tmp_path / "t0")
    ((_, info),) = wait_for_poll(ex, timeout=15.0)
    assert not info.ok and info.reason == "timeout"
    assert ex.in_flight == []


def test_local_executor_cannot_adopt(tmp_path):
    """Its children die with the driver, so a resume must re-queue them."""
    ex = LocalExecutor(max_parallel=1)
    assert ex.adopt(trial()) is None
