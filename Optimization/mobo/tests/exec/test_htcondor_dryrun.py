"""The Condor executor, without a scheduler.

Three things are worth testing without a cluster, and they are exactly the three
that break silently on one:

* the generated `job.sub` / `job.sh` (checked against golden files, so a change
  to either has to be deliberate);
* the parsing of `condor_submit` and `condor_q` output (recorded verbatim);
* the state machine — held, vanished, timed out — driven by sentinels written
  by hand.

Plus the parity check: the same wrapper run locally produces the same result as
the local executor, which is what makes "it worked locally" mean something.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import pytest

from mobo.core.types import Trial
from mobo.exec.base import ExitInfo, Handle
from mobo.exec.htcondor import CondorConfig, CondorExecutor, parse_cluster_id, parse_queue
from mobo.exec.local import LocalExecutor

GOLDEN = Path(__file__).resolve().parent / "golden"


def trial(n: int = 42, attempt: int = 0) -> Trial:
    return Trial(trial_id=n, params={"a": 1.0}, unit_x=[0.5], seed=7, attempt=attempt)


def executor(**kwargs) -> CondorExecutor:
    kwargs.setdefault("init_script", "/repo/init_key4ship.sh")
    return CondorExecutor(CondorConfig(**kwargs))


# ── generated job files ──────────────────────────────────────────────────────


def normalise(text: str, workdir: Path) -> str:
    """Golden files cannot contain a tmp_path, so put a placeholder back."""
    return text.replace(str(workdir), "/WORKDIR")


@pytest.mark.parametrize("name", ["job.sh", "job.sub"])
def test_generated_job_files_match_the_golden_copies(tmp_path, name, request):
    ex = executor()
    cmd = ["python3", "/pkg/payload.py", "--workdir", str(tmp_path)]
    ex.write_job_files(trial(), cmd, tmp_path)

    produced = normalise((tmp_path / name).read_text(), tmp_path)
    golden = GOLDEN / name
    if request.config.getoption("--update-golden", default=False):
        golden.parent.mkdir(parents=True, exist_ok=True)
        golden.write_text(produced)
    assert produced == golden.read_text(), (
        f"{name} changed; re-run with --update-golden if that was intended"
    )


def test_the_wrapper_is_executable_and_sources_the_stack(tmp_path):
    ex = executor()
    ex.write_job_files(trial(), ["python3", "payload.py"], tmp_path)
    wrapper = tmp_path / "job.sh"
    assert wrapper.stat().st_mode & 0o111
    text = wrapper.read_text()
    assert "init_key4ship.sh" in text
    assert f"cd {tmp_path}" in text
    # The safety net that turns a payload that never started into a classified
    # failure instead of a job the loop has to time out.
    assert "left no sentinel" in text


def test_submit_file_carries_the_credential_and_the_flavour(tmp_path):
    ex = executor(flavour="tomorrow", request_cpus=4, request_memory="8GB")
    ex.write_job_files(trial(), ["true"], tmp_path)
    text = (tmp_path / "job.sub").read_text()
    assert '+JobFlavour             = "tomorrow"' in text
    assert "request_cpus            = 4" in text
    assert "request_memory          = 8GB" in text
    # Without a forwarded credential the job cannot write into AFS at all.
    assert "MY.SendCredential       = true" in text


def test_extra_submit_lines_are_passed_through(tmp_path):
    ex = executor(
        extra_submit={"+AccountingGroup": '"group_u_SHIP"'}, send_credential=False
    )
    ex.write_job_files(trial(), ["true"], tmp_path)
    text = (tmp_path / "job.sub").read_text()
    assert '+AccountingGroup        = "group_u_SHIP"' in text
    assert "SendCredential" not in text


def test_arguments_with_spaces_are_quoted(tmp_path):
    ex = executor()
    ex.write_job_files(trial(), ["python3", "-c", "print('hi there')"], tmp_path)
    text = (tmp_path / "job.sh").read_text()
    assert """'print('\\''hi there'\\'')'""" in text


# ── parsing what condor says ─────────────────────────────────────────────────


def test_parse_cluster_id_from_recorded_output():
    verbose = "Submitting job(s).\n1 job(s) submitted to cluster 8271934.\n"
    assert parse_cluster_id(verbose) == "8271934.0"
    assert parse_cluster_id("8271935.0 - 8271935.0\n") == "8271935.0"
    assert parse_cluster_id("ERROR: something went wrong\n") is None


def test_parse_queue_from_recorded_output():
    recorded = "8271934.0 2\n8271935.0 1\n8271936.0 5\n"
    assert parse_queue(recorded) == {"8271934.0": 2, "8271935.0": 1, "8271936.0": 5}
    assert parse_queue("") == {}
    assert parse_queue("-- Schedd: bigbird18.cern.ch\n") == {}


# ── the state machine ────────────────────────────────────────────────────────


def adopt_handle(
    ex: CondorExecutor, workdir: Path, job_id="123.0", submitted_at=None
) -> Handle:
    """Put a job into the executor's book without submitting anything."""
    from mobo.exec.htcondor import _Job

    handle = Handle(
        trial_id=0,
        workdir=workdir,
        job_id=job_id,
        submitted_at=submitted_at or time.time(),
    )
    ex._jobs[0] = _Job(handle)
    return handle


def test_a_done_sentinel_ends_the_job_without_asking_condor(tmp_path, monkeypatch):
    ex = executor()
    adopt_handle(ex, tmp_path)
    (tmp_path / "DONE").write_text("ok\n")

    def explode(*_a, **_k):
        raise AssertionError("condor_q must not be consulted when a sentinel exists")

    monkeypatch.setattr(ex, "_queue_snapshot", explode)
    ((_, info),) = ex.poll()
    assert info.ok and info.reason == "done"
    assert ex.in_flight == []


def test_a_failed_sentinel_is_a_classified_failure(tmp_path, monkeypatch):
    ex = executor()
    adopt_handle(ex, tmp_path)
    (tmp_path / "FAILED").write_text(
        "job4_failed: k4run exited with code 1\n\ntraceback...\n"
    )
    monkeypatch.setattr(ex, "_queue_snapshot", lambda: {})
    ((_, info),) = ex.poll()
    assert not info.ok and info.reason == "failed"
    assert info.detail == "job4_failed: k4run exited with code 1"


def test_a_running_job_is_left_alone(tmp_path, monkeypatch):
    ex = executor()
    adopt_handle(ex, tmp_path)
    monkeypatch.setattr(ex, "_queue_snapshot", lambda: {"123.0": 2})
    assert ex.poll() == []
    assert len(ex.in_flight) == 1


def test_a_held_job_is_removed_and_reported(tmp_path, monkeypatch):
    ex = executor()
    handle = adopt_handle(ex, tmp_path)
    monkeypatch.setattr(ex, "_queue_snapshot", lambda: {"123.0": 5})
    removed = []
    monkeypatch.setattr(ex, "cancel", lambda h: removed.append(h.job_id))

    ((got, info),) = ex.poll()
    assert got is handle
    assert not info.ok and info.reason == "held"
    assert removed == ["123.0"]


def test_a_vanished_job_is_lost_only_after_the_grace_period(tmp_path, monkeypatch):
    """AFS can take seconds to show a file written on another machine.

    Declaring the job lost immediately would re-run trials that had in fact
    finished, which is both wasteful and a source of duplicate work.
    """
    ex = executor(grace_seconds=60.0)
    adopt_handle(ex, tmp_path)
    monkeypatch.setattr(ex, "_queue_snapshot", lambda: {})

    assert ex.poll() == [], "the first sighting only starts the clock"
    assert ex.poll() == [], "still inside the grace period"

    ex._jobs[0].missing_since = time.time() - 120.0
    ((_, info),) = ex.poll()
    assert not info.ok and info.reason == "lost"
    assert "without writing a sentinel" in info.detail


def test_a_sentinel_arriving_late_still_wins(tmp_path, monkeypatch):
    """The job left the queue, and only then did AFS show its DONE."""
    ex = executor(grace_seconds=0.0)
    adopt_handle(ex, tmp_path)
    monkeypatch.setattr(ex, "_queue_snapshot", lambda: {})
    assert ex.poll() == []

    (tmp_path / "DONE").write_text("ok\n")
    ((_, info),) = ex.poll()
    assert info.ok and info.reason == "done"


def test_an_unreachable_scheduler_does_not_kill_jobs(tmp_path, monkeypatch):
    """A condor_q that fails is a condor_q problem, not a trial problem."""
    ex = executor()
    adopt_handle(ex, tmp_path)
    monkeypatch.setattr(ex, "_queue_snapshot", lambda: None)
    assert ex.poll() == []
    assert len(ex.in_flight) == 1


def test_a_job_that_never_finishes_times_out(tmp_path, monkeypatch):
    ex = executor(timeout_hours=1.0)
    adopt_handle(ex, tmp_path, submitted_at=time.time() - 2 * 3600)
    monkeypatch.setattr(ex, "cancel", lambda h: None)
    monkeypatch.setattr(ex, "_queue_snapshot", lambda: {"123.0": 2})
    ((_, info),) = ex.poll()
    assert not info.ok and info.reason == "timeout"


def test_submit_failures_are_retried_then_raise(tmp_path, monkeypatch):
    ex = executor(submit_retries=3)
    calls = {"n": 0}

    def fake_run(cmd, **kwargs):
        calls["n"] += 1
        return subprocess.CompletedProcess(cmd, 1, "", "ERROR: schedd is busy")

    monkeypatch.setattr("mobo.exec.htcondor.subprocess.run", fake_run)
    monkeypatch.setattr("mobo.exec.htcondor.time.sleep", lambda _s: None)
    with pytest.raises(RuntimeError, match="condor_submit failed"):
        ex.submit(trial(), ["true"], tmp_path)
    assert calls["n"] == 3


def test_submit_records_the_cluster_id(tmp_path, monkeypatch):
    ex = executor()
    monkeypatch.setattr(
        "mobo.exec.htcondor.subprocess.run",
        lambda cmd, **kw: subprocess.CompletedProcess(
            cmd, 0, "1 job(s) submitted to cluster 991122.\n", ""
        ),
    )
    handle = ex.submit(trial(), ["true"], tmp_path)
    assert handle.job_id == "991122.0"
    assert (tmp_path / "job.sub").is_file()


def test_stale_sentinels_are_cleared_before_a_retry(tmp_path, monkeypatch):
    (tmp_path / "FAILED").write_text("the previous attempt\n")
    ex = executor()
    monkeypatch.setattr(
        "mobo.exec.htcondor.subprocess.run",
        lambda cmd, **kw: subprocess.CompletedProcess(cmd, 0, "cluster 5.\n", ""),
    )
    monkeypatch.setattr("mobo.exec.htcondor.parse_cluster_id", lambda _s: "5.0")
    ex.submit(trial(attempt=1), ["true"], tmp_path)
    assert not (tmp_path / "FAILED").exists()


# ── adoption after the driver died ───────────────────────────────────────────


def test_adopt_takes_back_a_job_still_in_the_queue(tmp_path, monkeypatch):
    ex = executor()
    monkeypatch.setattr(ex, "_queue_snapshot_uncached", lambda ids: {"777.0": 2})
    t = trial()
    t.workdir, t.job_id = str(tmp_path), "777.0"
    handle = ex.adopt(t)
    assert handle is not None and handle.job_id == "777.0"
    assert len(ex.in_flight) == 1


def test_adopt_takes_back_a_job_that_finished_while_we_were_away(tmp_path, monkeypatch):
    """Its result is already on disk; re-running it would pay twice."""
    ex = executor()
    (tmp_path / "DONE").write_text("ok\n")
    monkeypatch.setattr(ex, "_queue_snapshot_uncached", lambda ids: {})
    t = trial()
    t.workdir, t.job_id = str(tmp_path), "777.0"
    assert ex.adopt(t) is not None
    ((_, info),) = ex.poll()
    assert info.ok


def test_adopt_gives_up_on_a_job_that_is_really_gone(tmp_path, monkeypatch):
    ex = executor()
    monkeypatch.setattr(ex, "_queue_snapshot_uncached", lambda ids: {})
    t = trial()
    t.workdir, t.job_id = str(tmp_path), "777.0"
    assert ex.adopt(t) is None, "the loop must re-queue this one"


# ── RUN_LOCAL parity ─────────────────────────────────────────────────────────


def run_to_completion(ex, t, cmd, workdir, timeout=60.0) -> ExitInfo:
    ex.submit(t, cmd, workdir)
    deadline = time.time() + timeout
    while time.time() < deadline:
        finished = ex.poll()
        if finished:
            return finished[0][1]
        time.sleep(0.05)
    raise AssertionError("the job never finished")


def synthetic_cmd(workdir: Path) -> list[str]:
    payload = Path(__file__).resolve().parents[1] / "synthetic_payload.py"
    (workdir / "trial.json").write_text(json.dumps({"unit_x": [0.3, 0.7], "seed": 99}))
    return [sys.executable, str(payload), "--workdir", str(workdir)]


def test_the_condor_wrapper_gives_the_same_result_as_the_local_executor(tmp_path):
    """Same payload, same seed, two execution paths: the physics must not care.

    This is what makes the local executor a rehearsal rather than a separate
    code path that merely looks similar.
    """
    local_dir = tmp_path / "local"
    local_dir.mkdir()
    local = LocalExecutor(max_parallel=1)
    info = run_to_completion(local, trial(0), synthetic_cmd(local_dir), local_dir)
    assert info.ok, info.detail

    wrapper_dir = tmp_path / "wrapped"
    wrapper_dir.mkdir()
    # No key4hep on a test runner, so the wrapper's `source` must not be fatal;
    # pointing it at /dev/null exercises everything else in the script.
    condor = executor(run_local=True, init_script="/dev/null")
    info = run_to_completion(condor, trial(0), synthetic_cmd(wrapper_dir), wrapper_dir)
    assert info.ok, info.detail

    assert (
        json.loads((local_dir / "metrics.json").read_text())["f1"]
        == json.loads((wrapper_dir / "metrics.json").read_text())["f1"]
    )


def test_the_wrapper_reports_a_payload_that_never_started(tmp_path):
    """No sentinel from the payload -> the wrapper writes one. Not "lost"."""
    condor = executor(run_local=True, init_script="/dev/null")
    info = run_to_completion(condor, trial(0), ["false"], tmp_path)
    assert not info.ok and info.reason == "failed"
    assert "left no sentinel" in (info.detail or "")
