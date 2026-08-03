"""The store owns the ids and the truth. Both claims are tested here."""

from __future__ import annotations

import threading

import pytest

from mobo.core.store import TrialStore
from mobo.core.types import Result, Trial, TrialStatus


def test_ids_start_at_zero_and_increment(store):
    ids = [store.create({"a": i}, [0.1 * i], seed=i).trial_id for i in range(5)]
    assert ids == [0, 1, 2, 3, 4]
    assert store.next_id() == 5


def test_concurrent_creates_have_no_gaps_or_duplicates(tmp_path):
    """The bug this store exists to prevent: two writers, one id.

    The old controller derived the next id by globbing the filesystem, which
    hands the same id to two concurrent callers.
    """
    store = TrialStore(tmp_path / "trials.db")
    n_threads, per_thread = 8, 25
    got: list[int] = []
    lock = threading.Lock()

    def worker(k: int) -> None:
        mine = [store.create({"w": k}, [0.5], seed=k).trial_id for _ in range(per_thread)]
        with lock:
            got.extend(mine)

    threads = [threading.Thread(target=worker, args=(k,)) for k in range(n_threads)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert sorted(got) == list(range(n_threads * per_thread))
    store.close()


def test_reopening_reproduces_every_trial(tmp_path):
    path = tmp_path / "trials.db"
    store = TrialStore(path)
    t = store.create({"a": 1.5, "n": 3}, [0.25, 0.75], seed=42, tag="sobol")
    store.mark(t, TrialStatus.SUBMITTED, workdir="/tmp/t0", job_id="12345.0")
    store.complete(t, Result.ok({"nhits": 10.0, "cost": 2.5}))
    store.set_meta("git_sha", "deadbeef")
    store.close()

    again = TrialStore(path)
    back = again.get(0)
    assert back.params == {"a": 1.5, "n": 3}
    assert back.unit_x == [0.25, 0.75]
    assert back.seed == 42
    assert back.tag == "sobol"
    assert back.status is TrialStatus.COMPLETED
    assert back.metrics == {"nhits": 10.0, "cost": 2.5}
    assert back.workdir == "/tmp/t0"
    assert back.job_id == "12345.0"
    assert back.submitted_at is not None and back.finished_at is not None
    assert again.get_meta("git_sha") == "deadbeef"
    again.close()


def test_status_queries(store):
    trials = [store.create({"a": i}, [0.0], seed=i) for i in range(5)]
    store.complete(trials[0], Result.ok({"m": 1.0}))
    store.complete(trials[1], Result.failed("sim_failed"))
    store.complete(trials[2], Result.infeasible("geometry_error"))
    store.mark(trials[3], TrialStatus.RUNNING)
    # trials[4] stays PROPOSED

    assert [t.trial_id for t in store.completed()] == [0]
    assert [t.trial_id for t in store.in_flight()] == [3]
    assert [t.trial_id for t in store.pending()] == [3, 4]
    assert store.counts() == {
        "COMPLETED": 1,
        "FAILED": 1,
        "INFEASIBLE": 1,
        "RUNNING": 1,
        "PROPOSED": 1,
    }
    assert len(store) == 5
    assert store.get(1).error == "sim_failed"


def test_add_preserves_an_explicit_id(store):
    store.add(Trial(trial_id=7, params={"a": 1}, unit_x=[0.5], seed=1))
    assert store.get(7).params == {"a": 1}
    assert store.next_id() == 8


def test_missing_trial_raises(store):
    with pytest.raises(KeyError):
        store.get(99)


def test_dataframe_columns_do_not_depend_on_the_rows(store):
    """The failure mode of the old results.csv: a header derived per batch.

    A trial with a different parameter set, or one that failed before producing
    metrics, must not shift or drop anyone else's columns.
    """
    t0 = store.create({"a": 1.0, "b": 2.0}, [0.1, 0.2], seed=1, tag="baseline")
    store.complete(t0, Result.ok({"nhits": 5.0, "cost": 1.0}))
    t1 = store.create({"a": 3.0, "c": 9.0}, [0.3, 0.4], seed=2)
    store.complete(t1, Result.failed("boom"))
    t2 = store.create({"a": 4.0, "b": 5.0}, [0.5, 0.6], seed=3)
    store.complete(t2, Result.ok({"nhits": 7.0, "extra": 3.0}))

    df = store.to_dataframe()
    assert list(df["trial_id"]) == [0, 1, 2]
    for column in ("a", "b", "c", "nhits", "cost", "extra"):
        assert column in df.columns
    assert df.loc[0, "nhits"] == 5.0
    assert df.loc[1, "status"] == "FAILED"
    assert df["b"].isna()[1]  # trial 1 never had a `b`
    assert df["cost"].isna()[2]  # trial 2 never reported a cost
    assert df.loc[0, "tag"] == "baseline"


def test_dataframe_disambiguates_a_metric_named_like_a_parameter(store):
    t = store.create({"cost": 1.0}, [0.5], seed=1)
    store.complete(t, Result.ok({"cost": 99.0}))
    df = store.to_dataframe()
    assert df.loc[0, "cost"] == 1.0  # the parameter
    assert df.loc[0, "cost_metric"] == 99.0  # the metric
