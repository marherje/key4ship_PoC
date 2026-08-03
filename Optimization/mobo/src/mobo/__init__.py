"""Asynchronous multi-objective Bayesian optimization of a detector geometry.

Layering (enforced by tests/test_architecture.py):

    core/   generic ask/tell machinery — knows about a search space, an
            objective vector and a Trial, and nothing about any detector.
    exec/   how a trial is run (subprocess pool, HTCondor). Takes a command
            and a workdir; knows no physics either.
    ship/   everything SHiP/SND-specific, behind the Evaluator interface.
    viz/    plots and the HTML report, driven by the trial store.

Porting the loop to another detector means writing a new `ship/`-equivalent
Evaluator plus a config; `core/` and `exec/` are untouched.
"""

__version__ = "0.1.0"
