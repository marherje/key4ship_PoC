#!/usr/bin/env bash
#
# Lint, type-check and test, the way CI would if there were CI.
#
#   bash check.sh            # ruff + mypy + the fast tests
#   bash check.sh --slow     # ... and the slow ones (GP quality, real payload)
#
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${HERE}"

if [ -z "${VIRTUAL_ENV:-}" ]; then
    echo "Activate the venv first:"
    echo "    source \$(git rev-parse --show-toplevel)/init_key4ship.sh"
    echo "    source ${HERE}/.venv/bin/activate"
    exit 1
fi

# key4hep exports its own site-packages on PYTHONPATH, and PYTHONPATH is
# searched *before* the venv. Any package the stack also ships therefore wins
# over the venv's copy — which for mypy is fatal, because key4hep's pathspec is
# trimmed and lacks the submodule mypy imports. Putting the venv first fixes it
# and is harmless for a static checker.
export PYTHONPATH="${VIRTUAL_ENV}/lib/python3.13/site-packages:${PYTHONPATH:-}"

status=0
run() {
    echo
    echo "== $* =="
    "$@" || status=1
}

# pytest and mypy go through `python -m`: pytest also exists in key4hep, and
# pip does not create .venv/bin/<tool> for anything it considers already
# satisfied by --system-site-packages, so a bare `pytest` would run key4hep's
# copy with key4hep's interpreter — which cannot see the editable install.
#
# ruff is the exception and must be called directly: it is a standalone binary,
# and `python -m ruff` only looks for it inside the venv's bin (where, for the
# same reason, it was never installed).
run ruff check src tests
run ruff format --check src tests
# ship/ and viz/ are out of strict typing on purpose: they import ROOT, dd4hep
# and config.py dynamically by path, which mypy cannot follow.
run python -m mypy src/mobo/core src/mobo/exec

if [ "${1:-}" = "--slow" ]; then
    run python -m pytest -q --cov=src/mobo/core --cov=src/mobo/exec --cov-report=term-missing
else
    run python -m pytest -q -m "not slow" --cov=src/mobo/core --cov-report=term-missing
fi

echo
[ ${status} -eq 0 ] && echo "== all clean ==" || echo "== something failed =="
exit ${status}
