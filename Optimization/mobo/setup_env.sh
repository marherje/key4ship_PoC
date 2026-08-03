#!/usr/bin/env bash
#
# Create the virtualenv this package runs in.
#
#   bash setup_env.sh            # -> Optimization/mobo/.venv
#   MOBO_VENV_DIR=/eos/... bash setup_env.sh
#
# The venv is layered on top of the key4hep python with --system-site-packages,
# so it keeps ROOT / PyYAML / numpy / scipy / pandas / matplotlib from the stack
# and only adds what the stack lacks (torch, botorch, gpytorch, hydra, plotly).
# Torch is installed CPU-only: lxplus has no GPU and the CUDA wheels are ~5x
# larger, which an AFS home cannot absorb.
#
# Afterwards:
#   source /path/to/.venv/bin/activate      # key4hep must be sourced first
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${HERE}/../.." && pwd)"
VENV="${MOBO_VENV_DIR:-${HERE}/.venv}"

# Space the venv needs, in KiB. torch-cpu unpacks to ~450 MB and the rest
# (botorch/gpytorch/hydra/plotly + the editable install) to ~150 MB; the margin
# covers pip's temporary unpack directory.
NEED_KIB=$((1200 * 1024))

echo "== mobo setup =="
echo "repo:  ${REPO}"
echo "venv:  ${VENV}"

# ── 1. key4hep ───────────────────────────────────────────────────────────────
if [ -z "${KEY4HEP_STACK:-}" ]; then
    echo "-- sourcing ${REPO}/init_key4ship.sh"
    # key4hep's own setup.sh reads unset variables, so -u has to come off for
    # the duration of the source or it aborts on `compiler: unbound variable`.
    set +u
    # shellcheck disable=SC1091
    source "${REPO}/init_key4ship.sh"
    set -u
else
    echo "-- key4hep already in this shell (${KEY4HEP_STACK})"
fi
command -v python3 >/dev/null || { echo "ERROR: no python3 after sourcing key4hep"; exit 1; }
echo "-- python: $(command -v python3) ($(python3 --version))"

# ── 2. room for it ───────────────────────────────────────────────────────────
free_kib_of() {
    # AFS reports its own quota, which is what actually bites here; df on an AFS
    # path shows the whole cell and would happily let the install fail halfway.
    local dir="$1" vol
    while [ ! -d "${dir}" ]; do dir="$(dirname "${dir}")"; done
    if command -v fs >/dev/null && vol=$(fs listquota "${dir}" 2>/dev/null); then
        echo "${vol}" | awk 'NR>1 && NF>=3 { print $2 - $3; exit }'
    else
        df -Pk "${dir}" | awk 'NR==2 { print $4 }'
    fi
}

FREE_KIB="$(free_kib_of "${VENV}")"
if [ -n "${FREE_KIB}" ] && [ "${FREE_KIB}" -lt "${NEED_KIB}" ]; then
    cat >&2 <<EOF
ERROR: not enough space for the venv at
           ${VENV}
       free: $((FREE_KIB / 1024)) MiB, needed: $((NEED_KIB / 1024)) MiB

       Free some space, or put the venv elsewhere (the loop driver is the only
       thing that uses it — worker nodes run the payload with the bare key4hep
       python, so the venv does not have to live inside the repo):

           MOBO_VENV_DIR=/eos/user/\${USER:0:1}/\${USER}/mobo-venv bash setup_env.sh
EOF
    exit 1
fi
echo "-- free space: $((FREE_KIB / 1024)) MiB"

# ── 3. venv ──────────────────────────────────────────────────────────────────
if [ -d "${VENV}" ]; then
    echo "-- reusing existing venv"
else
    python3 -m venv --system-site-packages "${VENV}"
fi
# shellcheck disable=SC1091
source "${VENV}/bin/activate"

# pip's cache defaults to ~/.cache/pip, i.e. the same AFS quota we just checked;
# a torch download would eat 300 MB of it for nothing.
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-${TMPDIR:-/tmp}/${USER}/mobo-pip-cache}"
mkdir -p "${PIP_CACHE_DIR}"

python3 -m pip install --upgrade pip >/dev/null

# ── 4. torch first, from the CPU-only index ──────────────────────────────────
if python3 -c "import torch" 2>/dev/null; then
    echo "-- torch already present ($(python3 -c 'import torch; print(torch.__version__)'))"
else
    echo "-- installing torch (CPU-only wheel)"
    python3 -m pip install --index-url https://download.pytorch.org/whl/cpu torch
fi

# ── 5. the package itself ────────────────────────────────────────────────────
echo "-- installing mobo (editable) + dev extras"
python3 -m pip install -e "${HERE}[dev]"

# key4hep's pathspec 0.12.1 is trimmed: it has no patterns/gitignore.py, which
# modern mypy imports on startup. pip considers the requirement satisfied by
# that copy and installs nothing, so force a complete one into the venv. It only
# takes effect together with the PYTHONPATH reordering in check.sh, since
# key4hep's site-packages are searched first.
python3 -m pip install -q --force-reinstall --no-deps "pathspec>=0.12,<1"

echo
echo "== done =="
echo "Use it with:"
echo "    source ${REPO}/init_key4ship.sh"
echo "    source ${VENV}/bin/activate"
echo "    mobo-run --help"
