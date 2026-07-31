#!/bin/bash
# Run the FFT-enabled CPU maDGiCart build inside the CCR Apptainer/Singularity image.
#
# Defaults match Rachel's CCR layout. Override MADG_BASE, MADG_SRC, MADG_BUILD,
# or MADG_SIF if your OnDemand upload/build directories use different names.
set -euo pipefail

BASE="${MADG_BASE:-/projects/academic/kreyes3/BNL-M2DT/rlee/maDGiCart-github}"
SRC="${MADG_SRC:-$BASE/madgicart-fft}"
BUILD="${MADG_BUILD:-$BASE/madgicart-fft-gcc-build}"
SIF="${MADG_SIF:-/projects/academic/kreyes3/BNL-M2DT/maDGiCart/madg-gcc-amd64.sif}"

if [ -z "${FFTW_ROOT:-${EBROOTFFTW:-}}" ] && command -v module >/dev/null 2>&1; then
  module load "${MADG_FFTW_MODULE:-fftw/3.3.10}" >/dev/null 2>&1 || true
fi

FFTW_ROOT="${FFTW_ROOT:-${EBROOTFFTW:-}}"

[ -d "$SRC" ] || { echo "ERROR: no source directory $SRC" >&2; exit 1; }
[ -d "$BUILD" ] || { echo "ERROR: no build directory $BUILD" >&2; exit 1; }
[ -f "$SIF" ] || { echo "ERROR: no container image $SIF" >&2; exit 1; }
if [ -n "$FFTW_ROOT" ]; then
  [ -d "$FFTW_ROOT" ] || { echo "ERROR: FFTW root does not exist: $FFTW_ROOT" >&2; exit 1; }
  [ -f "$FFTW_ROOT/lib/libfftw3.so" ] || [ -f "$FFTW_ROOT/lib/libfftw3.so.3" ] || {
    echo "ERROR: FFTW shared library not found under $FFTW_ROOT/lib" >&2
    exit 1
  }
fi

APP=$(command -v apptainer || command -v singularity) || {
  echo "ERROR: no apptainer/singularity in PATH" >&2
  exit 2
}

container_args=(
  exec
  --env "LC_ALL=C"
  --bind "$SRC:/M2DT/maDGiCart/maDGiCart-CH"
  --bind "$BUILD:/M2DT/maDGiCart/maDGiCart-CH-gcc-build"
  --bind /scratch:/scratch
)

if [ -n "$FFTW_ROOT" ]; then
  container_args+=(--env "LD_LIBRARY_PATH=$FFTW_ROOT/lib")
  container_args+=(--bind "$FFTW_ROOT:$FFTW_ROOT")
else
  echo "WARNING: FFTW_ROOT/EBROOTFFTW is not set. Load an FFTW module before using dct_si." >&2
fi

exec "$APP" "${container_args[@]}" "$SIF" /M2DT/maDGiCart/maDGiCart-CH-gcc-build/maDGiCart "$@"
