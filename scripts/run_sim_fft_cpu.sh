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

[ -d "$SRC" ] || { echo "ERROR: no source directory $SRC" >&2; exit 1; }
[ -d "$BUILD" ] || { echo "ERROR: no build directory $BUILD" >&2; exit 1; }
[ -f "$SIF" ] || { echo "ERROR: no container image $SIF" >&2; exit 1; }

APP=$(command -v apptainer || command -v singularity) || {
  echo "ERROR: no apptainer/singularity in PATH" >&2
  exit 2
}

exec "$APP" exec --env LC_ALL=C \
  --bind "$SRC:/M2DT/maDGiCart/maDGiCart-CH" \
  --bind "$BUILD:/M2DT/maDGiCart/maDGiCart-CH-gcc-build" \
  --bind /scratch:/scratch \
  "$SIF" /M2DT/maDGiCart/maDGiCart-CH-gcc-build/maDGiCart "$@"
