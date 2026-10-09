#!/bin/bash
# The TRAFOLO bundle build in one command (MSYS2 UCRT64 shell, Microsoft MPI runtime installed):
# MUMPS 4.10.0, configure, compile, install, prune. Paths and overrides as in the scripts it runs.
set -e
SRC="$(cd "$(dirname "$0")" && pwd)"
export ELMER_SRC="$SRC"
export MUMPS410_INSTALL="${MUMPS410_INSTALL:-$SRC/../mumps410-install}"
export MUMPS_PREFIX="$MUMPS410_INSTALL"
bash "$SRC/mumps410/build_mumps410.sh"
bash "$SRC/build_msys2.sh"
ninja -C "${ELMER_BUILD:-$SRC/../elmer-build-win}"
bash "$SRC/deploy_msys2.sh"
bash "$SRC/prune_install.sh"
