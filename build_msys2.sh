#!/bin/bash
# Configure the solver-only (nogui) Windows build of the TRAFOLO Elmer fork,
# MSYS2 UCRT64 toolchain. This is the exact configuration of the solver bundle
# TRAFOLO distributes; it is published here so the distributed binaries can be
# rebuilt from this source (GPL-2 section 3).
#
# Usage (MSYS2 UCRT64 shell, from anywhere):
#   bash build_msys2.sh                        # configure
#   ninja -C "$ELMER_BUILD"                    # build
#   bash deploy_msys2.sh                       # install + runtime DLL closure
# Paths (env overrides; defaults are siblings of this checkout):
#   ELMER_SRC     this repo                    (default: directory of this script)
#   ELMER_BUILD   build dir, WIPED on each run (default: ../elmer-build-win)
#   ELMER_INSTALL install prefix               (default: ../elmer-install-win)
#   MUMPS_PREFIX  install prefix of mumps410/build_mumps410.sh (default: unset,
#                 no MUMPS); the distributed bundle sets it
#
# Toolchain (MSYS2 UCRT64):
#   pacman -S --needed mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-gcc-fortran \
#     mingw-w64-ucrt-x86_64-cmake mingw-w64-ucrt-x86_64-ninja \
#     mingw-w64-ucrt-x86_64-openblas mingw-w64-ucrt-x86_64-msmpi
#   with MUMPS also: mingw-w64-ucrt-x86_64-metis mingw-w64-ucrt-x86_64-scalapack
# plus the Microsoft MPI runtime ("C:\Program Files\Microsoft MPI").
#
# MUMPS: only the public-domain MUMPS 4.10.0 built by mumps410/build_mumps410.sh
#   (METIS 5.1 ordering linked statically, no PORD). MUMPS 5.x is CeCILL-C, whose
#   only declared compatibility is with CeCILL (its section 5.3.4), not the GPL:
#   never link it into a distributed build.
#
# WITH_UMFPACK=ON uses Elmer's VENDORED internal UMFPACK 4.4 (umfpack/ in this
# repo; LGPL-2.1-or-later, carries a user-documentation citation duty). It stays
# on with MUMPS: the "Circuit" preconditioner (VankaCreate.F90 CircuitPrec) uses
# UMFPACK unless a direct method is named, and switches to MUMPS only when
# UMFPACK is missing, which would change the iterative runs.
#
# BUNDLE_MSMPI_REDIST=ON stages redist/msmpisetup.exe (the Microsoft MPI
# redistributable installer) into the install tree. Requires the file
# msmpi_redist/msmpisetup.exe as a SIBLING of this checkout -- download
# MS-MPI 10.1.1 msmpisetup.exe from Microsoft and place it there. Drop the flag
# if you do not need the installer in the install tree.
set -e
SRC="${ELMER_SRC:-$(cd "$(dirname "$0")" && pwd)}"
BUILD="${ELMER_BUILD:-$SRC/../elmer-build-win}"
INSTALL="${ELMER_INSTALL:-$SRC/../elmer-install-win}"
PFX="${MSYSTEM_PREFIX:-/ucrt64}"

if [ -n "${MUMPS_PREFIX:-}" ]; then
  [ -s "$MUMPS_PREFIX/lib/libzmumps.a" ] && [ -s "$MUMPS_PREFIX/BUILD-INFO.txt" ] || {
    echo "ERROR: MUMPS_PREFIX=$MUMPS_PREFIX is not an install of mumps410/build_mumps410.sh"; exit 1; }
  MUMPS_FLAGS=(-DWITH_Mumps=ON "-DMUMPSROOT=$(cygpath -m "$MUMPS_PREFIX")"
    "-DMetis_LIBRARIES=$(cygpath -m "$PFX")/lib/libmetis.a" "-DMetis_INCLUDE_DIR=$(cygpath -m "$PFX")/include"
    "-DSCALAPACKROOT=$(cygpath -m "$PFX")")
else
  MUMPS_FLAGS=(-DWITH_Mumps=OFF)
fi

rm -rf "$BUILD"
mkdir -p "$BUILD"
cd "$BUILD"

cmake -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$INSTALL" \
  -DBLA_VENDOR=OpenBLAS \
  -DWITH_OpenMP=OFF \
  -DWITH_LUA=ON \
  -DWITH_MPI=ON \
  -DMPIEXEC_EXECUTABLE="/c/Program Files/Microsoft MPI/Bin/mpiexec.exe" \
  -DBUNDLE_MSMPI_REDIST=ON \
  -DWITH_UMFPACK=ON \
  "${MUMPS_FLAGS[@]}" \
  -DWITH_Hypre=OFF \
  -DWITH_Zoltan=OFF \
  -DWITH_ELMERGUI=OFF \
  -DWITH_ElmerIce=OFF \
  -DCMAKE_C_FLAGS="-Wno-error=incompatible-pointer-types" \
  "$SRC" > configure.log 2>&1

grep -E "MPI|BLAS|Mumps|Fortran compiler|Build type|install prefix|Configuring done|Generating done" configure.log | head -18

if [ -n "${MUMPS_PREFIX:-}" ]; then
  libs=$(grep "Mumps libraries:" configure.log) || { echo "ERROR: configure found no MUMPS"; exit 1; }
  case "$libs" in *"$(cygpath -m "$MUMPS_PREFIX")/lib/libzmumps.a"*) ;; *)
    echo "ERROR: configure did not pick MUMPS from $MUMPS_PREFIX"; exit 1;; esac
  case "$libs" in *"/lib/libmetis.a"*) ;; *) echo "ERROR: configure did not pick the static METIS"; exit 1;; esac
  grep -q "Checking if ParMetis library is needed by Mumps -- no" configure.log || {
    echo "ERROR: MUMPS would need ParMETIS"; exit 1; }
fi
