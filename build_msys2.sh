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
#
# Toolchain (MSYS2 UCRT64):
#   pacman -S --needed mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-gcc-fortran \
#     mingw-w64-ucrt-x86_64-cmake mingw-w64-ucrt-x86_64-ninja \
#     mingw-w64-ucrt-x86_64-openblas mingw-w64-ucrt-x86_64-msmpi
# plus the Microsoft MPI runtime ("C:\Program Files\Microsoft MPI").
#
# LICENSE CONSTRAINT - do NOT enable MUMPS for a distributed build:
#   MUMPS is CeCILL-C, whose only declared compatibility is with CeCILL (its
#   section 5.3.4), not the GPL. This bundle ships GPL-2 components (ElmerGrid,
#   most solver modules), so statically combining them with MUMPS is a GPL-
#   incompatibility risk for public distribution. Elmer's own
#   license_texts/ElmerLicensePolicy.md blesses only Umfpack/Hypre/(P)Arpack.
#
# WITH_UMFPACK=ON uses Elmer's VENDORED internal UMFPACK 4.4 (umfpack/ in this
# repo; Davis notice-retention license, carries a user-documentation citation
# duty). It is required: the "Circuit" preconditioner (VankaCreate.F90
# CircuitPrec) accepts ONLY umfpack or mumps and Fatals without one of them.
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
  -DWITH_Mumps=OFF \
  -DWITH_Hypre=OFF \
  -DWITH_Zoltan=OFF \
  -DWITH_ELMERGUI=OFF \
  -DWITH_ElmerIce=OFF \
  -DCMAKE_C_FLAGS="-Wno-error=incompatible-pointer-types" \
  "$SRC" > configure.log 2>&1

grep -E "MPI|BLAS|Fortran compiler|Build type|install prefix|Configuring done|Generating done" configure.log | head -15
