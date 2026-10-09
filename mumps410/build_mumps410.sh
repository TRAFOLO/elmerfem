#!/bin/bash
# Build the public-domain MUMPS 4.10.0 (arithmetics s, c, d, z) for the TRAFOLO Elmer bundle, MSYS2 UCRT64,
# with METIS 5.1 orderings through metis4_on_metis5.c and without PORD. build_msys2.sh links the result into
# libelmersolver (next to UMFPACK) when MUMPS_PREFIX points to the install prefix written here.
#
# Usage (MSYS2 UCRT64 shell, Microsoft MPI runtime installed):
#   bash mumps410/build_mumps410.sh
#   MUMPS_PREFIX=<install prefix> bash build_msys2.sh
# Paths (env overrides; defaults are siblings of this checkout):
#   MUMPS410_WORK     build area, its MUMPS tree is WIPED on each run   (default: ../mumps410-work)
#   MUMPS410_TARBALL  source archive, downloaded when missing          (default: $MUMPS410_WORK/mumps_4.10.0.dfsg.orig.tar.gz)
#   MUMPS410_INSTALL  install prefix, WIPED on each run                (default: ../mumps410-install)
#   MPIEXEC_EXE       Microsoft MPI launcher for the adapter test      (default: C:/Program Files/Microsoft MPI/Bin/mpiexec.exe)
# Toolchain packages in addition to those of build_msys2.sh:
#   pacman -S --needed make mingw-w64-ucrt-x86_64-metis mingw-w64-ucrt-x86_64-scalapack
#
# Source: Debian's mumps_4.10.0.dfsg.orig.tar.gz (upstream MUMPS 4.10.0 without the user guide), pinned by its
# SHA-256 and fetched from the hash-addressed snapshot.debian.org archive. MUMPS 4.10.0 is public domain (its
# LICENSE file, installed as LICENSE.MUMPS-4.10.0.txt). MUMPS 5.x is CeCILL-C, which is not GPL-compatible: never
# link it into a distributed build. TRAFOLO's changes are the source fix below and the METIS adapter, both
# LGPL-2.1-or-later. Every stage is checked: the script stops at the first failure and prints PASS only at the end.
set -uo pipefail
KIT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$(cd "$KIT/.." && pwd)"
WORK="${MUMPS410_WORK:-$SRC/../mumps410-work}"
TARBALL="${MUMPS410_TARBALL:-$WORK/mumps_4.10.0.dfsg.orig.tar.gz}"
I="${MUMPS410_INSTALL:-$SRC/../mumps410-install}"
MPIEXEC="${MPIEXEC_EXE:-/c/Program Files/Microsoft MPI/Bin/mpiexec.exe}"
PFX="${MSYSTEM_PREFIX:-/ucrt64}"
S="$WORK/mumps-4.10.0.dfsg"
TARBALL_SHA256=c76339bba516b96a3021af93d9a31b0fbf5a68cfcd02c9578d665ba8018e4b11
TARBALL_URL=https://snapshot.debian.org/file/5971cd9ccd8c2d789a222140ca6a3b71d6f5229d
ADAPTER_CFLAGS=(-O2 -Wall -Wextra -Werror "-I$PFX/include")

die() { echo "FAIL: $*" >&2; exit 1; }
step() { echo "== $*"; }
# A folder this script wipes: never a root, the checkout or one of its parents, nor one holding the archive.
check_wipe_target() { # path label
  local abs
  mkdir -p "$1" || die "mkdir $1"
  abs="$(cd "$1" && pwd -P)" || die "cannot resolve $1"
  case "$abs" in /|/[a-zA-Z]) die "$2=$1 is a root directory";; esac
  case "$SRC_ABS/" in "$abs"/*) die "$2=$1 contains the source checkout";; esac
  case "$TARBALL_ABS" in "$abs"/*) die "$2=$1 contains the source archive $TARBALL";; esac
}
# Symbol table of a library. A failing or empty nm run is an inspection failure, never "symbol absent".
symbols() {
  local out
  out=$(nm "$1" 2>&1) || die "nm failed on $1: $(printf '%s' "$out" | head -3)"
  [ -n "$out" ] || die "nm printed nothing for $1"
  printf '%s\n' "$out"
}

step "toolchain"
[ "${MSYSTEM:-}" = UCRT64 ] || die "run this in the MSYS2 UCRT64 shell"
for tool in make gcc gfortran ar nm ranlib curl sha256sum tar sed grep; do
  command -v "$tool" > /dev/null || die "$tool not found (make: pacman -S make)"
done
[ -s "$PFX/lib/libmetis.a" ] || die "static METIS missing: pacman -S mingw-w64-ucrt-x86_64-metis"
[ -s "$PFX/lib/libscalapack.dll.a" ] || die "ScaLAPACK missing: pacman -S mingw-w64-ucrt-x86_64-scalapack"
[ -f "$MPIEXEC" ] || die "mpiexec not found at $MPIEXEC (install the Microsoft MPI runtime or set MPIEXEC_EXE)"
# gfortran is a native program: give make a Windows temporary folder.
T="$(cygpath -m /tmp)"

step "source archive"
mkdir -p "$WORK" || die "mkdir $WORK"
if [ ! -f "$TARBALL" ]; then
  curl -fsSL --retry 3 -o "$TARBALL.part" "$TARBALL_URL" || die "download $TARBALL_URL"
  mv "$TARBALL.part" "$TARBALL" || die "move the download"
fi
echo "$TARBALL_SHA256 *$TARBALL" | sha256sum -c - || die "tarball checksum"
SRC_ABS="$(cd "$SRC" && pwd -P)"
TARBALL_ABS="$(cd "$(dirname "$TARBALL")" && pwd -P)/$(basename "$TARBALL")"
check_wipe_target "$S" "the MUMPS source tree"
check_wipe_target "$I" MUMPS410_INSTALL
rm -rf -- "$S" "$I" || die "cannot clean the old build"
mkdir -p "$I/lib" "$I/include" || die "mkdir $I"
tar xzf "$TARBALL" -C "$WORK" || die "extract"
[ -f "$S/Makefile" ] || die "extracted tree incomplete"

step "source fix: argument kind of MUMPS_731"
# ?mumps_part5.F passes the default-integer LRHS_CNTR_MASTER_ROOT to MUMPS_731, whose first argument is
# INTEGER(8); on the allocation-error path that reads 8 bytes from a 4-byte variable.
OLD='                    CALL MUMPS_731(LRHS_CNTR_MASTER_ROOT,INFO(2))'
for a in s c d z; do
  f="$S/src/${a}mumps_part5.F"
  n=$(grep -c -x -F "$OLD" "$f")
  [ "$n" = 1 ] || die "$f: expected the MUMPS_731 call once, found $n"
  sed -i 's/^                    CALL MUMPS_731(LRHS_CNTR_MASTER_ROOT,INFO(2))$/                    CALL MUMPS_731(\n     \&                   int(LRHS_CNTR_MASTER_ROOT,8),INFO(2))/' "$f" \
    || die "edit $f"
  grep -q -F '     &                   int(LRHS_CNTR_MASTER_ROOT,8),INFO(2))' "$f" || die "$f: fix not applied"
done

step "configuration (METIS only, no PORD)"
cat > "$S/Makefile.inc" <<'EOF' || die "write Makefile.inc"
# MUMPS 4.10.0, PORD-free: METIS 5.1 through the TRAFOLO METIS-4 adapter + MUMPS built-in orderings.
# PFX (the MSYS2 prefix) is set on the make command line.
LPORDDIR =
IPORD    =
LPORD    =
IMETIS   = -I$(PFX)/include
LMETIS   = -lmetis
ORDERINGSF  = -Dmetis
ORDERINGSC  = $(ORDERINGSF)
LORDERINGS = $(LMETIS)
IORDERINGSF =
IORDERINGSC = $(IMETIS)
PLAT    =
LIBEXT  = .a
OUTC    = -o $(empty)
OUTF    = -o $(empty)
RM = rm -f
CC = gcc
FC = gfortran
FL = gfortran
empty :=
AR = ar vr $(empty)
RANLIB  = ranlib
SCALAP  = -lscalapack
INCPAR = -I$(PFX)/include
LIBPAR = $(SCALAP) -lmsmpi
INCSEQ = -I$(topdir)/libseq
LIBSEQ  =  -L$(topdir)/libseq -lmpiseq
LIBBLAS = -lopenblas
LIBOTHERS = -lpthread
CDEFS   = -DAdd_
OPTF    = -O3 -fallow-argument-mismatch -std=legacy -DALLOW_NON_INIT
OPTL    = -O3
OPTC    = -O3
INCS = $(INCPAR)
LIBS = $(LIBPAR)
LIBSEQNEEDED =
EOF

step "make s c d z"
t0=$(date +%s)
make -C "$S" -j1 PFX="$PFX" TMP="$T" TEMP="$T" TMPDIR="$T" s c d z > "$S/build.log" 2>&1 \
  || { tail -15 "$S/build.log"; die "make failed, see $S/build.log"; }
echo "make ok ($(( $(date +%s) - t0 )) s, $(grep -c 'Warning:' "$S/build.log") compiler warning lines)"
for a in s c d z; do [ -s "$S/lib/lib${a}mumps.a" ] || die "lib${a}mumps.a missing"; done
[ -s "$S/lib/libmumps_common.a" ] || die "libmumps_common.a missing"
if compgen -G "$S/PORD/lib/*.o" > /dev/null; then die "PORD objects were compiled"; fi
if grep -q -- "-Dpord" "$S/build.log"; then die "a compile line used -Dpord"; fi
grep -q -- "-Dmetis -c mumps_orderings.c" "$S/build.log" || die "mumps_orderings.c was not compiled with -Dmetis"
if grep -q -E "INTEGER\(4\)/INTEGER\(8\)|INTEGER\(8\)/INTEGER\(4\)" "$S/build.log"; then
  die "an INTEGER(4)/INTEGER(8) argument mismatch is still reported by the compiler"
fi

step "METIS adapter"
gcc "${ADAPTER_CFLAGS[@]}" -c "$KIT/metis4_on_metis5.c" -o "$S/lib/metis4_on_metis5.o" || die "adapter compile"

step "adapter tests"
TD="$S/adapter-tests"
mkdir -p "$TD" || die "mkdir $TD"
gcc "${ADAPTER_CFLAGS[@]}" "$KIT/tests/test_adapter_ok.c" "$S/lib/metis4_on_metis5.o" "$PFX/lib/libmetis.a" -lmsmpi \
  -o "$TD/test_adapter_ok.exe" || die "compile test_adapter_ok"
"$TD/test_adapter_ok.exe" || die "test_adapter_ok"
gcc "${ADAPTER_CFLAGS[@]}" "$KIT/tests/test_adapter_fail.c" "$KIT/tests/fake_metis_fails.c" "$S/lib/metis4_on_metis5.o" -lmsmpi \
  -o "$TD/test_adapter_fail.exe" || die "compile test_adapter_fail"
t0=$(date +%s)
timeout 120 "$MPIEXEC" -n 2 "$(cygpath -w "$TD/test_adapter_fail.exe")" > "$TD/test_adapter_fail.out" 2>&1
rc=$?
[ "$rc" -ne 124 ] || die "test_adapter_fail: the MPI job hung after a METIS failure"
[ "$rc" -ne 0 ] || die "test_adapter_fail: the MPI job exited 0 after a METIS failure"
grep -q "stopping the run" "$TD/test_adapter_fail.out" || die "test_adapter_fail: no failure message"
if grep -q "test_adapter_fail: ERROR" "$TD/test_adapter_fail.out"; then
  die "test_adapter_fail: $(grep 'test_adapter_fail: ERROR' "$TD/test_adapter_fail.out" | head -2)"
fi
echo "test_adapter_fail: OK (both ranks stopped, exit code $rc, $(( $(date +%s) - t0 )) s)"

step "libraries"
ar rs "$S/lib/libmumps_common.a" "$S/lib/metis4_on_metis5.o" || die "add the adapter to libmumps_common.a"
cp "$S"/lib/lib*mumps*.a "$I/lib/" || die "copy libraries"
cp "$S"/include/*.h "$I/include/" || die "copy headers"
# Elmer's FindMumps insists on a PORD library file; an empty archive satisfies it without linking PORD.
ar rc "$I/lib/libpord.a" || die "create empty libpord.a"
ranlib "$I/lib/libpord.a" || die "ranlib libpord.a"
members=$(ar t "$I/lib/libpord.a") || die "ar t libpord.a"
[ -z "$members" ] || die "libpord.a is not empty"
cp "$S/LICENSE" "$I/LICENSE.MUMPS-4.10.0.txt" || die "copy the MUMPS licence"
grep -q -i "public domain" "$I/LICENSE.MUMPS-4.10.0.txt" || die "the MUMPS LICENSE file does not say public domain"

step "symbol checks"
common=$(symbols "$I/lib/libmumps_common.a") || exit 1
grep -q " T metis_nodend_$" <<< "$common" || die "adapter symbol metis_nodend_ missing"
grep -q " T metis_nodewnd_$" <<< "$common" || die "adapter symbol metis_nodewnd_ missing"
grep -q " U METIS_NodeND$" <<< "$common" || die "adapter does not reference METIS_NodeND"
for a in s c d z; do
  arith=$(symbols "$I/lib/lib${a}mumps.a") || exit 1
  grep -q " U metis_nodend_$" <<< "$arith" || die "lib${a}mumps.a does not call the METIS ordering"
done
all=$(for f in "$I"/lib/lib*mumps*.a; do symbols "$f" || exit 1; done) || exit 1
if grep -q -i "parmetis" <<< "$all"; then die "ParMETIS symbols present"; fi
if grep -q -i -E "pord|space_ordering" <<< "$all"; then die "PORD symbols present: $(grep -i -E 'pord|space_ordering' <<< "$all" | head -3)"; fi
echo "adapter defined and used, no ParMETIS or PORD symbols"

step "build record"
commit=$(git -C "$KIT" rev-parse HEAD 2>/dev/null) || commit=unknown
changed=$(git -C "$KIT" status --porcelain -- . 2>/dev/null | wc -l)
{
  echo "kit=$KIT"
  echo "kit_commit=$commit"
  echo "kit_uncommitted_files=$changed"
  echo "tarball=$(basename "$TARBALL")"
  echo "tarball_sha256=$TARBALL_SHA256"
  echo "tarball_url=$TARBALL_URL"
  echo "metis=$(pacman -Q mingw-w64-ucrt-x86_64-metis 2>/dev/null || echo unknown), linked statically by build_msys2.sh"
  echo "scalapack=$(pacman -Q mingw-w64-ucrt-x86_64-scalapack 2>/dev/null || echo unknown)"
} > "$I/BUILD-INFO.txt" || die "write BUILD-INFO.txt"
cat "$I/BUILD-INFO.txt"
echo "PASS: PORD-free MUMPS 4.10.0 with the METIS adapter installed in $I"
