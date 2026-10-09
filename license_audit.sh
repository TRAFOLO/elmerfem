#!/bin/bash
# Phase 4 -- license audit gate: every shipped binary's import closure must be
# in the allowlist; forbidden names fail the build. Also stages license texts
# and SOURCE.txt into the install root.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="${ELMER_SRC:-$HERE}"
INSTALL="${ELMER_INSTALL:-$SRC/../elmer-install-win}"
# Direct solver is Elmer's VENDORED UMFPACK 4.4 (in-tree under $SRC/umfpack);
# MUMPS (CeCILL-C, GPL-incompatible) is never part of a distributed build.
# Limitation: audit detects a mumps DLL import but cannot detect a statically
# linked MUMPS — the build flags (build_msys2.sh WITH_Mumps=OFF) are the real guard.
# git may be absent from the MSYS2 login shell; fall back to Git-for-Windows.
# safe.directory avoids dubious-ownership refusals when the checkout belongs to
# another Windows user.
GIT="$(command -v git || echo '/c/Program Files/Git/cmd/git.exe')"
git() { "$GIT" -c safe.directory='*' "$@"; }
# Prefer the commit the BINARY was built from (banner "Rev:") over the checkout
# HEAD: after doc-only commits, HEAD drifts past the build commit and SOURCE.txt
# would point at a descendant. Fallback to HEAD (still a valid superset pointer).
BANNER_REV=$(env -i PATH="$(cygpath -w "$INSTALL/bin");C:\\Windows\\System32" "$INSTALL/bin/ElmerSolver.exe" 2>/dev/null </dev/null | grep -oE "Rev: [0-9a-f]+" | awk '{print $2}')
if [ -n "$BANNER_REV" ] && [ "$BANNER_REV" != "unknown" ]; then
  FORK_COMMIT=$(git -C "$SRC" rev-parse "$BANNER_REV" 2>/dev/null)
fi
[ -z "$FORK_COMMIT" ] && FORK_COMMIT=$(git -C "$SRC" rev-parse HEAD 2>/dev/null)
if [ -z "$FORK_COMMIT" ]; then
  echo "ERROR: could not resolve fork commit (git failed in this shell) -- SOURCE.txt would be incomplete"
  exit 1
fi
echo "SOURCE commit: $FORK_COMMIT (banner Rev: ${BANNER_REV:-n/a})"
FAIL=0

cd "$INSTALL"
ALL_IMPORTS=$(for f in bin/*.exe bin/*.dll share/elmersolver/lib/*.dll; do
  objdump -x "$f" 2>/dev/null | grep "DLL Name" | awk '{print $3}'; done | sort -u)

echo "=== aggregated import closure ==="; echo "$ALL_IMPORTS"

# forbidden anywhere in the closure (umfpack-as-DLL forbidden too: only the
# vendored in-tree copy may be used, which links statically and never imports)
for bad in mumps parmetis scotch umfpack suitesparse cholmod fftw readline gsl mkl pardiso hypre zoltan; do
  if echo "$ALL_IMPORTS" | grep -qi "$bad"; then echo "FORBIDDEN IMPORT: $bad"; FAIL=1; fi
done
# msmpi may be imported (system prerequisite) but must not be SHIPPED
if ls bin/msmpi.dll 2>/dev/null; then echo "FORBIDDEN FILE: bin/msmpi.dll must not be zipped"; FAIL=1; fi

# every shipped DLL file must be in the allowlist (catches accidental copies)
# libmsmpifec: since MSYS2 msmpi 10.1.1-23 the MS-MPI Fortran interface is in
# this DLL. MIT (Microsoft-MPI), so it ships with its license text;
# msmpi.dll itself stays a system prerequisite (checked above).
ALLOW='^(libelmersolver|libmatc|libfhuti|libmpi_stubs|libgcc_s_seh-1|libgfortran-5|libquadmath-0|libwinpthread-1|libgomp-1|libstdc\+\+-6|libopenblas|libscalapack|libarpack|libparpack|zlib1|libmsmpifec)\.dll$'
for f in bin/*.dll; do
  b=$(basename "$f")
  echo "$b" | grep -qE "$ALLOW" || { echo "UNLISTED SHIPPED DLL: $b (classify before shipping)"; FAIL=1; }
done

# --- stage license texts + SOURCE note ---
mkdir -p licenses
cp -r "$SRC/license_texts"/* licenses/ 2>/dev/null
# UMFPACK 4.4 notice (vendored in the fork; Davis notice-retention license with a
# user-documentation citation duty — see THIRD_PARTY_NOTICES.md in licenses/)
U_LIC=$(find "$SRC/umfpack" -maxdepth 4 -iname "*license*" 2>/dev/null | head -1)
[ -n "$U_LIC" ] && cp "$U_LIC" licenses/LICENSE.UMFPACK.txt \
  || echo "note: no standalone UMFPACK license file found in-tree (covered by THIRD_PARTY_NOTICES.md)"
rm -f licenses/LICENSE.MUMPS-CeCILL-C.txt  # never ship a MUMPS license text
for pkg in scalapack openblas zlib arpack-ng; do
  L=$(ls /c/msys64/ucrt64/share/licenses/$pkg/* 2>/dev/null | head -1)
  [ -n "$L" ] && cp "$L" "licenses/LICENSE.$pkg.txt"
done
# The msmpi package carries no license file of its own; the MIT text is kept here.
if [ -f bin/libmsmpifec.dll ]; then
  cp "$HERE/trafolo_bundle/licenses/LICENSE.msmpi.txt" licenses/LICENSE.msmpi.txt \
    || { echo "ERROR: MS-MPI license text missing for bin/libmsmpifec.dll"; FAIL=1; }
fi

cat > SOURCE.txt <<EOF
ElmerFEM solver bundle for TRAFOLO -- source and license notes
==============================================================
Built: $(date -u +%Y-%m-%d)

Elmer sources (GPL-2.0+/LGPL-2.1, see licenses/):
  TRAFOLO fork of ElmerFEM, branch 'trafolo'
  commit $FORK_COMMIT
  Repository: https://github.com/TRAFOLO/elmerfem
  Includes TRAFOLO modules ProcessFields, LoadFields, htc_udf
  (fem/src/modules/, GPL-2.0+). Corresponding sources are available
  from the repository at the commit above.

UMFPACK 4.4 (direct solver; vendored in-tree at umfpack/ of the fork,
  Davis notice-retention license with a user-documentation citation duty):
  part of the corresponding sources above - see licenses/ and
  licenses/THIRD_PARTY_NOTICES.md. No MUMPS in this bundle (CeCILL-C is
  GPL-incompatible and banned for distribution).

Dynamic third-party runtimes in bin/:
  OpenBLAS (BSD-3), GCC runtimes libgcc/libgfortran/libquadmath
  (GPL with runtime exception), libwinpthread (MIT/BSD), the MS-MPI
  Fortran interface libmsmpifec (MIT, built by MSYS2 from Microsoft-MPI)
  when present, and any further DLLs enumerated by the audit above (see
  licenses/ for their texts).

System prerequisite (NOT included): Microsoft MPI runtime (msmpi.dll).
EOF
echo "SOURCE.txt + licenses/ staged."

[ $FAIL = 0 ] && echo "LICENSE AUDIT: CLEAN" || echo "LICENSE AUDIT: FAILED"
exit $FAIL
