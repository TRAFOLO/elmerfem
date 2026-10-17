#!/bin/bash
# Finish install + make the install tree self-contained (runtime DLL closure).
# Run after build_msys2.sh + ninja; same env overrides and defaults.
#
# CPACK_BUNDLE_EXTRA_WINDOWS_DLLS=OFF: that option's install rules expect
# prebuilt sibling dirs (e.g. ../bundle_msys2/bin) that this layout does not
# have; the DLL closure is copied from /ucrt64/bin below instead.
set -e
SRC="${ELMER_SRC:-$(cd "$(dirname "$0")" && pwd)}"
BUILD="${ELMER_BUILD:-$SRC/../elmer-build-win}"
INSTALL="${ELMER_INSTALL:-$SRC/../elmer-install-win}"

cd "$BUILD"
cmake -DCPACK_BUNDLE_EXTRA_WINDOWS_DLLS=OFF . > reconfigure.log 2>&1
ninja install > install2.log 2>&1
echo INSTALL_OK

echo "=== key binaries ==="
ls "$INSTALL/bin" | grep -iE "ElmerSolver|ElmerGrid|libelmersolver" || true
echo "=== solver modules of interest ==="
ls "$INSTALL/share/elmersolver/lib" | grep -iE "MagnetoDynamics|CircuitsAndDynamics|ProcessFields|LoadFields|SaveData" || true
echo "=== transhom keywords in installed SOLVER.KEYWORDS ==="
grep -c "Transient Homogenization\|Nu 11 y0\|Homogenization Ladder Order" "$INSTALL/share/elmersolver/lib/SOLVER.KEYWORDS"

echo "=== copying UCRT64 runtime DLL closure into bin ==="
for target in "$INSTALL"/bin/*.exe "$INSTALL"/bin/*.dll "$INSTALL"/share/elmersolver/lib/*.dll; do
  ldd "$target" 2>/dev/null | grep -io "/ucrt64/bin/[a-z0-9_.+-]*\.dll" || true
done | sort -u > /tmp/elmer_dll_deps.txt
while read dll; do
  cp -n "$dll" "$INSTALL/bin/" && echo "copied $(basename $dll)" || true
done < /tmp/elmer_dll_deps.txt

echo "=== self-check: solver banner ==="
cd /tmp
PATH="$INSTALL/bin:$PATH" ELMER_HOME="$INSTALL" "$INSTALL/bin/ElmerSolver.exe" 2>&1 | head -14 || true
