#!/bin/bash
# Phase 3 -- prune the install tree to what TRAFOLO uses (17-module whitelist).
set -e
INSTALL="${ELMER_INSTALL:-$(cd "$(dirname "$0")" && pwd)/../elmer-install-win}"
LIB="$INSTALL/share/elmersolver/lib"

echo "=== size before prune ==="
du -sh "$INSTALL"

# --- solver module whitelist (16 from the SIF-builder grep + htc_udf) ---
KEEP="CircuitsAndDynamics CoilSolver CoordinateTransform DirectionSolver HeatSolve \
LoadFields MagnetoDynamics MagnetoDynamics2D Poisson ProcessFields \
ReloadData ResultOutputSolve SaveData StatCurrentSolve StatElecSolve StatElecSolveVec \
WPotentialSolver htc_udf"

cd "$LIB"
mkdir -p /tmp/elmer_keep
for m in $KEEP; do
  [ -f "$m.dll" ] && mv "$m.dll" /tmp/elmer_keep/ || echo "MISSING WHITELIST MODULE: $m.dll"
done
mv SOLVER.KEYWORDS /tmp/elmer_keep/
# everything left in lib/ is deletable physics modules
DELETED=$(ls *.dll 2>/dev/null | wc -l)
rm -f *.dll
mv /tmp/elmer_keep/* .
echo "deleted $DELETED non-whitelist modules; kept: $(ls *.dll | wc -l) + SOLVER.KEYWORDS"

# --- dev tools and build-time files ---
cd "$INSTALL/bin"
rm -f Mesh2D.exe ViewFactors.exe Radiators.exe GebhardtFactors.exe \
      elmerf90 elmerf90.bat elmerld elmerld.bat ElmerSolver_mpi.bat \
      *.dll.a *.a
# keep: ElmerSolver.exe ElmerSolver_mpi.exe ElmerGrid.exe + runtime DLLs

rm -rf "$INSTALL/include" "$INSTALL/share/elmersolver/include"
rm -rf "$INSTALL/lib"          # cmake export dirs / import libs if present
rm -rf "$INSTALL/share/elmersolver/m4" "$INSTALL/share/doc" 2>/dev/null || true

# --- libmpi_stubs: drop only if no shipped binary imports it ---
cd "$INSTALL"
NEED_STUBS=0
for f in bin/*.exe bin/*.dll share/elmersolver/lib/*.dll; do
  if objdump -x "$f" 2>/dev/null | grep -q "DLL Name: libmpi_stubs.dll"; then
    echo "libmpi_stubs needed by: $f"; NEED_STUBS=1
  fi
done
if [ "$NEED_STUBS" = "0" ]; then rm -f bin/libmpi_stubs.dll && echo "libmpi_stubs.dll dropped (no importers)"; fi

echo "=== closure check: imports not present in bin\\ or system ==="
for f in bin/*.exe bin/*.dll share/elmersolver/lib/*.dll; do
  objdump -x "$f" 2>/dev/null | grep "DLL Name" | awk '{print $3}'
done | sort -u | while read d; do
  case "$d" in
    KERNEL32.dll|msvcrt.dll|ucrtbase.dll|api-ms-*|ADVAPI32.dll|WS2_32.dll|USER32.dll|SHELL32.dll|ole32.dll|RPCRT4.dll|bcrypt.dll|CRYPT32.dll|dbghelp.dll|ntdll.dll|msmpi.dll|MSMPI.DLL) ;; # system / prerequisite
    *) [ -f "bin/$d" ] || echo "MISSING FROM bin\\: $d" ;;
  esac
done
echo "=== size after prune ==="
du -sh "$INSTALL"
