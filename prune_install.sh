#!/bin/bash
# Last build step: removes from the deploy_msys2.sh install tree what the TRAFOLO bundle does not
# ship (unused solver modules, developer tools, headers, import libraries). Same env overrides.
set -e
SRC="${ELMER_SRC:-$(cd "$(dirname "$0")" && pwd)}"
INSTALL="${ELMER_INSTALL:-$SRC/../elmer-install-win}"
LIB="$INSTALL/share/elmersolver/lib"
KEEP="CircuitsAndDynamics CoilSolver CoordinateTransform DirectionSolver HeatSolve LoadFields \
MagnetoDynamics MagnetoDynamics2D Poisson ProcessFields ReloadData ResultOutputSolve SaveData \
StatCurrentSolve StatElecSolve StatElecSolveVec WPotentialSolver htc_udf"

for m in $KEEP; do
  [ -f "$LIB/$m.dll" ] || echo "MISSING WHITELIST MODULE: $m.dll"
done
for f in "$LIB"/*.dll; do
  case " $KEEP " in *" $(basename "$f" .dll) "*) ;; *) rm -f "$f" ;; esac
done
rm -f "$INSTALL"/bin/{Mesh2D.exe,ViewFactors.exe,Radiators.exe,GebhardtFactors.exe} \
      "$INSTALL"/bin/{elmerf90,elmerf90.bat,elmerld,elmerld.bat,ElmerSolver_mpi.bat} \
      "$INSTALL"/bin/*.dll.a "$INSTALL"/bin/*.a
rm -rf "$INSTALL/include" "$INSTALL/share/elmersolver/include" "$INSTALL/lib" \
       "$INSTALL/share/elmersolver/m4" "$INSTALL/share/doc"
# libmpi_stubs.dll only when a shipped binary imports it
if ! objdump -x "$INSTALL"/bin/*.exe "$INSTALL"/bin/*.dll "$LIB"/*.dll 2>/dev/null \
     | grep -q "DLL Name: libmpi_stubs.dll"; then
  rm -f "$INSTALL/bin/libmpi_stubs.dll"
fi
