#!/bin/bash
# Phase 5 -- validation gates 1-9 against the (pruned) install.
#
# Direct-solver (UMFPACK) gate design notes:
#  * MUMPS is never shipped (CeCILL-C, GPL-incompatible); the bundled direct
#    solver is Elmer's vendored UMFPACK 4.4. Gates 3/4 exercise it.
#  * 3D Whitney AV (IBC test) is UNGAUGED with Piola basis (tree gauge unavailable)
#    -> a direct factorization of that singular system is ill-posed by formulation,
#    so it is NOT a direct-solver gate. The app uses iterative solvers for 3D AV.
#  * TRAFOLO's production direct-solver use case is the 2D nodal systems
#    (MagnetoDynamics2D + circuits): unique solutions, reference norms hold.
#    real    <- transient case (circuits2D_transient_stranded)
#    complex <- harmonic  case (circuits2D_harmonic_stranded_homogenization)
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="${ELMER_SRC:-$HERE}"
INSTALL="${ELMER_INSTALL:-$SRC/../elmer-install-win}"
VAL="${ELMER_VAL:-$HERE/trafolo_bundle/validation}"
GATES="${ELMER_GATES:-$SRC/../elmer-gates}"
MPIEXEC="/c/Program Files/Microsoft MPI/Bin/mpiexec.exe"

export ELMER_HOME="$INSTALL"
export ELMER_LIB="$INSTALL/share/elmersolver/lib"
export PATH="$INSTALL/bin:$PATH"

rm -rf "$GATES"; mkdir -p "$GATES"
declare -A RESULT

check_passed() { # dir label
  if [ "$(cat "$1/TEST.PASSED" 2>/dev/null)" = "1" ]; then RESULT[$2]=PASS; else RESULT[$2]=FAIL; fi
  echo "[$2] ${RESULT[$2]}"
}

# Parallel leg of the direct-solver gates: Elmer implements parallel direct only
# for MUMPS/CPardiso ("CheckLinearSolverOptions: Only MUMPS and CPardiso...");
# UMFPACK is serial-only, and MUMPS is not shippable. The app emits iterative
# solvers everywhere, so the honest parallel gate is the ITERATIVE 4-rank run:
# it must complete and its PHYSICAL solver norm must match the serial reference
# (the circuit-coupled master norm is a serial-only artifact).
# The harmonic case additionally exercises the UMFPACK-backed Circuit
# preconditioner (VankaCreate CircuitPrec) in parallel.
check_mpi_solver() { # logfile solverN label
  if grep -q "ALL DONE" "$1" && grep -qE "Solver $2 PASSED" "$1" \
     && ! grep -qiE "umfpack.*(error|fail)|Fatal" "$1"; then
    RESULT[$3]=PASS; else RESULT[$3]=FAIL; fi
  echo "[$3] ${RESULT[$3]}"
}

# The ungauged 3D AV norm is solver/partitioning dependent (gauge freedom).
# For MPI runs of the IBC case, pass on completion + the PHYSICAL norm
# (Solver 5 SaveScalars) + no errors -- not on the raw AV norm.
check_physics() { # logfile label
  if grep -q "ALL DONE" "$1" && grep -q "CompareToReferenceSolution: Solver 5 PASSED" "$1" \
     && ! grep -qiE "(MUMPS|UMFPACK|umfpack).*(ERROR|fail)|singular matrix|Fatal" "$1"; then
    RESULT[$2]=PASS; else RESULT[$2]=FAIL; fi
  echo "[$2] ${RESULT[$2]}"
}

# Switch MagnetoDynamics* solver blocks to Direct/UMFPACK, leaving auxiliary
# solvers (direction, coil, w-potential) on their reference iterative setup.
# Notes: "Solver N :: Reference Norm" trailer lines are NOT block headers
# (the :: guard) and the buffer is flushed at EOF -- losing the reference
# trailers silently disables CompareToReferenceSolution.
direct_main_solvers() { # sif-file
  awk '
    function flushbuf(  i, l, hasiter) {
      hasiter = 0
      for (i = 1; i <= n; i++)
        if (buf[i] ~ /^[^!]*Linear System Iterative Method = [A-Za-z0-9]/) hasiter = 1
      for (i = 1; i <= n; i++) {
        l = buf[i]
        if (hit) {
          if (l ~ /^[^!]*Linear System Iterative Method = [A-Za-z0-9]/) {
            sub(/Linear System Iterative Method = [A-Za-z0-9]+/, "Linear System Direct Method = UMFPACK", l)
          } else if (l ~ /^[^!]*Linear System Solver = "?[Ii]terative"?/) {
            sub(/Linear System Solver = "?[Ii]terative"?/, "Linear System Solver = Direct", l)
            print l
            if (!hasiter) print "   Linear System Direct Method = UMFPACK"
            continue
          }
        }
        print l
      }
      n = 0; inblk = 0; hit = 0
    }
    /^Solver [0-9]+/ && $0 !~ /::/ { if (inblk) flushbuf(); inblk = 1 }
    {
      if (inblk) {
        buf[++n] = $0
        if ($0 ~ /Procedure = "MagnetoDynamics/) hit = 1
        if ($0 ~ /^End/) flushbuf()
      } else print
    }
    END { if (inblk) flushbuf() }
  ' "$1" > "$1.tmp" && mv "$1.tmp" "$1"
}

# ---------- gate 1: impedance BC (3D Whitney harmonic, iterative), serial + 4 ranks ----------
D="$GATES/g1_ibc"; mkdir -p "$D"
cp "$SRC/fem/tests/mgdyn_harmonic_wire_impedanceBC"/{IBC.sif,wire.grd,ELMERSOLVER_STARTINFO} "$D/"
( cd "$D" && ElmerGrid.exe 1 2 wire > grid.log 2>&1 && ElmerSolver.exe > serial.log 2>&1 )
check_passed "$D" "1a-IBC-serial"
( cd "$D" && rm -f TEST.PASSED && ElmerGrid.exe 2 2 wire -metis 4 4 -removeunused > grid4.log 2>&1 \
  && "$MPIEXEC" -n 4 ElmerSolver_mpi.exe > mpi.log 2>&1 )
check_physics "$D/mpi.log" "1b-IBC-mpi4"

# ---------- gate 2: transient homogenization, serial (iterative reference) ----------
D="$GATES/g2_transhom"; mkdir -p "$D"
cp -r "$SRC/fem/tests/circuits_transient_stranded_homogenization"/{sif,2241,ELMERSOLVER_STARTINFO} "$D/"
mkdir -p "$D/2241/dat"
( cd "$D" && ElmerSolver.exe > serial.log 2>&1 )
check_passed "$D" "2-transhom-serial"

# ---------- gate 3: 2D transient stranded circuits, UMFPACK serial + iterative 4 ranks ----------
# (NOT the 3D transhom case: its WhitneyAVSolver system is ungauged 3D edge-element
#  AV -- ill-posed for direct factorization by formulation, and its pinned norms are
#  gauge-dependent. The app's direct-solver use case is the 2D nodal systems.)
D="$GATES/g3_trans2d_direct"; mkdir -p "$D"
cp -r "$SRC/fem/tests/circuits2D_transient_stranded"/{sif,1381,ELMERSOLVER_STARTINFO} "$D/"
mkdir -p "$D/1381/dat"
SIF=$(ls "$D"/sif/*.sif | head -1)
direct_main_solvers "$SIF"
( cd "$D" && ElmerSolver.exe > serial.log 2>&1 )
check_passed "$D" "3a-trans2D-UMFPACK-serial"
# parallel leg: iterative (no parallel direct without MUMPS); physical solver = 6
DI="$GATES/g3_iter_mpi"; mkdir -p "$DI"
cp -r "$SRC/fem/tests/circuits2D_transient_stranded"/{sif,1381,ELMERSOLVER_STARTINFO} "$DI/"
mkdir -p "$DI/1381/dat"
( cd "$DI" && ElmerGrid.exe 2 2 1381 -metis 4 4 -removeunused > grid4.log 2>&1 \
  && "$MPIEXEC" -n 4 ElmerSolver_mpi.exe > mpi.log 2>&1 )
check_mpi_solver "$DI/mpi.log" 6 "3b-trans2D-iter-mpi4"

# ---------- gate 4: 2D harmonic homogenization, UMFPACK serial + iterative 4 ranks ----------
D="$GATES/g4_harm2d_direct"; mkdir -p "$D"
cp -r "$SRC/fem/tests/circuits2D_harmonic_stranded_homogenization"/{sif,2077,ELMERSOLVER_STARTINFO} "$D/"
mkdir -p "$D/2077/dat"
SIF=$(ls "$D"/sif/*.sif | head -1)
direct_main_solvers "$SIF"
( cd "$D" && ElmerSolver.exe > serial.log 2>&1 )
check_passed "$D" "4a-harm2D-UMFPACK-serial"
# parallel leg: iterative with UMFPACK-backed Circuit preconditioner; physical solver = 3
DI="$GATES/g4_iter_mpi"; mkdir -p "$DI"
cp -r "$SRC/fem/tests/circuits2D_harmonic_stranded_homogenization"/{sif,2077,ELMERSOLVER_STARTINFO} "$DI/"
mkdir -p "$DI/2077/dat"
( cd "$DI" && ElmerGrid.exe 2 2 2077 -metis 4 4 -removeunused > grid4.log 2>&1 \
  && "$MPIEXEC" -n 4 ElmerSolver_mpi.exe > mpi.log 2>&1 )
check_mpi_solver "$DI/mpi.log" 3 "4b-harm2D-CircuitPrec-mpi4"

# ---------- gate 5: ProcessFields attach (app-faithful solver block) ----------
D="$GATES/g5_procfields"; mkdir -p "$D"
cp "$SRC/fem/tests/mgdyn_harmonic_wire_impedanceBC"/{IBC.sif,wire.grd,ELMERSOLVER_STARTINFO} "$D/"
sed -i 's/Active Solvers(3) = 1 2 3/Active Solvers(4) = 1 2 3 6/' "$D/IBC.sif"
cat >> "$D/IBC.sif" <<'EOF'

Solver 6
  Equation = post process data
  Procedure = "ProcessFields" "ProcessFields"
  Update Exported Variables = Logical True
End
EOF
( cd "$D" && ElmerGrid.exe 1 2 wire > grid.log 2>&1 && ElmerSolver.exe > serial.log 2>&1 )
check_passed "$D" "5-ProcessFields-attach"
grep -q "Core loss will not be computed" "$D/serial.log" && echo "  (expected core-loss notice present)"

# ---------- gate 6: htc_udf smoke ----------
D="$GATES/g6_htc"; mkdir -p "$D"
cp -r "$SRC/fem/tests/heateq/Mesh" "$D/"
cp "$VAL/htc_smoke.sif" "$D/"
echo "htc_smoke.sif" > "$D/ELMERSOLVER_STARTINFO"
( cd "$D" && ElmerSolver.exe > serial.log 2>&1 )
HTC=$(awk '{print $1}' "$D/htc_smoke.dat" 2>/dev/null | tail -1)
TMX=$(awk '{print $2}' "$D/htc_smoke.dat" 2>/dev/null | tail -1)
echo "  htc max=$HTC  T max=$TMX"
# h: vertical 0.1 m plate band; T: peak is conduction-dominated (k=1 W/mK), allow up to 400 C
if awk -v h="$HTC" -v t="$TMX" 'BEGIN{exit !(h>=1.5 && h<=15 && t>25 && t<400)}'; then
  RESULT[6-htc_udf]=PASS; else RESULT[6-htc_udf]=FAIL; fi
echo "[6-htc_udf] ${RESULT[6-htc_udf]}"

# ---------- gate 7: no Lua errors anywhere + backslash-path regression ----------
D="$GATES/g7_lua"; mkdir -p "$D"
cp -r "$SRC/fem/tests/heateq/Mesh" "$D/"
cp "$VAL/htc_smoke.sif" "$D/"
WINPATH="$(cygpath -w "$D/htc_smoke.sif" 2>/dev/null || echo "$D/htc_smoke.sif")"
( cd "$D" && ElmerSolver.exe "$WINPATH" > abs.log 2>&1 )
LUAERR=$(grep -rl "Caught LUA error" "$GATES" 2>/dev/null | wc -l)
READSIF_WARN=$(grep -rl "WARNING: readsif" "$GATES" 2>/dev/null | wc -l)
if [ "$LUAERR" = "0" ] && [ "$READSIF_WARN" = "0" ]; then RESULT[7-lua-clean]=PASS; else RESULT[7-lua-clean]=FAIL; fi
echo "[7-lua-clean] ${RESULT[7-lua-clean]} (LUA errors: $LUAERR, readsif warnings: $READSIF_WARN)"

# ---------- gate 8: StatElecSolveVec thin-layer BC responds to permittivity ----------
# Regression gate for the Eps0 fix (fork f16dc1896). Before the fix the layer Robin
# term was ~1/Eps0 too stiff, so ANY eps collapsed to a rigid Dirichlet and the three
# norms below were identical. Pass = insulating (eps=1e-6) and Dirichlet (eps=1e6)
# limits differ by >10x, both runs complete. (Limitation: 2D scalar check only.)
D="$GATES/g8_layerbc"; mkdir -p "$D"
cp -r "$SRC/fem/tests/heateq/Mesh" "$D/"
cat > "$D/base.sif" <<'EOF'
Header
  Mesh DB "." "Mesh"
End
Simulation
  Coordinate System = Cartesian 2D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
End
Constants
  Permittivity Of Vacuum = 8.8542e-12
End
Body 1
  Equation = 1
  Material = 1
End
Equation 1
  Active Solvers(1) = 1
End
Solver 1
  Equation = "es"
  Variable = Potential
  Procedure = "StatElecSolveVec" "StatElecSolver"
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 1000
  Linear System Convergence Tolerance = 1e-10
  Linear System Preconditioning = ILU0
End
Material 1
  Relative Permittivity = 1.0
End
Boundary Condition 1
  Target Boundaries(1) = 1
  Potential = 0.0
End
Boundary Condition 2
  Target Boundaries(1) = 3
  Layer Relative Permittivity = Real LAYEREPS
  Layer Thickness = Real 1.0e-4
  Electrode Potential = Real 1.0
End
EOF
echo "es.sif" > "$D/ELMERSOLVER_STARTINFO"
N_LO=""; N_HI=""
for EPS in 1.0e-6 1.0e6; do
  sed "s/LAYEREPS/$EPS/" "$D/base.sif" > "$D/es.sif"
  ( cd "$D" && ElmerSolver.exe > "r_$EPS.log" 2>&1 )
  N=$(grep -oE "NRM,RELC\): \( *[0-9.E+-]+" "$D/r_$EPS.log" | tail -1 | grep -oE "[0-9.E+-]+$")
  [ "$EPS" = "1.0e-6" ] && N_LO="$N" || N_HI="$N"
done
echo "  layer-BC norms: eps=1e-6 -> $N_LO | eps=1e6 -> $N_HI"
if awk -v lo="$N_LO" -v hi="$N_HI" 'BEGIN{lo+=0; hi+=0; exit !(lo>0 && hi>0 && hi/lo>10)}'; then
  RESULT[8-StatElecVec-LayerBC]=PASS; else RESULT[8-StatElecVec-LayerBC]=FAIL; fi
echo "[8-StatElecVec-LayerBC] ${RESULT[8-StatElecVec-LayerBC]}"

# ---------- gate 9: nonlinear lumped circuit elements (per-iteration export refresh) ----------
# Both cases assert analytic roots of a solution-dependent Resistance (Picard fixed
# point within a timestep); reference norms in the SIFs. Feature is opt-in
# (Export Circuit Variables + Steady State Max Iterations > 1) - these gate it ON.
for tc in circuits2D_transient_nonlinear_resistor circuits2D_transient_nonlinear_resistor_picard; do
  D="$GATES/g9_${tc##*_}"; mkdir -p "$D"
  cp -r "$SRC/fem/tests/$tc"/{sif,1381,ELMERSOLVER_STARTINFO} "$D/"
  mkdir -p "$D/1381/dat"
  ( cd "$D" && ElmerSolver.exe > serial.log 2>&1 )
  check_passed "$D" "9-${tc##*transient_}"
done

# ---------- gate 10: MUMPS 4.10 with the TRAFOLO A-V block, serial + 4 ranks ----------
# 3D harmonic A-V with a massive circuit coil, solved with the block TRAFOLO writes when the
# user selects MUMPS (tree gauge off, null pivots at 1e-10, sequential root, working space 200),
# against the iterative run of the same case on the same partitioning. The A-V norm of an
# ungauged direct solve is gauge-dependent, so the reference norms of the test cannot be used:
# pass = eddy current power, Joule loss and the source current and voltage within 1e-4, no
# NOT CONVERGED or tree-gauge warning, and the banner reports MUMPS.
mumps_av_block() { # sif: the WhitneyAVHarmonicSolver block gets the TRAFOLO MUMPS lines
  awk '
    function flushbuf(  i, l, k) {
      for (i = 1; i <= n; i++) {
        l = buf[i]; k = tolower(l)
        if (hit && k ~ /^ *(linear system (solver|iterative method|preconditioning|convergence tolerance|max iterations|residual output|abort not converged|robust|direct method)|bicgstabl polynomial degree|use tree gauge) *=/) continue
        print l
        if (hit && k ~ /^ *procedure *= *"magnetodynamics" *"whitneyavharmonicsolver"/) {
          print "   Linear System Solver = Direct"
          print "   Linear System Direct Method = MUMPS"
          print "   Use Tree Gauge = Logical False"
          print "   Mumps Null Pivot Detection = Logical True"
          print "   Mumps Null Pivot Tolerance = Real 1e-10"
          print "   Mumps Sequential Root = Logical True"
          print "   Mumps Percentage Increase Working Space = Integer 200"
        }
      }
      n = 0; inblk = 0; hit = 0
    }
    /^Solver [0-9]+/ && $0 !~ /::/ { if (inblk) flushbuf(); inblk = 1 }
    {
      if (inblk) {
        buf[++n] = $0
        if (tolower($0) ~ /whitneyavharmonicsolver/) hit = 1
        if ($0 ~ /^End/) flushbuf()
      } else print
    }
    END { if (inblk) flushbuf() }
  ' "$1" > "$1.tmp" && mv "$1.tmp" "$1"
}
# The test keeps SaveScalars switched off; switch it on (that block only).
enable_scalars() {
  awk '
    /^Solver [0-9]+/ && $0 !~ /::/ { inblk = 1; n = 0; ss = 0 }
    inblk {
      buf[++n] = $0
      if (tolower($0) ~ /"savescalars"/) ss = 1
      if ($0 ~ /^End/) {
        for (i = 1; i <= n; i++) { l = buf[i]; if (ss) sub(/Exec Solver = Never/, "Exec Solver = Always", l); print l }
        inblk = 0
      }
      next
    }
    { print }
  ' "$1" > "$1.tmp" && mv "$1.tmp" "$1"
}
scalar_col() { # names-file column-name -> 1-based column
  tr -d '\r' < "$1" | grep -i -m1 ": $2 *\$" | sed 's/^ *\([0-9]*\):.*/\1/'
}
value_of() { # names-file dat-file column-name -> value in the last row
  local c; c=$(scalar_col "$1" "$3"); [ -n "$c" ] || return 1
  tail -1 "$2" | awk -v c="$c" '{print $c}'
}
# Powers relative to themselves; source current and voltage as complex phasors, relative to
# the magnitude (a part that is zero to round-off must not fail the gate).
compare_scalars() { # names-file ref.dat run.dat label
  local ok=1 q a b ar ai br bi
  for q in "res: eddy current power" "res: joule loss"; do
    a=$(value_of "$1" "$2" "$q") && b=$(value_of "$1" "$3" "$q") || { echo "  $4: column '$q' missing"; ok=0; continue; }
    awk -v a="$a" -v b="$b" 'BEGIN{d=a-b; if(d<0)d=-d; s=(a<0?-a:a); if(s<1e-30)s=1e-30; exit !(d/s<=1e-4)}' \
      || { echo "  $4: $q differs: iterative $a, MUMPS $b"; ok=0; }
  done
  for q in "res: i_testsource" "res: v_testsource"; do
    ar=$(value_of "$1" "$2" "$q re") && ai=$(value_of "$1" "$2" "$q im") && br=$(value_of "$1" "$3" "$q re") \
      && bi=$(value_of "$1" "$3" "$q im") || { echo "  $4: columns '$q re/im' missing"; ok=0; continue; }
    awk -v ar="$ar" -v ai="$ai" -v br="$br" -v bi="$bi" \
      'BEGIN{d=sqrt((ar-br)^2+(ai-bi)^2); s=sqrt(ar^2+ai^2); if(s<1e-30)s=1e-30; exit !(d/s<=1e-4)}' \
      || { echo "  $4: $q differs: iterative $ar + j $ai, MUMPS $br + j $bi"; ok=0; }
  done
  return $((1 - ok))
}
G10="$SRC/fem/tests/circuits_harmonic_massive"
for leg in iter mumps; do
  for np in 1 4; do
    D="$GATES/g10_${leg}_np$np"; mkdir -p "$D"
    cp -r "$G10"/{sif,1962,ELMERSOLVER_STARTINFO} "$D/"
    mkdir -p "$D/1962/dat"
    enable_scalars "$D/sif/1962.sif"
    [ "$leg" = mumps ] && mumps_av_block "$D/sif/1962.sif"
    if [ "$np" = 1 ]; then
      ( cd "$D" && ElmerSolver.exe > run.log 2>&1 )
    else
      ( cd "$D" && ElmerGrid.exe 2 2 1962 -metis 4 3 -removeunused > grid4.log 2>&1 \
        && "$MPIEXEC" -n 4 ElmerSolver_mpi.exe > run.log 2>&1 )
    fi
  done
done
for np in 1 4; do
  L="g10-MUMPS-np$np"; RESULT[$L]=PASS
  DI="$GATES/g10_iter_np$np"; DM="$GATES/g10_mumps_np$np"
  grep -q "Mumps Null Pivot Detection" "$DM/sif/1962.sif" || { echo "  $L: the MUMPS block was not inserted"; RESULT[$L]=FAIL; }
  grep -q "MUMPS library linked in" "$DM/run.log" || { echo "  $L: the solver reports no MUMPS"; RESULT[$L]=FAIL; }
  grep -q "ALL DONE" "$DM/run.log" && grep -q "ALL DONE" "$DI/run.log" || { echo "  $L: a run did not finish"; RESULT[$L]=FAIL; }
  if grep -q "NOT CONVERGED\|tree gauge on by itself" "$DM/run.log"; then echo "  $L: warning in the MUMPS run"; RESULT[$L]=FAIL; fi
  n=0
  for dat in "$DI"/1962/dat/1962.dat*; do
    case "$dat" in *.names) continue;; esac
    n=$((n + 1))
    compare_scalars "$DI/1962/dat/1962.dat.names" "$dat" "$DM/1962/dat/$(basename "$dat")" "$L" || RESULT[$L]=FAIL
  done
  [ "$n" -gt 0 ] || { echo "  $L: no scalars written"; RESULT[$L]=FAIL; }
  echo "[$L] ${RESULT[$L]}"
done

echo ""
echo "================ GATE SUMMARY ================"
FAILED=0
for k in $(echo "${!RESULT[@]}" | tr ' ' '\n' | sort); do
  printf "  %-30s %s\n" "$k" "${RESULT[$k]}"
  [ "${RESULT[$k]}" = "FAIL" ] && FAILED=1
done
exit $FAILED
