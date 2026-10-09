# TRAFOLO fork of ElmerFEM

This is [TRAFOLO](https://trafolo.eu)'s fork of [ElmerFEM](https://github.com/ElmerCSC/elmerfem),
maintained as the solver source for the TRAFOLO magnetics-simulation application. It is licensed
under **GPL-2.0+/LGPL-2.1**, identical to upstream Elmer (see `LICENSE.md`, `license_texts/`).

## What differs from upstream

TRAFOLO-authored additions/fixes on the `trafolo` branch (all GPL-2.0+, in
`fem/src/modules/` and `fem/src/`):

- `htc_udf` — temperature-dependent convective heat-transfer-coefficient BC.
- `ProcessFields`, `LoadFields` — transient winding-homogenization post-processing.
  `ProcessFields` scales the core loss (harmonic and transient GSE/iGSE) per node by the optional
  Material keyword `Harmonic Loss Temperature Factor` (e.g. `Variable Temperature` LUA), evaluated
  at the current Temperature field. The Steinmetz thermal correction can't go into the
  `Harmonic Loss * Coefficient` keywords: CalcFields reads them with `ListGetFun(..., Freq)`, which
  evaluates any dependency at the frequency.
- `LoadFields` logs the harmonic result files as a count and first / last index; the full list
  overran the message line at 64 frequencies ("End of record", DEV-1157). Test:
  `loadfields_many_harmonics`. Test `circuits_transient_stranded_foster_ladder` pins the stranded
  reluctivity Foster ladder of order 2 and 4 (equal ladders agree, BDF2, trapezoid edges, start-up,
  two coils of different order in one model).
- `StatElecSolveVec` — thin-layer Robin coefficient fix (missing `Eps0`).
- Transient winding homogenization. The skin ladder of a stranded coil (`Sigma 33 y0`, `Sigma 33
  alpha`, `Sigma 33 Sigma(1,1)`, n = 1) runs in impedance form, R0 + L s/(1 + s tau) with R0 =
  1/G0, L = sigma alpha/G0^2, tau = sigma y0/G0, with one state per component
  (`CircuitsAndDynamics`, DEV-1548, DEV-1549): the diagonal keeps the ladder's conductivity at the
  step, `G_skin`, and the history of the state goes to the right hand side, with the BDF weights of
  the field solver. A DC current now sees K/G0 at any step, and the transient loss converges to the
  harmonic one as the step shrinks; before, the history was left out and the litz loss grew with
  the number of steps (Hybrid Foil-Litz: +25 % at 24, +50 % at 48 points per period). `y0 = 0`,
  where TRAFOLO's fit sits when it hits its cap, is the exact limit: a resistor K/G0 in series with
  the inductance K sigma/G0 (DEV-1549; DEV-1548 stopped there with a Fatal). `CalcFields` spreads
  the ladder heat K [R0 I^2 + L tau q^2] over the winding like the DC loss, instead of
  |J|^2/G_skin. New circuit scalars `r_dc_skin_component(k)` (K/G0), `p_skin_component(k)` (that
  heat) and `v_skin_component(k)` (the ladder voltage: the coil voltage without its flux linkage);
  `r_component` stays the diagonal K/G_skin. Needs `CircuitsOutput` in every time step. Stops with
  a Fatal on a negative `Sigma 33 y0` or `Sigma 33 alpha`, on an explicit `Resistance` with `Sigma
  33 alpha` /= 0, and on a restart that continues a transient or does not set `Restart Time = 0`
  (the state is not in restart files). The stranded coil's flux coupling takes the variable-step
  BDF weights (`TransientLadderBDF`) instead of a fixed 1.5/-2/0.5. Tests:
  `circuits_transient_stranded_skinladder_dc`, `_sine`, `_rl`, `_energy`, `_fatal`.
- Windows SIF path handling fix in the Lua layer.
- `CircuitsAndDynamics` — nonlinear (solution-dependent) lumped circuit elements in
  transient runs: the exported `crt i`/`crt v` variables are refreshed from the latest
  coupled iterate on every execution instead of once per timestep, so component keywords
  such as `Resistance = Variable "crt i 2"` converge as a Picard fixed point within each
  timestep (needs `Export Circuit Variables = True` and `Steady State Max Iterations > 1`;
  optional damping via `Circuit Variable Relaxation Factor`).
- `Coil Type = flat wire` — rectangular (flat) wire windings as a block of turn cells
  (`CircuitsAndDynamics`, 3D, transient and harmonic). The block is split by the normalized
  `Alpha`/`Beta` direction fields into one cell per turn along `Stacking Direction = beta|alpha`
  (default `beta`, optionally a grid with `Turns Across`); every cell owns a voltage dof and
  carries the full component current, so `V = sum_k V_k` and the per-turn current is exact.
  Unlike the foil polynomial model the FEM resolves skin effect inside thick turns. Same
  solver chain as foil (DirectionSolver Alpha/Beta, RotMSolver, Wsolve). Component keywords:
  `Number of Turns` (integer), `Electrode Boundaries(2)`, `Circuit Equation Voltage Factor`,
  optional `Conductor Thickness` [+ `Conductor Width`] for a geometric fill factor or `Fill
  Factor` directly. Per-turn voltages are exported as `v_component(i) dof k`. Tests:
  `circuits_transient_flatwire`, `circuits_harmonic_flatwire`.
- DC / low-frequency foil and flat wire coils (3D harmonic circuits). The circuit source
  `sigma*V(x)*grad W` is discretely solenoidal only when `V` is constant over the component, so
  for a foil polynomial of order > 0 or a flat wire with more than one cell the ungauged edge-only
  A system at `omega = 0` is singular *and* inconsistent and the linear solver floors on a residual
  (measured 1.4e-4 for the foil test, 6e-7 for flat wire) instead of converging. The fix couples
  the circuit voltage dofs to the nodal scalar potential `v` of the AV solver in both directions
  (`CircuitsAndDynamics`, harmonic and transient kernels; matrix structure in `CircuitUtils`), so
  `v` absorbs the non-solenoidal part and `J = sigma*(-i*omega*a - grad v - V(x) grad W)` is
  divergence free. `CalcFields` adds the matching `-grad v` to `E` for `massive`, `foil winding`
  and `flat wire`. Gated on the component keyword `Activate Constraint = True` (which also turns on
  the AV solver's automatic `v = 0` electrode BC and therefore wants `Electrode Boundaries`); the
  harmonic circuits solver switches it on by itself when the angular frequency is exactly zero, and
  warns if it is explicitly off or if `Electrode Boundaries` is missing. Nothing changes at
  `omega > 0` unless the keyword is set. Tests: `circuits_harmonic_foil_dc`,
  `circuits_harmonic_flatwire_dc`.
- `Coil Type = foil sheet` — foil and flat wire windings as conducting sheets (`CircuitsAndDynamics`,
  3D, harmonic and transient), split into `Sheet Cells` turn cells along the stack and `Sheet
  Segments` strands across the turn width. `Sheet Cells = 0` gives each turn a cell above 0 Hz while
  the cells x segments of all foil sheet components fit `FOIL_SHEET_AUTO_CELL_STRANDS` (500); past
  that the budget left by components with given cells is shared in proportion to the turns, and each
  takes the largest divisor of its turn count that fits its share, but never fewer cells than the stack
  skin depth rule of before (cells at most `FOIL_SHEET_CELL_SKIN_DEPTHS`, 0.55, stack skin depths
  thick); one at DC, one per turn in transient or several-frequency runs (DEV-1554). Lumped turns
  lose the gap fringing loss of the turn next to the gap (10-turn pot core at 50 Hz: R -1.3 %, inner
  turn -20 %), and the budget keeps many-turn windings at the cost of before (Oil-Immersed template,
  6 x 120 turns: one cell per turn took 25 x the time). The count is not tied to the mesh (DEV-1545).
  `Sheet Sublayers = 0`: a winding stacked along Alpha (foil) takes three sub-layers for
  0.5 < t/delta <= 1.5 and one otherwise, on any mesh (DEV-1554); stacked along Beta (edgewise flat
  wire) it keeps its sub-layers on elements up to one skin depth and, while an element spans at most
  2.3 sub-layers, up to two (DEV-1545). The coil body may mix
  linear tetrahedra and linear wedges (prism boundary layers): each wedge is cut into three
  tetrahedra for the strand geometry and its strand source is mapped to the wedge edges, so it stays
  exactly divergence free (DEV-1545; it stopped with a Fatal before). Tests:
  `circuits_harmonic_foilsheet*`, `circuits_transient_foilsheet*`.
- Touching foil sheet and flat wire windings (DEV-1541). `Direction Method = distance` takes the
  pair of boundary values of each body from its own faces, so windings that share a face are
  numbered as a chain (winding k spans k-1..k, any integer shift) and the shared face carries one
  value. A node shared by several bodies gets the value of each: the elemental copy of the field
  (`Alpha Direction`, `Beta Direction`) holds each body's own value, the nodal field that of the
  first body, so windings of a different height or radial build may touch; on a face with a
  boundary value of both bodies the values must agree (Fatal otherwise). The `foil sheet` and
  `flat wire` kernels read that copy (for distance fields only; Laplace fields are read as before)
  and map each component onto its own 0..1 of `Alpha` and `Beta` (published as
  `Foil Sheet Alpha Range` / `Foil Sheet Beta Range`). When coils of different components share
  nodes, `Wsolve` solves W one component at a time, each solve starting from its own boundary
  values, so the electrode values of one winding do not leak into the other and plain CG
  converges. Tests: `circuits_harmonic_foil_distance_touching`,
  `circuits_harmonic_foilsheet_touching*`.
- `WPotentialSolver` (`Wsolve`) grows its rotation matrix buffer to the element it saves; on a
  mesh mixing element types (tetrahedra and wedges in a coil) it wrote past the buffer left by the
  last assembled element and corrupted the heap (DEV-1545).
- `Wsolve` of touching windings (one W solve per component) no longer crashes on a partition whose
  W rows all belong to the component being solved and that holds no electrode face of it, such as
  a partition owning the first component's winding but none of the other one: the Dirichlet
  bookkeeping of the matrix is allocated before the start values are read from it (DEV-1557; the
  solver stopped with a segmentation fault after "Coils touch at ... nodes"). Reached with
  ElmerGrid `-metis` partitions on 4 or more ranks. Test:
  `circuits_harmonic_foilsheet_touching_partition`.
- The Lagrange multiplier separated from the collection solution (`SolveWithLinearRestriction`)
  follows `Nonlinear System Relaxation Factor` (and `Nonlinear System Relaxation After`) unless
  `Lagrange Multiplier Relaxation Factor` is given, so the exported circuit voltages and currents
  are the same iterate as the relaxed field. Before, `MagnetoDynamicsCalcFields` combined a relaxed
  `A` with an unrelaxed `V` in `E = -i*omega*A - V*grad W` of massive coils, and the current density
  and Joule loss of a nonlinear loop that stopped before convergence came out orders of magnitude
  too large. Test: `circuits_harmonic_massive_relaxation`.
- MUMPS direct solver interface (`DirectSolve`): the status of every MUMPS initialization,
  analysis, factorization and solve is checked, so a failure stops all ranks with
  `INFOG(1)`/`INFOG(2)` (`INFO` in the rank-local wrappers) and a hint instead of returning an
  unusable solution. New options `Mumps Null Pivot Detection` (`ICNTL(24)`), `Mumps Null Pivot
  Tolerance` (`CNTL(3)`) and `Mumps Sequential Root` (`ICNTL(13)`); the ordering, factor size and
  null-pivot threshold used are logged. `Mumps Residual Check = True` reports the residual of the
  system handed to MUMPS before and after each solve, separately for the field rows and the
  constraint (circuit) rows. Tested with the public-domain MUMPS 4.10.0 (DEV-1462).
  The residual is also checked after every solve of the distributed double-precision wrappers
  (real and complex) by default (`Mumps Residual Check = False` switches it off): MUMPS can report
  success on a singular system it factorized without null pivot detection, so a normwise residual
  max|r| / (max|A||x| + max|b|) above `Mumps Residual Tolerance` (default 1e-8) writes
  `WARNING:: MumpsResidualCheck: NOT CONVERGED: direct solver=...`. On 156 TRAFOLO A-V runs with
  the tree gauge off, correct solutions stayed at or below 1.7e-10 and failed solves reached 5e-2;
  the componentwise indicator, which reached 0.9 on correct runs, is only logged. The check
  detects failed solves, it does not verify a result: a solution can satisfy the linear system
  and still be wrong, and on rows shared by ranks |A| sums the moduli of the partial entries, so
  both ratios can understate the residual of the assembled system when partial entries cancel.
  The A-V solvers warn when a direct solver switches the tree gauge on by itself in a model with
  circuits (and, harmonic, with the electrodynamics model): on `circuits_harmonic_massive` that
  setup gave 0.11 J field energy instead of 0.6165 J with exit code 0, while `Use Tree Gauge =
  False` with `Mumps Null Pivot Detection = True` (tested with tolerance 1e-13) and `Mumps
  Sequential Root = True` matches the iterative solver on 1 and 4 ranks.

## Building

The TRAFOLO product is a **Windows-only** solver bundle built with the standard Elmer CMake
system (MSYS2 UCRT64 toolchain). The direct solvers are Elmer's vendored UMFPACK 4.4 and, when
`MUMPS_PREFIX` is set, the public-domain MUMPS 4.10.0 built by `mumps410/build_mumps410.sh`
(METIS 5.1 ordering, no PORD; see `mumps410/README.md`). MUMPS 5.x is CeCILL-C and is never
linked into a distributed build. The scripts that compile and install the distributed binaries
are in the repo root; see the header of `build_msys2.sh` for the toolchain packages and path
overrides:

```bash
# MSYS2 UCRT64 shell, Microsoft MPI runtime installed
bash build_bundle_msys2.sh            # all steps below, in this order
bash mumps410/build_mumps410.sh       # MUMPS 4.10.0 into ../mumps410-install
MUMPS_PREFIX=../mumps410-install \
bash build_msys2.sh                   # configure (the distributed flag set; wipes the build dir)
ninja -C ../elmer-build-win           # build
bash deploy_msys2.sh                  # install to ../elmer-install-win + copy runtime DLLs
bash prune_install.sh                 # remove what the bundle does not ship
```

The distributed bundle is this install tree after `prune_install.sh`. There is no CI here —
validation is run locally.

Corresponding source for any distributed TRAFOLO Elmer binary is this repository at the commit
recorded in the bundle's `SOURCE.txt` (GPL compliance).
