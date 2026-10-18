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
- Magnetic energy of homogenized litz windings (`CalcFields`, DEV-1590). `res: Magnetic Field
  Energy` (with `Separate Magnetic Energy`) now holds what a stranded coil with `Homogenization
  Model` stores, so that a harmonic 2 W/(VF |I|^2), I the peak current, is the circuit's inductance.
  Harmonic: 1/2 Re(nu) |B|^2 of the homogenized tensor, rotated at the integration point as
  WhitneyAVHarmonicSolver assembles it, and 1/2 Im(1/sigma 33)/omega |J|^2 for the reactance of a
  complex `Sigma 33` (none at 0 Hz or with an explicit `Resistance`); with or without `Calculate
  Homogenization Loss`, which before decided whether the coil counted at all (12.7 % low in the test
  with it, the material reluctivity without it). Transient: the Foster ladder storage y0 B^2/2 + sum
  r xi (B - xi/2) in the local Alpha/Beta frame from the published `Xi Alpha`/`Xi Beta` states (3D),
  and the skin inductor's L chi^2/2, chi = I - tau q; before, the material reluctivity (6 % low).
  The skin inductor energy takes the circuit's weight N_j^2 |w|^2, as K does, instead of the
  normalized winding vector (0.125 % off the circuit's reactive power in the test, 2.5e-7 now); the
  skin heat keeps the normalized vector, like the DC and the harmonic loss it converges to. These
  blocks report coenergy = energy. Harmonic energies are peak-phasor values, twice the time average.
  Foil sheet blocks: the harmonic skin reactance of the strands, 1/2 Im(1/sigma_s)/omega |J|^2, counts
  too (1.5 % of the inductance in `circuits_harmonic_foilsheet_energy`). Transient (DEV-1594): the
  in-plane reluctivity ladder stores as the litz one does, from the `Foil Sheet Xi` states of the two
  in-plane directions, and strand j's skin ladder stores w_j (sum_k tau_k x_kj^2 + L_tail y_j^2/2), w_j
  its DC resistance weight, which `CircuitsAndDynamics` publishes every step as `Foil Sheet Skin
  Energy Quadratic`/`Linear`/`Constant`; before, the DC reluctivity alone (9.3 % low with both
  ladders, 4.6 % with the skin ladder alone). The ladder order limit is `MaxHomogLadderOrder`
  (`CircuitUtils`). Tests: `circuits_harmonic_stranded_homogenization_energy` (W against the coil's
  reactive power, complex `Sigma 33`, no loss output), `circuits_transient_stranded_homogenization_energy`
  (mean W over the last period against the BDF2-corrected reactive power),
  `circuits_transient_foilsheet_energy` (the same with both foil sheet ladders, the skin ladder alone and
  both ladders stacked along Beta, on 1 and 4 ranks, within 1e-4: the currents circulating between the
  sheet segments decay slowly, so each case runs 3 or 5 periods).
- Windows SIF path handling fix in the Lua layer.
- `CircuitsAndDynamics` — nonlinear (solution-dependent) lumped circuit elements in
  transient runs: the exported `crt i`/`crt v` variables are refreshed from the latest
  coupled iterate on every execution instead of once per timestep, so component keywords
  such as `Resistance = Variable "crt i 2"` converge as a Picard fixed point within each
  timestep (needs `Export Circuit Variables = True` and `Steady State Max Iterations > 1`;
  optional damping via `Circuit Variable Relaxation Factor`).

  Built on the same machinery, `CircuitsAndDynamics` also gains a first-class **diode**
  circuit element, `Component Type = String Diode`, for transient runs. It is a lumped
  element (no `Master Bodies`) stamped on its own voltage row as `R(V)*I - V = 0` with a
  smoothed, voltage-controlled resistance

  ```
  R(V) = Ron + 0.5*(Roff - Ron)*(1 - tanh((V - Vf)/Vs))
  ```

  where `V` is the component's **own** branch voltage, taken internally from the relaxed
  latest coupled iterate. Forward bias is `V > Vf`; the sign convention is handled by the
  component. Because the voltage is read internally, a diode needs **no**
  `Export Circuit Variables` and no MATC expression.

  | Component keyword | Symbol | Default | Notes |
  | --- | --- | --- | --- |
  | `Diode Forward Voltage = Real` | Vf | `0.0` | Threshold; forward bias is `V > Vf`. |
  | `Diode On Resistance = Real` | Ron | *(required)* | Fatal if absent. |
  | `Diode Off Resistance = Real` | Roff | `1.0e5` | Keep <= ~1e6 to protect conditioning. |
  | `Diode Smoothing Voltage = Real` | Vs | `0.1` | Must be > 0. |

  Solve control is Picard only — `R` is frozen at the previous iterate, there is no `dR/dV`
  Newton term — so a stiff diode needs `Steady State Max Iterations > 1` and typically
  `Circuit Variable Relaxation Factor < 1`. The validated recipe for
  `Ron/Roff/Vs = 0.1/1e5/0.1` is a relaxation factor of `0.10` with
  `Steady State Max Iterations = 100` and `Steady State Min Iterations = 3` (both are
  **Simulation**-section keywords), a direct linear solve, and
  `Steady State Convergence Absolute = True`; that converges every timestep of the
  reference rectifier case in a median of 7 (max 30) coupled iterations. Undamped Picard
  diverges, and the measured stability limit is a relaxation factor between 0.15 and 0.20.
  Diodes are **transient-only**: `CircuitsAndDynamicsHarmonic` aborts with a Fatal if one
  is present. Tests: `fem/tests/circuits2D_transient_diode_component_{rectifier,blocking}`
  (component) and `circuits2D_transient_diode_{rectifier,blocking}` (the equivalent model
  built from `Component Type = Resistor` plus a MATC expression, which covers the generic
  Picard mechanism).
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
- `Wsolve` takes `Without RotM` per element from its own component (DEV-1559). Before, the flag
  that a stranded coil with `Without RotM = True` set stayed on for the foil winding, flat wire and
  foil sheet elements visited after it, so their W was solved with the foil conductivity tensor
  left in the global axes instead of the foil's frame. TRAFOLO writes `Without RotM` for every
  non-homogenized stranded coil, so models mixing such a coil with foil windings could get foil
  losses up to ~30 % too high, depending on element order and partitioning. Test:
  `circuits_harmonic_foil_stranded_norotm`.
- The Lagrange multiplier separated from the collection solution (`SolveWithLinearRestriction`)
  follows `Nonlinear System Relaxation Factor` (and `Nonlinear System Relaxation After`) unless
  `Lagrange Multiplier Relaxation Factor` is given, so the exported circuit voltages and currents
  are the same iterate as the relaxed field. Before, `MagnetoDynamicsCalcFields` combined a relaxed
  `A` with an unrelaxed `V` in `E = -i*omega*A - V*grad W` of massive coils, and the current density
  and Joule loss of a nonlinear loop that stopped before convergence came out orders of magnitude
  too large. Test: `circuits_harmonic_massive_relaxation`.
- `residual` convergence measure of a constrained system (DEV-1567). With `Nonlinear System
  Consistent Norm` or `Nonlinear System Convergence Without Constraints`,
  `SolveWithLinearRestriction` hands `ComputeChange` the field rows only, so `Nonlinear System
  Convergence Measure = residual` left out the constraint forces. Where the constraints (a
  circuit) drive the field, the right-hand side of those rows is zero: the first iterate had a
  zero residual and the nonlinear loop stopped after one iteration. The residual is now that of
  the whole constrained system for the previous field and the new multipliers
  (`LaggedResidual` in `SolverUtils`), the value `ComputeChange` gets without those options.
  Unchanged with `Eliminate Linear Constraints` and `Restriction System Residual Mode`. Test:
  `circuits_harmonic_massive_residual`.
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
- `Mumps Drop Zeros = True` (complex systems, `ZMumps_SolveSystem`, DEV-1564) leaves the
  off-diagonal entries that are exactly zero out of the matrix given to MUMPS and logs how many.
  The couplings of Dirichlet dofs are zeroed but stay in the sparsity pattern, where they enter
  the ordering and the fill of the factors. The equations do not change; in the tests the
  results agree with the solve without the option. In the A-V system the edge-potential
  couplings of non-conducting regions are zero as well: in `mgdyn_harmonic_wire_mumps_drop_zeros`
  a third of the entries are left out, the factors have 21 to 27 % fewer entries and their
  computation 36 % fewer flops. The effect depends on the ordering MUMPS picks: on small systems
  ordered by AMF the factors can grow. Tests: `mgdyn_harmonic_wire_mumps_drop_zeros` (with
  `Mumps Diagonal Regularization`), `VectorHelmholtzWaveguideMumpsDropZeros` (PEC walls), both on
  1 and 2 ranks.
- `Mumps Ordering` (`amd`, `amf`, `scotch`, `pord`, `metis`, `qamd` or `auto`; `DirectSolve`, real
  and complex, DEV-1564) selects the fill-reducing ordering of the MUMPS analysis (`ICNTL(7)`).
  Unset keeps MUMPS's automatic choice, which on systems with quasi-dense rows such as circuit
  equations can pick QAMD and much larger factors: in `circuits_harmonic_foilsheet_mumps_ordering`
  the automatic choice is QAMD with 4.96e6 entries in the factors and 2.6e9 flops, METIS gives
  3.20e6 and 9.0e8, and both give the reference solution. An unknown value stops the run. An
  ordering the MUMPS build lacks (the shipped MUMPS 4.10.0 has no SCOTCH and no PORD) falls back
  to MUMPS's own choice; the ordering used (`INFOG(7)`) is logged. Tests:
  `circuits_harmonic_foilsheet_mumps_ordering` (complex), `WinkelBmPoissonMumpsOrdering` (real).
- Failures a caller can see (DEV-1560). A missing solver procedure and a MATC error in the SIF
  exit with status 1 instead of 0 (`Load.c`), and the MATC message no longer prints the
  unterminated Fortran buffer. `Fatal` aborts all MPI ranks (`MPI_Abort`, logged as
  `StopOnError: rank <r> aborts all <n> MPI ranks`) when more than one rank runs, so an error
  raised on some ranks only cannot leave the others waiting. A solve that is kept although it
  did not converge writes one warning in a fixed form at the default output level:
  `WARNING:: <caller>: NOT CONVERGED: linear|nonlinear|coupled ... key=value ...`
  (linear: `solver`, `iterations`, `tolerance`, also when `Global Abort Not Converged = False`
  turns the numerical error into a warning; nonlinear: `solver`, `iterations`, `change`,
  `tolerance`; coupled: `iterations`, `max`; the MUMPS residual check above writes `direct`).
  No warning comes from a failed `Linear System Trialing` attempt that the next strategy
  replaces, from a single nonlinear iteration (`Nonlinear System Max Iterations` 1 or unset,
  a plain linear solve) or from a fixed number of coupled passes (`Steady State Min Iterations`
  equal to the maximum, as in the homogenization scans, whose passes solve different problems).
  A relaxed single nonlinear iteration warns once per solver with `iterations=1 relaxation=...`,
  because its result is the solution blended with the initial guess. Not covered: the outer
  iterations of a steady simulation (`Steady State Max Iterations` > 1 there is a timestep
  count) and the inner iterative solves of a preconditioner, which still warn. The NaN check of
  the solution norm also catches Inf, and its message and the abort messages of the linear and
  nonlinear solvers name the solver's `Equation`. UMFPACK errors are `Fatal` with the solver and
  the status instead of a bare `PRINT`/`STOP`, and a singular matrix (status 1) warns.
  Tests (label `robustness`, macro `RUN_ELMER_EXPECT` in `test_macros.cmake`):
  `robustness_procedure_not_found`, `robustness_matc_syntax_error`,
  `robustness_linear_not_converged`, `robustness_nonlinear_not_converged`,
  `robustness_fixed_coupled_passes`, `robustness_fatal_on_some_ranks`.
- Gauge-free nonlinear convergence measure (DEV-1564). `Nonlinear System Convergence Measure =
  "flux density"` measures the L2 change of B = curl(x) between nonlinear iterations of an edge
  element solver, real or complex, relative to the L2 norm of B (absolute with `Nonlinear System
  Convergence Absolute`) (`CurlChange` in `SolverUtils`). The default `norm` measures the
  potential, whose gradient part a direct solve of the ungauged system (MUMPS, tree gauge off)
  changes on every solve without changing B, so `norm` may stay above the tolerance and run to
  `Nonlinear System Max Iterations` although the field has converged. Not available with
  `Nonlinear System Linesearch` or `Nonlinear System Compute Change in Scaled System` (both stop
  with an error). Tests: `mgdyn_bh_flux_density` (static, iterative),
  `mgdyn_harmonic_flux_density` (harmonic, linear: stops at iteration 2).
- Thread-safe `GaussPointsAdapt` (`SolverUtils`, DEV-1595), the function as upstream ElmerCSC fixed
  it in `f5374b755` (July 2026), unchanged. Its rule cache was shared by the OpenMP threads: threads
  entering the first threaded assembly of a solver together re-ran its initialization, which
  briefly reset the cached p-element flag, so an element could get the plain Gauss rule. With
  `Element = p:1 b:4` that gave a singular bubble block (`LUDecomp: Matrix is singular`) in about
  6 % of `Step_stokes_heat_vec` runs at 16 or 32 threads, none in 600 runs after the fix. The flag
  and the point counts are now set per call, and the cache is initialized once, under a lock.
  TRAFOLO's thermal runs use `HeatSolve`, which assembles serially and does not call it;
  `StatElecSolve` calls it in a threaded loop, but with first-order elements and no integration
  keywords the race only rewrote equal values. The distributed bundle is built without OpenMP,
  so shipped TRAFOLO never ran these loops threaded; the failure showed in threaded test builds.
- MUMPS with a gauge-free A-V system (DEV-1564). Without a gauge the A-V system is singular.
  `Mumps Null Pivot Detection` can then pin a physical pivot that falls below its threshold,
  which gives a wrong solution without an error; without it MUMPS has to delay the singular
  pivots, which can cost much time and memory or stop with `INFOG(1) = -17`, and in a Newton
  loop the gauge part of the solution can grow from iteration to iteration.
  `Mumps Diagonal Regularization = alpha` adds alpha*|a_ii| to the diagonal of the matrix
  MUMPS factorizes (real and complex, distributed: each rank shifts its partial entry), and
  `Mumps Iterative Refinement` (on with a shift) refines against the unshifted matrix, x <- x +
  F^-1 (b - A x): at least two corrections, then on while the componentwise backward error,
  measured with the row scales of the first solution, halves (at most `Mumps Refinement Max
  Iterations`, 10). It must end below `Mumps Refinement Tolerance` (1e-11; in the comparison
  tests below it ends between 2e-16 and 2e-14): a correction that raises the error tenfold, or
  a later one that raises it above the tolerance, is taken back, and a stop above the tolerance
  warns. An inconsistent system (e.g. a coil current that is not divergence free) cannot be
  refined below the level of its inconsistency. A row with next to no activity, below 1e-12 of
  the largest one of its block, is measured against that largest activity (a block-relative
  fallback after the normwise part of Arioli, Demmel and Duff's sparse backward error): its
  residual can be rounding noise in a part of the solution that is zero for practical purposes,
  e.g. a circuit branch without current, which measured against the old floor of 1e-12 of the
  block could exceed the tolerance and warn about a correct solution. Such a row is accepted up
  to the tolerance times its block's largest activity. Test: `mumps_refinement_row_scale`.
  `Mumps Regularize Nodal Dofs = False` leaves the rows of the nodal dofs unshifted. The
  harmonic A-V solver sets it with its electrodynamics or Darwin model: there the potential
  carries the capacitive field, which a shift of its rows next to a conductor would disturb
  (the conductive part of such a row is many orders of magnitude larger), so the potential must
  be fixed somewhere (a boundary value, an electrode); without them the potential exists in
  conductors only and the shift of its rows also fixes the potential constant of a floating
  conductor. In the tests the shifted solutions equal the unshifted direct solve or the
  iterative solver converged to 1e-12 or below, in 2 to 5 corrections. Tests:
  `mgdyn_harmonic_wire_mumps_regularization` (and `_ed`, the edge-only shift),
  `circuits_harmonic_massive_mumps`, `WinkelBmPoissonMumpsRegularized` (real systems),
  `robustness_mumps_refinement` (the warnings).
- Harmonic core loss independent of the source phase (DEV-1599). `MagnetoDynamicsCalcFields`
  evaluated the `Harmonic Loss` terms on Re B and Im B separately and summed them, so for a field
  exponent other than 2 the loss changed when every source was turned by the same angle (unit cube,
  45 degrees: +14 % at exponent 1.6, -16 % at 2.5). The terms are now evaluated on the major and
  minor semi-axes of the elliptical B, which do not depend on the phase; a unidirectional field
  gives the Steinmetz law at its peak B, exponent 2 and real fields (transient runs) give the
  previous values exactly. `ComponentLoss` and its log lines refer to the two axes. Upstream Elmer
  (devel) and `FourierLoss` (for each harmonic's cos/sin pair) have the same per-mode sum. Test:
  `mgdyn_harmonic_loss_phase`.
- `layered` laminate stack model of the harmonic A-V solver (DEV-1600), for stacks of insulated
  laminations whose thickness d exceeds the skin depth. A body with `Laminate Stack Model =
  "layered"`, `Laminate Thickness` d, `Laminate Stack Conductivity` sigma, `Laminate Stacking
  Factor` F (new, default 1) and `Laminate Normal(3)` (new, optional) gets the reluctivity tensor
  nu_par = 1/(F mu g + (1-F) mu0) in the lamination plane and nu_perp = F nu + (1-F)/mu0 along the
  normal, g = tanh(x)/x, x = sqrt(i omega mu sigma) d/2 (`LaminateLayeredReluctivity` in `Utils`);
  without a normal (corner bodies, where the laminations bend) nu_par in every direction. The body
  must have zero conductivity: the net current of every lamination flows in a `Thin Sheet` boundary
  on its mid-surface with the new keyword `Thin Sheet Laminate Thickness` d (and `Thin Sheet
  Electric Conductivity`, `Thin Sheet Relative Permeability` of the lamination), which gives the
  sheet the complex conductance sigma d sinh(kd)/(kd) and no magnetic term (`LaminateSheetFactor`):
  the body already holds i omega mu g d/4 of the net-current impedance, and with this factor the
  impedance between the lamination faces is the exact x coth(x)/(sigma d). `CalcFields` uses the
  same tensor for H and the energy, adds the sheets' reactive part to the field energy and
  publishes `res: Laminate stack loss`: the Im(nu) loss of the bodies plus the sheet loss; only the
  sum is meaningful, the sheet part alone can be negative (d of about 3-5 skin depths). The loss is
  not in the `Joule Heating` field (no thermal coupling yet). An H-B curve, a tensor or complex
  permeability, a nonzero conductivity or a coil component in a `layered` body, and `layered` in
  the transient solver, stop with a Fatal. Validated against a resolved model of three 8 um NiFe
  laminations (mu_r 190, 1.67 MS/m) around a copper bar, 1 MHz - 1 GHz, single bar and a pair 20 um
  apart in and out of phase: R -3.2 .. +0.9 %, L -0.3 .. +2.5 % (the tensor alone: L up to +9.8 %).
  Test: `mgdyn_laminate_layered_harmonic` (MUMPS).
- `res: Thin sheet current power` (`CalcFields`) evaluates the sheet field with the DOF counts and
  the edge basis of the potential solver, on the same mesh face as the assembly (DEV-1600). Before,
  it took its own solver's counts and a Piola basis and reported about 1e9 W instead of 1e-2 W in a
  3-D harmonic model with the default edge elements.

## Building the distributed bundle

The TRAFOLO product ships a **Windows-only** solver zip built with the standard Elmer CMake
system (MSYS2 UCRT64 toolchain). The direct solvers are Elmer's vendored UMFPACK 4.4 and the
public-domain MUMPS 4.10.0 built by `mumps410/build_mumps410.sh` (see `mumps410/README.md`). Every script
that turns this source into that zip is in this repo, and TRAFOLO's release job runs exactly these
scripts (`bash build_all.sh mumps package 0 0`), so following this section reproduces the
distributed archive.

### One-time setup

1. Install [MSYS2](https://www.msys2.org) to `C:\msys64` and the
   [Microsoft MPI](https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi)
   runtime (the SDK is not needed).
2. In an MSYS2 shell:

   ```bash
   pacman -S --needed git make \
     mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-gcc-fortran \
     mingw-w64-ucrt-x86_64-cmake mingw-w64-ucrt-x86_64-ninja \
     mingw-w64-ucrt-x86_64-openblas mingw-w64-ucrt-x86_64-msmpi \
     mingw-w64-ucrt-x86_64-metis mingw-w64-ucrt-x86_64-scalapack
   ```

   Do not install `mingw-w64-ucrt-x86_64-mumps`: that is MUMPS 5.x, CeCILL-C, not GPL-compatible, and is
   never part of the distributed build (see `build_msys2.sh`).
3. Put the MS-MPI 10.1.1 redistributable installer next to this checkout, at
   `../msmpi_redist/msmpisetup.exe`
   ([download](https://download.microsoft.com/download/7/2/7/72731ebb-b63c-4170-ade7-836966263a8f/msmpisetup.exe)).
   The build installs it into the bundle as `redist\msmpisetup.exe`, and packaging refuses to
   run without it.

### Build

From any shell (the driver re-runs itself under MSYS2 UCRT64):

```bash
bash build_all.sh                  # all seven stages, ~60 min
```

Everything is created next to this checkout: `../elmer-build-win` (build tree, wiped on each
build), `../elmer-install-win` (install tree, pruned in place), `../elmer-gates` (validation
scratch) and the zip `../ElmerFEM-nogui-mpi-Windows-AMD64.zip`. Override with `ELMER_BUILD`,
`ELMER_INSTALL`, `ELMER_GATES`, `MUMPS410_WORK`, `MUMPS410_INSTALL`, `CODE_DIR` (zip location), `BUNDLE_NAME` (zip file name) and
`FOLDER_NAME` (top folder inside the zip, default `ElmerFEM-nogui-mpi-Windows-AMD64`). The release
job sets `BUNDLE_NAME=ElmerFEM-nogui-mpi-Windows-trafolo-<date>`; the folder inside is the same.

Done when the driver prints `ALL STAGES OK`, the audit printed `LICENSE AUDIT: CLEAN`, and every
gate in the summary is `PASS`.

| Stage | Script | What it does |
|---|---|---|
| `mumps` | `mumps410/build_mumps410.sh` | Builds the public-domain MUMPS 4.10.0 (METIS orderings, no PORD) into `MUMPS410_INSTALL` |
| `build` | `build_msys2.sh` + `ninja` | Configures with the distributed flag set (UMFPACK on, MUMPS 4.10.0 linked statically, MS-MPI, OpenBLAS, Lua, OpenMP off), builds |
| `deploy` | `deploy_msys2.sh` | `ninja install` + copies the UCRT64 runtime DLL closure into `bin\` |
| `prune` | `prune_install.sh` | Keeps the 18 solver modules the TRAFOLO app uses, drops dev tools, headers and import libs |
| `audit` | `license_audit.sh` | Fails on any GPL-incompatible or unlisted DLL in the import closure, and on an Elmer binary that imports `libgomp-1.dll`; writes `licenses\` and `SOURCE.txt` |
| `gates` | `run_gates.sh` | Runs the validation cases against the pruned install (inputs in `trafolo_bundle/validation/`); all must pass |
| `package` | `package_bundle.sh` | Refuses an unpruned or unaudited install, zips it |

The bundle's Elmer is built without OpenMP (`-DWITH_OpenMP=OFF`, DEV-1596). TRAFOLO runs Elmer
with MPI ranks only (`mpiexec -n <ranks>`). Threaded, some assemblies have races (DEV-1595), and on
Windows libgfortran's `setlocale` can crash when several threads write formatted output.
`libgomp-1.dll` still ships because OpenBLAS imports it. Local test builds use the same flag, so
that ctest runs what ships: unless `OMP_NUM_THREADS` is set, `RUN_ELMER_TEST` sets it to the
physical cores divided by the ranks, so an OpenMP build runs most tests with several threads.

Validation gates (`run_gates.sh`, against the pruned install):

| Gate | What |
|---|---|
| 1a/1b | Impedance BC (3D Whitney harmonic), serial + 4-rank MPI (physical norm; AV norm is gauge-dependent) |
| 2 | Transient homogenization, serial |
| 3a/3b | 2D transient circuits: UMFPACK direct serial + iterative 4-rank |
| 4a/4b | 2D harmonic circuits: UMFPACK direct serial + iterative 4-rank with UMFPACK-backed Circuit preconditioner |
| 5 | ProcessFields attached to gate 1 must not crash (module/core ABI check) |
| 6 | htc_udf smoke (`trafolo_bundle/validation/htc_smoke.sif`): h in [1.5, 15] W/m²K on a 0.1 m vertical plate |
| 7 | No `Caught LUA error` / readsif warnings in any log, incl. an absolute backslash-path run |
| 8 | StatElecSolveVec thin-layer BC responds to permittivity (Eps0 fix regression) |
| 9 | Nonlinear lumped circuit elements: `circuits2D_transient_nonlinear_resistor` + `_picard` |
| 10 | MUMPS 4.10.0 with the TRAFOLO A-V block on `circuits_harmonic_massive`, serial + 4-rank, against the iterative run |

The parallel direct solver is MUMPS 4.10.0 (public domain). MUMPS 5.x (CeCILL-C) and CPardiso
cannot be shipped.

`ninja install` reinstalls every module, so the stages after `deploy` must always run after it.
Stages take positional arguments `FROM STOP_AFTER PROD DRYRUN`:

```bash
bash build_all.sh gates              # resume from a stage after a failure
bash build_all.sh build gates        # stop after a stage (validate without packaging)
bash build_all.sh build package 0 1  # print the stage plan, run nothing
```

`package_bundle.sh` also copies the bundle into the TRAFOLO app's folders under
`%LOCALAPPDATA%` (`Trafolo_third_party_dev`, and `Trafolo_third_party` with `PROD=1`) and stages
`../choke_externals_staging`. Those are conveniences for TRAFOLO developers and do not change the
zip.

The zip matches the distributed one in content, not byte for byte: the build date in the solver
banner, `SOURCE.txt` and file timestamps differ.

Corresponding source for any distributed TRAFOLO Elmer binary is this repository at the commit
recorded in the bundle's `SOURCE.txt` (GPL compliance).
