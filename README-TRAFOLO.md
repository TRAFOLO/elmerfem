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
  fixed the same dependence in `CalcFields` and `FourierLoss` on 2026-10-07 (5216c95d2), with one
  amplitude, |B|^2 = |Re B|^2 + |Im B|^2: a circular field gives 2^(n/2) times the alternating loss
  at its peak, where the fork gives 2; the two agree for a unidirectional field, and neither is
  checked against measured rotational loss. The fork's `FourierLoss` keeps the per-mode sum
  (TRAFOLO does not run it). Test: `mgdyn_harmonic_loss_phase` (upstream has a different test of
  the same name).
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
  impedance between the lamination faces is the exact x coth(x)/(sigma d) (floored for thick
  laminations since the DEV-1606 lamination sheets, below). `CalcFields` uses the
  same tensor for H and the energy, adds the sheets' reactive part to the field energy and
  publishes `res: Laminate stack loss`: the Im(nu) loss of the bodies plus the sheet loss; only the
  sum is meaningful, the sheet part alone can be negative (d of about 3-5 skin depths). The loss is
  not in the `Joule Heating` field (no thermal coupling yet). A tensor or complex permeability, a
  nonzero conductivity or a coil component in a `layered` body, and `layered` in the transient
  solver, stop with a Fatal (an H-B curve is accepted since DEV-1604, the transient solver since
  DEV-1606). Validated against a
  resolved model of three 8 um NiFe laminations (mu_r 190, 1.67 MS/m) around a copper bar, 1 MHz -
  1 GHz, single bar and a pair 20 um apart in and out of phase: R -3.2 .. +0.9 %, L -0.3 .. +2.5 %
  (the tensor alone: L up to +9.8 %). Test: `mgdyn_laminate_layered_harmonic` (MUMPS).
- `res: Thin sheet current power` (`CalcFields`) evaluates the sheet field with the DOF counts and
  the edge basis of the potential solver, on the same mesh face as the assembly (DEV-1600). Before,
  it took its own solver's counts and a Piola basis and reported about 1e9 W instead of 1e-2 W in a
  3-D harmonic model with the default edge elements.
- Transient core loss of `ProcessFields` without the switch-on step (DEV-1607). The first step
  starts from the zero initial field, not from the excitation at t = 0; its dB/dt was taken against
  that zero field and the iGSE DeltaB bounds started at B = 0, and without `Averaging Start Time`
  (TRAFOLO with no skipped periods) both went into the averages. A field of 0.5 + 0.1 cos(wt) T over
  one period gave 25 times (GSE) and 109 times (iGSE) the steady core loss. Now the first step has
  no dB/dt, so no core loss, and the bounds start at its field; the Joule averages are unchanged. A
  waveform that does start from zero loses that one step, up to about 2/K of the core loss averaged
  over K steps. `Averaging Start Time` is also read without `Calculate Harmonic Loss`: models
  without core loss ignored it and averaged from t = 0. Tests: `mgdyn_transient_loss_start_gse`,
  `mgdyn_transient_loss_start_igse`; `circuits_transient_foilsheet_skinloss` and
  `circuits_transient_stranded_skinladder_rl` / `_sine` now get their start time without core loss.
- H-B curves in the `layered` laminate stack model (DEV-1604). A `layered` body may take an `H-B
  Curve` instead of `Relative Permeability`: the curve of the lamination steel (as data sheets and
  Epstein tests give it). It is read at the flux density of the laminations, B_s = B_par/F + B_perp
  (the laminations carry the stack's in-plane flux on the fraction F of the cross-section and all of
  the normal flux), and nu = H(|B_s|)/|B_s| is the material reluctivity of the stack tensor, skin
  factor g included (`LaminateMetalFlux`, `LaminateMetalFluxDensity` in `Utils`). The Newton term
  carries dNu/dnu B (relative finite difference, `LaminateLayeredReluctivityDerivative`) and the
  weight B_par/F^2 + B_perp of d|B_s|/dB, with the complex-linear approximation of the isotropic H-B
  term; the converged solution does not depend on it (Newton and Picard agree to 5e-12) and the
  iteration converges like the isotropic term (same iteration history at 1 MHz). `CalcFields` builds
  the same tensor for H and the energy (secant: 1/2 Re(B* Nu B), twice the time average) and
  evaluates the `Harmonic Loss` law of a layered body on B_s over the metal fraction F of its
  volume, the way per-kg loss data of the steel are defined (F^(1-beta) more than on the stack
  average: +3.4 % for F 0.96 and beta 1.8; this also applies to linear layered bodies). Validated
  against resolved laminations (M300-35A, M310-50A at 500 Hz; M270-50A, M235-35A, M165-35S at 50 Hz;
  F 0.96, mean flux density 0.5 - 1.7 T, flux-driven stack): magnetizing field and eddy loss within
  0.5 % at 50 Hz (grain-oriented M165-35S: 1.5 %); at 500 Hz -11 .. +4 % and -8 .. +7 %, largest for
  0.5 mm laminations at 1.0 - 1.3 T, where the skin effect makes the saturation uneven across the
  lamination. The low-frequency model in the same runs: -86 .. -2 % and -6 .. +8 %, mostly because
  it reads the curve at the stack average. A lamination sheet next to a layered body with an H-B
  curve takes the steel's reluctivity there instead of its `Thin Sheet Relative Permeability`: H/B at
  B_s averaged over the neighbouring elements, from the previous iterate, in the solver and in
  `CalcFields` (`LaminateSheetSteelReluctivity`), so sheet and body keep the same permeability, which
  the sheet factor needs. In parallel the sheet must lie inside a partition: ElmerGrid keeps the two
  parents of a boundary face together by default, and a sheet that a partition interface cuts (with
  `-partnobcoptim`, for example) stops with a Fatal, as its parent across the interface is not at
  hand there (every lamination sheet since the DEV-1606 lamination sheets, below). With the fixed
  sheet permeability a saturated 8 um NiFe lamination at
  100 MHz turned active and the iteration did not converge; it now converges, linearly in deep
  saturation since the sheet's dependence on B is not in the Newton term. Maxwell stress has no
  magnetic part in layered bodies (scalar reluctivity 0, as before), and `Relative Permeability`
  output is 1/(mu0 Nu_11) of the stack tensor, not the steel's. Tests:
  `mgdyn_laminate_layered_harmonic_hb_linear` (a linear curve gives the DEV-1600 norm),
  `mgdyn_laminate_layered_harmonic_hb` (saturating curve at 1 MHz; stack loss, loss law and
  magnetic energy; converges at the 12th solve; also on 3 ranks of the default partitioning, with
  the same norms), `mgdyn_laminate_layered_harmonic_hb_split_sheet` (that model on 3 ranks with
  sheets cut by the partition interfaces: the Fatal).
- `layered` laminate stack model in the transient A-V solver (DEV-1606), so that TRAFOLO can move
  harmonic and transient runs off the low-frequency model at once. `WhitneyAVSolver` takes a
  `layered` body (before: a Fatal) with `Relative Permeability` or an `H-B Curve`: the static stack
  tensor of the harmonic model at zero frequency (an H-B curve read at B_s, with the Newton term of
  DEV-1604, real) and, in the lamination plane, the eddy currents of the laminations as a Foster
  ladder per element, H_dyn = (nu(s) - nu(0)) B_par (`LaminateNuFoster` in `CircuitUtils`). The
  ladder holds the lamination modes of tanh(x)/x, the first N-1 of them and one pole that carries the
  rest of their DC weight and first moment, so that the DC value, the low-frequency eddy term and
  the infinite-frequency limit (1-F) mu0 are those of the lamination; the eddy term is
  sigma d^2/(12 F) (F mu_r/(1 - F + F mu_r))^2 dB/dt, sigma d^2/(12 F) dB/dt for F mu_r >> 1 - F
  (exact at F = 1). Its poles are the zeros of the modal permeability (bisection), N = `Laminate
  Ladder Order` (new, 1-6, default 6).
  With the constant DC tail of the foil sheet ladder instead, the eddy loss of a square flux ripple
  (1 % edges) came out up to 8 % low; with the moment-matched pole within 1 % (1-D check, 0.10-0.35
  mm, 1-20 kHz, triangle and square flux; d^2/12 alone: up to +486 %). The ladder follows the steel's
  permeability: `Relative Permeability`, or with a curve the mean incremental reluctivity at the
  laminations' flux density and at that of their core (the slowest ladder state) of the previous
  step, frozen over the step so that Newton stays exact. The incremental permeability is what the
  eddy currents of a ripple on a bias see: the secant, a fixed permeability or d^2/12 miss its loss
  by 20-100 % and more. The states (three components per pole) follow the element mean of B_par,
  advance with the field's own BDF weights and are eliminated per element and step; the solver
  declares them on the layered bodies (`Laminate Xi`, `Laminate Ladder`, about 0.5 kB per element
  at order 6). With F = 1, order 1 and no normal the ladder is sigma d^2/12 dB/dt alone and the model
  is the old `Low-Frequency Model` (same solution to 12 digits), which however reports no eddy loss.
  The eddy loss is the power the ladder takes less the change of its stored energy, the power the
  block absorbs also while the ladder follows the curve (the ladder's own dissipation, sum -r_k T_k
  xi_k'^2, then over-reports it, by up to 30 % at 500 Hz in saturation). It goes into `res: Eddy
  current power` and `res: Laminate stack loss`, and into an elemental `Proximity Loss` when one is
  declared (`ProcessFields` adds that to the core load: a Steinmetz eddy term on top counts the
  classical loss twice). The transient core loss of `ProcessFields` (GSE, iGSE) takes, in a layered
  body, the steel's flux density B_s = B_par/F + B_perp (its dB/dt, its magnitude and the iGSE
  bounds) on the metal fraction F of the volume, as the harmonic loss of `CalcFields` does; the
  exported B and BT stay the stack's. `CalcFields` adds the eddy current field to H and the ladder's stored
  energy to the magnetic energy; the static energy of an H-B body is the potential of its law along
  the path from zero, integrated along the curve (`LaminateLayeredEnergy`). Where the flux has
  in-plane and normal parts the static law has no potential and its Newton term is not symmetric:
  with `Linear System Symmetric` the solver gets its symmetric part, in the matrix and the right
  hand side alike, so that a solver reading one triangle still converges to the law. A Fatal stops
  `Constant (Bulk) System`, `Adaptive Timestepping`, `Convergence Control`, a predictor-corrector
  AV solver, a timestepping method other than BDF of order 1 - 2 or implicit Euler,
  `Electrodynamics Model` (second order in time) and a second WhitneyAVSolver with layered bodies
  (and lamination sheets, until the DEV-1606 lamination sheets below). The states start from the
  initial field as steady, also after a
  restart, whose files carry no ladder history, and in each new simulation of the same process (a
  `Run Control` sweep, a reloaded model, also when the AV solver first runs at a later step of it,
  a solver initialised again). A mesh swap during a run (`Mesh Name Index`, remeshing) initialises
  the solver again, so the states start again from the new mesh's field there.
  Validated on a flux-driven stack in 3-D
  (a quarter of a lamination width, `hfhomog/stack3d.py` transient mode,
  `studies/run_stack_transient.py` in CORE HOMOGENIZATION):
  - linear steel (mu_r 1000 - 5000, 0.35 / 0.50 mm, F 0.96 and 1), sine 50 Hz - 50 kHz, against the
    harmonic layered model, which is exact there: complex permeability and loss within 0.07 % up to
    20 kHz, 0.8 / 0.9 % at 50 kHz (38 - 49 times the first lamination mode); energy balance 1.0005.
  - saturating sine at 50 Hz, 0.5 - 1.7 T, against meshed laminations (transient, same mesh and
    curve; M300-35A, M310-50A, M270-50A, M235-35A, grain-oriented M165-35S; F 0.96): non-oriented
    steels magnetizing field |H1| within 0.2 %, peak H within 0.3 %, eddy loss -0.9 .. +0.1 % (the
    old model: 0.4 %, 0.4 %, -0.8 .. +0.2 %); M165-35S |H1| within 0.5 %, peak H 1.3 %, loss +0.6 /
    -0.9 / -3.9 / -4.2 % at 0.5 / 1.0 / 1.5 / 1.7 T (the old model: |H1| to -1.9 %, loss +1.0 / +0.3
    / -3.2 / -3.5 %). There the skin depth below the knee is about the lamination thickness and the
    surface layers saturate first, so the eddy currents across the lamination are nonlinear, which
    neither a linear ladder nor d^2/12 follows: a 1-D check (resolved lamination against the block)
    gives -3.5 .. -4.4 % for every choice of the ladder's permeability and for d^2/12 alone.
  - 500 Hz, 0.5 - 1.7 T: M300-35A |H1| -1.5 .. +0.9 %, peak H -1.8 .. +1.0 %, eddy loss -5.1 ..
    -2.5 %; M310-50A (0.5 mm) |H1| -4.5 .. +2.2 %, peak H -3.8 .. +16.7 % (0.5 T), loss -9.4 ..
    +0.4 %, the largest misses at 1.0 - 1.7 T, where saturation makes the eddy currents across the
    lamination nonlinear, which a linear ladder does not follow (the 1-D check gives the same for
    every choice of its permeability; the harmonic layered model against its meshed reference: -11
    .. +4 % and -8 .. +7 %). The old model: |H1| down to -11 %, peak H to -17 %, loss -6.8 .. +9.5 %.
  - triangular flux ripple +-0.05 T at 5 kHz on a bias of 0.5 / 1.0 / 1.5 T, M300-35A: ripple |H1|
    -0.5 / +0.2 / +0.2 %, peak H -3.8 / -3.0 / +0.2 %, eddy loss -2.2 / +0.0 / +0.3 %; the old model
    +184 / +55 / +0.9 % loss, +84 / +6 / -0.0 % |H1|.
  - energy balance: reported eddy loss over the absorbed power 1.0005 in the sines; in the ripple
    1.0000 against the absorbed power in the run's own BDF2 quadrature, 1.007 - 1.024 against the
    trapezoidal one (dB/dt jumps at the triangle's corners; the comparisons above use the
    trapezoidal power for both models).
  - against the harmonic layered model at 50 Hz in saturation (3 steels): eddy loss within 0.8 %;
    the field amplitudes differ by method (the harmonic model's H is the curve at the peak flux,
    the transient one the fundamental of the true waveform: -53 .. -79 % at 1.5 - 1.7 T).
  Known limits, accepted: the eddy loss of grain-oriented M165-35S at 50 Hz and 1.5 / 1.7 T comes
  out -3.9 / -4.2 % low, and that of 0.5 mm laminations (M310-50A) at 500 Hz and 1.0 - 1.7 T
  -6.9 .. -9.4 % low. In both the eddy currents across a saturating lamination are nonlinear, which
  a model with a linear ladder on the lamination's mean flux cannot follow, whatever permeability
  the ladder takes (1-D check above); the old model misses them too (-3.2 / -3.5 %, down to
  -6.8 %).
  Tests: `mgdyn_laminate_layered_transient` (linear, 2 kHz, order 6: stack loss, energy with the
  ladder's share, surface field; also on 2 ranks), `mgdyn_laminate_layered_transient_lf` (F = 1,
  order 1: the old model's solution and its eddy loss; also on 4 ranks),
  `mgdyn_laminate_layered_transient_hb` (M300-35A curve, 500 Hz, 1.3 T; also on 3 ranks cut along
  the stack, one of them without a layered element), `mgdyn_laminate_layered_transient_history`
  (a start from a nonzero field, two step sizes, two solves per step that give the result of one,
  the elemental `Proximity Loss`), `mgdyn_laminate_layered_transient_runs` (two Run Control runs
  that give the result of one, tilted laminations with an H-B curve, ending in saturation; with
  the states of the first run kept the loss is 0.02 % off), `mgdyn_laminate_layered_transient_reload`
  (two one-step Run Control runs, the second at 0.8 of the flux: the norms of a single run at 0.8,
  also of an elemental `Proximity Loss`), `mgdyn_laminate_layered_transient_reinit` (a second
  one-step simulation after `RUN` with the AV solver initialised again and 0.8 of the flux, whose
  step, time, run counter and solver step count repeat those of the first: the norms of a restart
  from the first one's result; with the ladders of the first kept the loss is 11 % off),
  `mgdyn_laminate_layered_transient_delayed` (a second simulation after `RUN` of two steps whose
  first the AV solver skips: the norms of a restart from the first one's result; with the ladders
  of the first kept the loss is 48 % off), `mgdyn_laminate_layered_transient_second_av` (a
  second, air-only AV solver whose step count differs from the layered one's after only that is
  initialised again: the norms of a restart; with the other solver resetting the ladders the loss
  is 98 % low), `mgdyn_laminate_layered_transient_cubic` (an ordinary and a monotone cubic spline
  of the same
  curve in two stacks: each energy from its own table), `mgdyn_laminate_layered_transient_spd`
  (MUMPS reading one triangle, deep saturation at stacking factor 0.5: the symmetric Newton term
  gives the norms of the full one with a general solver; with the full term MUMPS's residual check
  stops the run), `mgdyn_laminate_layered_transient_sheet` (a lamination sheet: the Fatal; now the
  transient sheet test of the DEV-1606 lamination sheets, below),
  `mgdyn_transient_loss_layered_gse` / `_igse` (the fields of `mgdyn_transient_loss_start_gse` /
  `_igse` in a stack of stacking factor 0.5, the iGSE one with tilted laminations so that the flux
  has normal and in-plane parts: the GSE loss, the iGSE base and bounds of B_s on the metal
  fraction, as the law gives them to 9 digits; without layered bodies the loss is the same to the
  last bit) and `mgdyn_transient_loss_layered_mixed` (that field in a stack without a normal in a
  quarter of the cube next to plain steel: the loss of each body as the law gives it, to 9
  digits). Elmer as it is: an edge initial
  condition needs the scalar `AV {e}` entry for `SetInitialConditions` to read its components,
  the masked means of `SaveScalars` need `Parallel Reduce` on ranks without that body, and the
  `int` operator of `SaveScalars` skips elements of an elemental (`-elem`) field.
- Lamination sheets for laminations of any thickness and in the transient A-V solver (DEV-1606,
  second part), so that a core whose insulated laminations close around a conductor has one model in
  both solvers. In the TRAFOLO tube case (a copper bar through square CoNiFe tubes, each insulated,
  which carry currents circling around them of zero net value) the `layered` block alone is up to
  8.8 % off in L and 12 % in R, the block with one lamination sheet per tube within 1 % (CORE
  HOMOGENIZATION, research notes section 7g).
  - Harmonic: `LaminateSheetFactor` floors the sheet's own share of the face impedance at 1e-5 of
    the block's share, f = 1/(kd/sinh(kd) + 1e-5 (kd/2) tanh(kd/2)), with kd/sinh(kd) through
    exp(-kd) above Re(kd) = 20 and tanh = 1 above Re(kd/2) = 20 (the complex SINH and TANH of
    gfortran overflow to NaN). The impedance between the lamination faces stays within 1e-5
    |tanh(x)|^2 (at most 1.3e-5, at 2.4 skin depths) of the exact x coth(x)/(sigma d) at any
    thickness, and a lamination much thicker than the skin depth gives a near-perfect conductor
    (sigma d times 2/(1e-5 kd)) instead of the Fatal above k d = 600. Thin laminations change by
    1e-5 x tanh(x) of their face impedance (4.6e-5 of the sheet factor at 2.8 skin depths).
  - Transient: `WhitneyAVSolver` assembles lamination sheets (before: a Fatal) as resistive sheets
    of conductance sigma d, with the scalar potential of their nodes and no magnetic term as in the
    harmonic solver, integrated with the solver's own time steps like the volume conductors (BDF of
    order 1 or 2; implicit Euler, which the TRAFOLO app writes, is BDF1). With the block's share x
    tanh(x) their face impedance is 1 + x tanh(x) for the exact x coth(x) (over sigma d): 2x^2/3
    high at low and 1/x at high frequency. The sheet's own share of the exact impedance has poles of
    both signs, so no passive ladder of its own carries it; the block's ladder carries the eddy
    currents of zero net value. A steady solve leaves the sheets out (an insulated lamination
    carries no current at DC). A Fatal stops lamination sheets in an eigen analysis, with BDF orders
    above 2 or other methods, `Adaptive Timestepping`, `Convergence Control`, a predictor-corrector
    solver, `Electrodynamics Model` and a `Timestep Scale` or `Timestep Over Intervals` on the AV
    solver or `CalcFields` (the AV solver keeps a timestep of its own only during its solve).
  - Loss: `CalcFields` evaluates the transient sheet loss with E = -(a1 A^{n+1} + a2 A^n + a3
    A^{n-1})/dt - grad V, the BDF weights of the AV solver's step (`LaminateStepBDFWeights`, which
    the DEV-1606 ladder now shares): the dissipation of the discrete equations, free of the gauge of
    A. The first-order difference CalcFields takes for volume conductors is not (not changed here):
    in the gauge-free MUMPS runs of the tube case it gave the copper 52 - 20000 times its loss by
    the energy balance. A harmonic field takes the harmonic law, also in a transient simulation
    (next to a transient heat solver). As in harmonic the sheet loss is in `res: Thin sheet current
    power` and in `res: Laminate stack loss`, not in `res: Eddy current power`. Each sheet face's
    loss also goes as one density into the `Proximity Loss` of its layered parent elements, from
    which `ProcessFields` takes the heat load: in harmonic into the DG `Proximity Loss E` of
    `CalcFields` (`Calculate Homogenization Loss`), whose integral is then the laminate stack loss,
    in transient into the elemental `Proximity Loss` of the AV solver, on top of the element's
    ladder loss. The nodal `Proximity Loss` of a harmonic run leaves the sheets out.
  - Conductivity: `Laminate Stack Conductivity` and the sheets' `Thin Sheet Electric Conductivity`
    are the means of their nodal values over each element and face (`LaminateElementMean`), so that
    they may follow the temperature (`Variable Temperature`): the harmonic tensor, `CalcFields` and
    the transient ladder (built at the start of each step) take the element's, the sheets the
    face's. Constants give the same values to the last bit. `Laminate Stack Conductivity = 0` gives
    the static stack tensor (an empty ladder); the sheets keep their own `Thin Sheet Electric
    Conductivity` (0 there: no sheet current).
  - Partitions: every lamination sheet that a partition interface cuts stops with a Fatal
    (`CheckLaminateSheetPartition`), in both solvers; DEV-1604 stopped only next to layered H-B
    bodies. ElmerGrid writes such a face to both partitions, each copy with its local parent, and
    both were assembled in full: the DEV-1600 test model on 3 ranks with `-partnobcoptim` held 15 of
    the 36 faces of one sheet twice, its sheet power 2.3 % off. ElmerGrid's default partitioning
    (what TRAFOLO's `-metis N` uses) keeps the sheets inside partitions. A sheet face on the outer
    boundary whose nodes all lie on a partition interface stops too: the partitioned mesh does not
    tell it from a cut one.
  - Validated on the tube case (`hfhomog/tube_sheets.py`, `studies/run_tube_sheets.py` in CORE
    HOMOGENIZATION, results `results/trafolo_tube/tube_sheets_dev1606sheets_*.csv`). Harmonic,
    default mesh, 50 Hz - 1 GHz: the floored law within 1.7e-5 of the plain law in R up to 10 MHz,
    within 3.3e-5 of perfect-conductor sheets from 10 MHz to 1 GHz (100 MHz 1.65005 Ohm, 4.95708e-8
    H), energy balance 1 to 2e-6; against the converged meshed tubes L +0.1 .. +1.0 % and R -0.5 ..
    +0.1 % from 1 kHz to 100 MHz, R -2.2 % at 50 Hz. Transient against harmonic runs with the same
    resistive sheets (one coarse mesh per frequency, the same uniform terminal current density, BDF2
    with 200 steps per period, 3 ramp and 3 steady periods) at 50 Hz, 1, 5 and 20 kHz: R +0.02 ..
    +0.08 %, L +0.01 .. +0.16 %, stack loss +0.02 .. +0.07 %, sheet loss -0.9 .. +0.07 %, settled to
    0.005 %; the input power less the stack loss (sheets included) and the harmonic copper loss
    leaves 1.6e-4 of it at 50 Hz and at most 3.3e-5 from 1 kHz. The resistive sheet is within 0.39 %
    of the floored law in R, L and stack loss. A triangular current at 5 kHz (400 steps per period)
    against the sum of harmonic runs at its odd harmonics up to the 15th: input power +0.04 %, stack
    loss +0.04 %, sheet loss -0.27 %, residual 2.2e-5 of the input power. Implicit Euler, the
    TRAFOLO app's time stepping, at 1, 5 and 20 kHz with 400 and 800 steps per period and its
    first-order error extrapolated out: R and stack loss within 0.05 %, L within 0.15 %, the copper
    loss within 0.1 % (with implicit Euler the first-order dA/dt of CalcFields is the solver's),
    energy balance to 4e-5 of the input power. At 400 steps per period implicit Euler itself puts R
    0.3 - 0.6 % high, L 0.2 - 0.5 % low and 0.3 - 0.9 % of the input power into numerical damping.
    Limit: the transient sheet loss is 0.4 % (5 kHz) and 0.9 % (20 kHz) below the harmonic one with
    the same sheets, in both time steppings, nearly independent of the ladder order and also with a
    static block (`Laminate Stack Conductivity = 0`); both solvers' energy balances above close, so their
    discrete solutions differ there, not the loss evaluation (the sheets carry 0.3 % of the input
    power).
  - Tests: `mgdyn_laminate_layered_harmonic_thick` (100 MHz, 746 skin depths: the norm of
    perfect-conductor sheets to 7e-6; the Fatal before),
    `mgdyn_laminate_layered_harmonic_sheet_heat` (the integral of `Proximity Loss E` is the stack
    loss), `mgdyn_laminate_layered_harmonic_split_sheet` (the linear DEV-1600 model with sheets cut
    by partitions: the Fatal), `mgdyn_laminate_layered_harmonic_temperature` /
    `_transient_temperature` (conductivity following a temperature of 20 - 140 C; at 20 C everywhere
    the constant's norms to the last digit), `mgdyn_laminate_layered_harmonic_sigma0` /
    `_transient_sigma0`, `mgdyn_laminate_layered_transient_sheet` (implicit Euler as the TRAFOLO
    app: the loss, the heat field and the terminal potential; also on 2 ranks with the same norms,
    and BDF order 1 gives them too), `_transient_sheet_dt` (the variable-step BDF2 weights),
    `_transient_sheet_restart` (the BDF1 step after a restart), `_transient_sheet_split` (the Fatal)
    and `_transient_sheet_eigen` (the Fatal of an eigen analysis), all on one 616-element wedge of
    the tube case. Each fails when its feature is broken (checked by breaking the BDF weights, the
    implicit Euler weights, the step-size and start logic, the heat field and the sheet assembly in
    turn). Models without lamination sheets give the same results to the last printed digit (six
    transient, harmonic and thin-sheet tests compared).
- Random crashes and locale switches in formatted output on Windows (DEV-1612). libgfortran (MSYS2
  UCRT64, GNU Fortran 16.1, the `libgfortran-5.dll` of the bundle) wraps every formatted transfer in
  `old = setlocale(LC_NUMERIC, NULL)`, `setlocale(LC_NUMERIC, "C")` and `setlocale(LC_NUMERIC,
  old)`. The UCRT rebuilds its locale data on every set, also "C" over "C", so `old` points to freed
  memory when it is restored. Once that memory has been released the restore crashes: a SIGSEGV at
  random, the more likely the more the heap is used (the TRAFOLO tube case on 8 ranks with MUMPS
  refinement crashed in each of 5 runs, within 60 - 1957 steps, at a `WRITE` of `ComputeChange`).
  When the freed bytes start with a NUL it sets the user's Windows locale instead of "C" (in a plain
  gfortran program after 2 - 8 % of the formatted writes), and C code between two transfers (MATC,
  Lua) then parses and prints numbers with that locale's decimal separator, a comma on many European
  systems. `ElmerSolver` now points libgfortran's `setlocale` import at a wrapper that returns a
  copy of the queried name (`GfortranLocale.c`, MinGW only, needs the shared libgfortran), first
  thing in the program and in the `ElmerSolver` routine, keeps the DLL of the wrapper loaded, and
  says so in the banner (`libgfortran locale workaround installed.`; a warning, also in silent mode,
  if it cannot or finds no import of libgfortran that resolves to `setlocale`). `ViewFactors` and
  `Radiators` are not covered. With it the tube case ran through in 5 of 5 runs (2400 steps; the
  unchanged build crashed at step 1135), the plain gfortran program kept "C" in all 200 000 writes,
  and six tests give outputs identical to the last bit. Test: `robustness_gfortran_locale` (MinGW
  builds: 200 000 formatted writes inside ElmerSolver leave the C runtime's numeric locale "C", and
  the banner line; it also fails once a libgfortran no longer imports `setlocale`).

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
