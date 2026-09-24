# TRAFOLO fork of ElmerFEM

This is [TRAFOLO](https://trafolo.eu)'s fork of [ElmerFEM](https://github.com/ElmerCSC/elmerfem),
maintained as the solver source for the TRAFOLO magnetics-simulation application. It is licensed
under **GPL-2.0+/LGPL-2.1**, identical to upstream Elmer (see `LICENSE.md`, `license_texts/`).

## What differs from upstream

TRAFOLO-authored additions/fixes on the `trafolo` branch (all GPL-2.0+, in
`fem/src/modules/` and `fem/src/`):

- `htc_udf` — temperature-dependent convective heat-transfer-coefficient BC.
- `ProcessFields`, `LoadFields` — transient winding-homogenization post-processing.
- `StatElecSolveVec` — thin-layer Robin coefficient fix (missing `Eps0`).
- Transient winding homogenization (CalcFields effective conductivity).
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
- Circuit preconditioner (`Linear System Preconditioning = Circuit`, new module
  `CircuitAuxPrec`): the circuit block is summed over the ranks and LU-factored densely with
  LAPACK on every rank, replacing the UMFPACK/MUMPS solve. By default the field block keeps its
  ILU and the preconditioner stays block diagonal, as before. For complex (harmonic) systems,
  `Circuit Prec Block Method = schur` uses the Schur complement (one field preconditioner
  application per circuit column coupled to the field), and for `WhitneyAVHarmonicSolver`
  `Circuit Prec Field Method = "auxiliary space"` wraps ILU smoothing around a PRESB correction
  with the real edge operator `K + omega*sigma*M`, solved by Hypre AMS in the helper solver
  `HarmonicEdgePrecSolver` (needs `WITH_Hypre`; solver sections as in the test sif).
  `Linear System Scaling Row Norm Fallback = True` scales the zero-diagonal circuit rows of
  complex systems by their row norm, so they do not dominate the stopping test. A Hypre solve
  nested in a parallel Krylov solve now restores the caller's `GlobalData` (`SParIterSolver`).
  Tests: `circuits_harmonic_massive_schur`, `circuits_harmonic_massive_ams` (Hypre only).

## Building

The TRAFOLO product is a **Windows-only** solver bundle built with the standard Elmer CMake
system (MSYS2 UCRT64 toolchain). The direct solver is Elmer's vendored UMFPACK 4.4. Hypre,
built from source without SuperLU_dist and ParMETIS, is needed for the harmonic
auxiliary-space circuit preconditioner (`WITH_Hypre=ON`). Build and
packaging tooling is maintained internally by TRAFOLO and is **not** part of this repository;
this repo carries source only. There is no CI here — validation is run locally.

Corresponding source for any distributed TRAFOLO Elmer binary is this repository at the commit
recorded in the bundle's `SOURCE.txt` (GPL compliance).
