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

## Building

The TRAFOLO product is a **Windows-only** solver bundle built with the standard Elmer CMake
system (MSYS2 UCRT64 toolchain). The direct solver is Elmer's vendored UMFPACK 4.4. Build and
packaging tooling is maintained internally by TRAFOLO and is **not** part of this repository;
this repo carries source only. There is no CI here — validation is run locally.

Corresponding source for any distributed TRAFOLO Elmer binary is this repository at the commit
recorded in the bundle's `SOURCE.txt` (GPL compliance).
