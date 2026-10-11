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
| `build` | `build_msys2.sh` + `ninja` | Configures with the distributed flag set (UMFPACK on, MUMPS 4.10.0 linked statically, MS-MPI, OpenBLAS, Lua), builds |
| `deploy` | `deploy_msys2.sh` | `ninja install` + copies the UCRT64 runtime DLL closure into `bin\` |
| `prune` | `prune_install.sh` | Keeps the 18 solver modules the TRAFOLO app uses, drops dev tools, headers and import libs |
| `audit` | `license_audit.sh` | Fails on any GPL-incompatible or unlisted DLL in the import closure; writes `licenses\` and `SOURCE.txt` |
| `gates` | `run_gates.sh` | Runs the validation cases against the pruned install (inputs in `trafolo_bundle/validation/`); all must pass |
| `package` | `package_bundle.sh` | Refuses an unpruned or unaudited install, zips it |

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
