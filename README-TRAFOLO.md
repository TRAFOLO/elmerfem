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

## Building

The TRAFOLO product is a **Windows-only** solver bundle built with the standard Elmer CMake
system (MSYS2 UCRT64 toolchain). The direct solver is Elmer's vendored UMFPACK 4.4. Build and
packaging tooling is maintained internally by TRAFOLO and is **not** part of this repository;
this repo carries source only. There is no CI here — validation is run locally.

Corresponding source for any distributed TRAFOLO Elmer binary is this repository at the commit
recorded in the bundle's `SOURCE.txt` (GPL compliance).
