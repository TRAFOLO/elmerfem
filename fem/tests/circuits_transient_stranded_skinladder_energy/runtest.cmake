include(test_macros)
# Energy balance of the Sigma 33 skin ladder of the Hybrid Foil-Litz template
# (tau = 7.0 ns) at a switch-on and a switch-off, implicit Euler at 2 and 1 ns
# steps. Over the on-time the heat is the energy put in minus the energy left
# stored, K L chi^2/2; after the switch-off it is the stored energy released.
# The stage heat K L tau q^2 cannot be negative. Each case pins the discrete
# heat, and the continuous one within the numerical dissipation of the step,
# which halves with the step (12.5 % and 6.7 % of the stage energy).
FOREACH(case on_h2e-9 off_h2e-9 on_h1e-9 off_h1e-9)
  FILE(READ sif/${case}.params CASE_PARAMS)
  FILE(WRITE sif/case.params "${CASE_PARAMS}")
  RUN_ELMER_TEST()
ENDFOREACH()
