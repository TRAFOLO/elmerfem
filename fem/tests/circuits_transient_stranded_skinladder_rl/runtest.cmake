include(test_macros)
# y0 = 0, the limit of the Sigma 33 skin ladder: a resistor K/G0 in series with
# the inductance K sigma/G0, whose heat is K I^2/G0 at any frequency. DC at two
# step sizes, a ramp with implicit Euler at two step sizes and with BDF2 at a
# fixed and at a changing step, a sine at 12 and 24 points per period. The
# sweep cases take y0 = 1e-2, 1e-4 and 1e-6 G0 at the same G0 and converge to
# the y0 = 0 sine linearly in y0/G0.
FOREACH(case dc_ie_dt2e-7 dc_ie_dt1e-6 ramp_ie_dt2e-7 ramp_ie_dt1e-7 ramp_bdf2_dt2e-7
    ramp_bdf2_vardt sine_ie_N12 sine_ie_N24 sweep_eps1e-2 sweep_eps1e-4 sweep_eps1e-6)
  FILE(READ sif/${case}.params CASE_PARAMS)
  FILE(WRITE sif/case.params "${CASE_PARAMS}")
  RUN_ELMER_TEST()
ENDFOREACH()
