include(test_macros)
# DC current through a stranded coil whose Sigma 33 skin ladder has a memory,
# at two step sizes with BDF2 and once with implicit Euler. Every run must give
# the DC resistance K/G0, whatever the step.
FOREACH(case bdf2_dt2e-7 bdf2_dt1e-6 ie_dt2e-7)
  FILE(READ sif/${case}.params CASE_PARAMS)
  FILE(WRITE sif/case.params "${CASE_PARAMS}")
  RUN_ELMER_TEST()
ENDFOREACH()
