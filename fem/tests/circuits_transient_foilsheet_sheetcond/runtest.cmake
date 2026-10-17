include(test_macros)
# One sheet at two winding temperatures and two loss factors. Each run pins
# its loss at the ratio the sheet conductivity law predicts, so the test
# fails with the first run that is off.
FOREACH(case t20_lf2 t100_lf2 t20_lf25 t100_lf25)
  FILE(READ sif/${case}.params CASE_PARAMS)
  FILE(WRITE sif/case.params "${CASE_PARAMS}")
  RUN_ELMER_TEST()
ENDFOREACH()
