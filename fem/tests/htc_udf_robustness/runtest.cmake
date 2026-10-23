include(test_macros)
# DEV-1622: the convection UDF (htc_udf) never stops a run, whatever its input.
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 block)
# Losses that heat the surface far above the 1500 K range of the model: one warning, the run finishes.
RUN_ELMER_EXPECT(SIF hot.sif EXIT ZERO MATCH_ONCE "above the range of the convection model"
  MATCH "PASSED all 1 tests" NOMATCH "out of plausible range")
# Unphysical or missing parameters on every face and a heat sink far below the range: the run finishes.
RUN_ELMER_EXPECT(SIF hostile.sif EXIT ZERO MATCH "below the range of the convection model"
  "Ambient temperature .* is outside the range" "characteristic length" "Air velocity"
  "convection parameter .* is missing" "PASSED all 1 tests")
# The function itself: on the same faces at temperatures from -HUGE to +HUGE C, NaN and +-Infinity,
# and with each parameter set to NaN, +-Infinity, +-HUGE, -1, 0, 1e-300 or 1e300 for every orientation.
# The coefficient must come from the corrected inputs, never from the last-resort fallback.
RUN_ELMER_EXPECT(SIF sweep.sif EXIT ZERO MATCH "102 face calls, 0 not finite or negative"
  "74970 parameter calls, 0 not finite or negative" "Reference h matches"
  "not a number, h taken at the ambient temperature" "Inclination .* is not a number"
  NOMATCH "gave no finite h")
