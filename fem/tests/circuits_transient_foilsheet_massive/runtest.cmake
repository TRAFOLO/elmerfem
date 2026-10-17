include(test_macros)
# The massive coil sits in a copy of the model 500 mm along x, bodies +3 and
# boundaries +12, so both coils share one mesh without touching.
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 6479 -clone 2 1 1 -clonesize 500 0 0 -cloneinds -out 6479x2)
RUN_ELMER_TEST()
