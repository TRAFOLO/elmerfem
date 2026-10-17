include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 rings_touch.grd)
RUN_ELMER_TEST()
