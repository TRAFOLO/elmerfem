include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 mesh -partition 2 1 1)
RUN_ELMER_TEST()
