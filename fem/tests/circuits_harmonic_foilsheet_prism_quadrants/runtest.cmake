include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 potcore -partition 2 2 1 -nooverwrite)
RUN_ELMER_TEST()
