include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 rings -metiskway ${MPIEXEC_NTASKS} -nooverwrite)
RUN_ELMER_TEST()
