include(test_macros)
# on several ranks: Metis partitions, each with stack and air elements
IF(MPIEXEC_NTASKS GREATER 1)
  execute_process(COMMAND ${ELMERGRID_BIN} 2 2 mesh -metiskway ${MPIEXEC_NTASKS} -nooverwrite)
ENDIF()
RUN_ELMER_TEST()
