include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 wire -nooverwrite)
IF(MPIEXEC_NTASKS GREATER 1)
  execute_process(COMMAND ${ELMERGRID_BIN} 2 2 wire -partdual -metis ${MPIEXEC_NTASKS} 3 -nooverwrite)
ENDIF()
RUN_ELMER_EXPECT(EXIT ZERO
  MATCH "MUMPS drop zeros: [1-9][0-9]* of [0-9]+ entries left out"
        "CompareToReferenceSolution: PASSED all")
