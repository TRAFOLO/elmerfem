include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 winkel -partdual -metis ${MPIEXEC_NTASKS} 3 -nooverwrite)
# the newline keeps the complex interface (ZMumps_SolveSystem) from matching
RUN_ELMER_EXPECT(EXIT ZERO
  MATCH "\nMumps_SolveSystem: MUMPS ordering used [(]INFOG[(]7[)][)]: 0,"
        "CompareToReferenceSolution: PASSED all")
