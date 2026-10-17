include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 wire)
RUN_ELMER_EXPECT(SIF converged.sif EXIT ZERO
  MATCH_ONCE "MUMPS iterative refinement: [0-9]+ corrections, componentwise backward error [^(]*[(]converged[)]"
  NOMATCH "WARNING:: [A-Za-z]*MumpsRefineSolution")
RUN_ELMER_EXPECT(SIF limit.sif EXIT ZERO
  MATCH_ONCE "MUMPS iterative refinement: 1 corrections, componentwise backward error [^(]*[(]iteration limit[)]"
  MATCH_ONCE "WARNING:: ZMumpsRefineSolution: MUMPS iterative refinement iteration limit at a componentwise backward error")
RUN_ELMER_EXPECT(SIF stagnated.sif EXIT ZERO
  MATCH_ONCE "WARNING:: ZMumpsRefineSolution: MUMPS iterative refinement (stagnated|diverging) at a componentwise backward error")
