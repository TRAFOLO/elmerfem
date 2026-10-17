include(test_macros)
RUN_ELMER_EXPECT(EXIT ZERO NOMATCH "NOT CONVERGED: coupled")
