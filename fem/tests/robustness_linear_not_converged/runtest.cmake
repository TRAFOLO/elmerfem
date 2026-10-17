include(test_macros)
RUN_ELMER_EXPECT(EXIT ZERO MATCH "NOT CONVERGED: linear solver=\"mgdynamics\" iterations=3 tolerance=")
