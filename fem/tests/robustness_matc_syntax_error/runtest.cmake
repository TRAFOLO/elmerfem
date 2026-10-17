include(test_macros)
RUN_ELMER_EXPECT(EXIT NONZERO MATCH "Solver input file error: MATC ERROR: [^ ]")
