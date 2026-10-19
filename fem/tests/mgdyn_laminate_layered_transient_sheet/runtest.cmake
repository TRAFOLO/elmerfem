include(test_macros)
RUN_ELMER_EXPECT(EXIT NONZERO MATCH "available in the harmonic solver only")
