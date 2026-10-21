include(test_macros)
RUN_ELMER_EXPECT(EXIT NONZERO MATCH "support neither .Timestep Scale. nor .Timestep Over Intervals. on the AV solver")
