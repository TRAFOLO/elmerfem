include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 1962 -partition 1 4 1 -nooverwrite)
RUN_ELMER_EXPECT(EXIT NONZERO MATCH "Non existent Coil Type")
