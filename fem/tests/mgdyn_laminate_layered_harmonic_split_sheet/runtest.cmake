include(test_macros)
# The linear model of mgdyn_laminate_layered_harmonic on 3 ranks cut across x without the
# boundary ownership optimisation of ElmerGrid, which then writes 15 of the 36 faces of one sheet
# to two partitions: assembled in both, their conductance doubled (sheet power -2.3 % before this
# stop), so the run stops.
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 mesh -partition 3 1 1 -partnobcoptim -nooverwrite)
RUN_ELMER_EXPECT(EXIT NONZERO MATCH "lies on a partition interface" "aborts all 3 MPI ranks")
