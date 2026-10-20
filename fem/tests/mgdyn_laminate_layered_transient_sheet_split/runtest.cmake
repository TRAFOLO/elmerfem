include(test_macros)
# The model of mgdyn_laminate_layered_transient_sheet on 3 ranks cut across x without the
# boundary ownership optimisation of ElmerGrid, which then writes sheet faces whose parents it
# puts in two partitions to both of them: the sheet would count twice there, so the run stops.
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 mesh -partition 3 1 1 -partnobcoptim -nooverwrite)
RUN_ELMER_EXPECT(EXIT NONZERO MATCH "lies on a partition interface" "aborts all 3 MPI ranks")
