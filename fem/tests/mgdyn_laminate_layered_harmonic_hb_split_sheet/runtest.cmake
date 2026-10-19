include(test_macros)
# The model of mgdyn_laminate_layered_harmonic_hb on 3 ranks cut across x without the
# boundary ownership optimisation of ElmerGrid, which then splits lamination sheets between
# partitions (by default it keeps them whole). A rank holds only one parent of such a sheet,
# so the steel's reluctivity there would differ from the serial one: the run stops instead.
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 mesh -partition 3 1 1 -partnobcoptim -nooverwrite)
RUN_ELMER_EXPECT(EXIT NONZERO MATCH "lies on a partition interface" "aborts all 3 MPI ranks")
