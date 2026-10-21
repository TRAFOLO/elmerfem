include(test_macros)
# The model of mgdyn_laminate_layered_transient_sheet as an eigen analysis of WhitneyAVSolver,
# which would leave the lamination sheets out of the equations: the run stops.
RUN_ELMER_EXPECT(EXIT NONZERO MATCH "not supported in an eigen analysis")
