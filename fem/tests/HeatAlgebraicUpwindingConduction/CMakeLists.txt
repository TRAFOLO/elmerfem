include(test_macros)
configure_file(case.sif case.sif COPYONLY)
file(COPY ELMERSOLVER_STARTINFO mesh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
ADD_ELMER_TEST(HeatAlgebraicUpwindingConduction LABELS cht)
