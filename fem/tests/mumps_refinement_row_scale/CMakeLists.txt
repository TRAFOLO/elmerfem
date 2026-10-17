IF(HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(case.sif case.sif COPYONLY)
  ADD_ELMERTEST_MODULE(mumps_refinement_row_scale MumpsRowScaleTest MumpsRowScaleTest.F90)
  file(COPY ELMERSOLVER_STARTINFO square.grd MumpsRowScaleTest.F90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(mumps_refinement_row_scale LABELS mumps)
ENDIF()
