IF(HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(converged.sif converged.sif COPYONLY)
  CONFIGURE_FILE(limit.sif limit.sif COPYONLY)
  CONFIGURE_FILE(stagnated.sif stagnated.sif COPYONLY)
  file(COPY wire.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(robustness_mumps_refinement LABELS robustness mumps whitney)
ENDIF()
