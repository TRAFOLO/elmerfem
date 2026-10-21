IF(HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  # The linear model and the mesh of mgdyn_laminate_layered_harmonic.
  SET(base "${CMAKE_CURRENT_SOURCE_DIR}/../mgdyn_laminate_layered_harmonic")
  CONFIGURE_FILE(${base}/case.sif case.sif COPYONLY)
  FILE(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  FILE(COPY ${base}/mesh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(mgdyn_laminate_layered_harmonic_split_sheet NPROCS 3 LABELS whitney harmonic)
ENDIF()
