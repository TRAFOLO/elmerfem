IF(HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE( case.sif case.sif COPYONLY)

  SET(base "${CMAKE_CURRENT_SOURCE_DIR}/../mgdyn_laminate_layered_transient")
  FILE(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  FILE(COPY ${base}/mesh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(mgdyn_laminate_layered_transient_spd LABELS whitney transient)
ENDIF()
