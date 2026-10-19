IF(HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(case.sif case.sif COPYONLY)
  file(COPY ELMERSOLVER_STARTINFO cube2.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(mgdyn_transient_loss_layered_mixed LABELS whitney transient quick)
ENDIF()
