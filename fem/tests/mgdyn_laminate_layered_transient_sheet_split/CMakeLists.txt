IF(HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  # The model of mgdyn_laminate_layered_transient_sheet and the mesh of mgdyn_laminate_layered_harmonic_thick.
  SET(base "${CMAKE_CURRENT_SOURCE_DIR}/..")
  CONFIGURE_FILE(${base}/mgdyn_laminate_layered_transient_sheet/case.sif case.sif COPYONLY)
  FILE(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  FILE(COPY ${base}/mgdyn_laminate_layered_harmonic_thick/mesh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(mgdyn_laminate_layered_transient_sheet_split NPROCS 3 LABELS whitney transient)
ENDIF()
