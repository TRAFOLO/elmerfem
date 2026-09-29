IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/case.sif sif/case.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY mesh/mesh.boundary mesh/mesh.header mesh/mesh.nodes mesh/mesh.elements mesh/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/mesh")
  file(COPY dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/dat")

  ADD_ELMER_TEST(circuits_harmonic_foil_round LABELS circuits harmonic 3D mgdyn whitney dc)
ENDIF()
