IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/potcore.sif sif/potcore.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY sif/potcore-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY potcore/mesh.boundary potcore/mesh.header potcore/mesh.nodes potcore/mesh.elements potcore/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/potcore")
  file(COPY potcore/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/potcore/dat")

  ADD_ELMER_TEST(circuits_harmonic_foilsheet_prism NPROCS 1 4 LABELS circuits harmonic 3D mgdyn whitney)
ENDIF()
