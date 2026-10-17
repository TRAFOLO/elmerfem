IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/6459x2.sif sif/6459x2.sif COPYONLY)
  SET(base "${CMAKE_CURRENT_SOURCE_DIR}/../circuits_harmonic_foil")
  FILE(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  FILE(COPY sif/6459x2-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  FILE(COPY ${base}/6459/mesh.boundary ${base}/6459/mesh.header ${base}/6459/mesh.nodes
    ${base}/6459/mesh.elements ${base}/6459/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/6459")
  FILE(COPY ${base}/6459/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/6459/dat")

  ADD_ELMER_TEST(circuits_harmonic_foil_stranded_norotm LABELS circuits harmonic 3D mgdyn whitney)
ENDIF()
