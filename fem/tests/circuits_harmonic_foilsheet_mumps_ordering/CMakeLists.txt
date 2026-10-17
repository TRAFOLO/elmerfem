# METIS must be in the MUMPS build: the metis case requires INFOG(7) = 5
IF(HAVE_MUMPS AND METIS_NEEDED)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  # The mesh and the circuit of circuits_harmonic_foilsheet
  SET(FOILSHEET ${CMAKE_CURRENT_SOURCE_DIR}/../circuits_harmonic_foilsheet)
  CONFIGURE_FILE(sif/metis.sif sif/metis.sif COPYONLY)
  CONFIGURE_FILE(sif/qamd.sif sif/qamd.sif COPYONLY)
  CONFIGURE_FILE(sif/unknown.sif sif/unknown.sif COPYONLY)
  file(COPY ${FOILSHEET}/sif/6459-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY ${FOILSHEET}/6459/mesh.boundary ${FOILSHEET}/6459/mesh.header ${FOILSHEET}/6459/mesh.nodes
    ${FOILSHEET}/6459/mesh.elements ${FOILSHEET}/6459/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/6459")
  file(COPY ${FOILSHEET}/6459/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/6459/dat")

  ADD_ELMER_TEST(circuits_harmonic_foilsheet_mumps_ordering LABELS circuits harmonic 3D mgdyn whitney mumps)
ENDIF()
