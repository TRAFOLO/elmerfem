IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/touching_partition.sif sif/touching_partition.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  SET(RINGS ${CMAKE_CURRENT_SOURCE_DIR}/../circuits_harmonic_foilsheet_touching)
  file(COPY ${RINGS}/sif/touching-circuits.definitions sif/order_ab.params sif/order_ba.params
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY ${RINGS}/rings/mesh.boundary ${RINGS}/rings/mesh.header ${RINGS}/rings/mesh.nodes
    ${RINGS}/rings/mesh.elements ${RINGS}/rings/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/rings")
  file(COPY ${RINGS}/rings/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/rings/dat")

  ADD_ELMER_TEST(circuits_harmonic_foilsheet_touching_partition NPROCS 1 4 LABELS circuits harmonic 3D mgdyn whitney)
ENDIF()
