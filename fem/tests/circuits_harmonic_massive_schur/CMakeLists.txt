INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(sif/1962.sif sif/1962.sif COPYONLY)

SET(MASSIVE ${CMAKE_CURRENT_SOURCE_DIR}/../circuits_harmonic_massive)
file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
file(COPY ${MASSIVE}/sif/1962-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
file(COPY ${MASSIVE}/1962/mesh.boundary ${MASSIVE}/1962/mesh.header ${MASSIVE}/1962/mesh.nodes
  ${MASSIVE}/1962/mesh.elements ${MASSIVE}/1962/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/1962")

ADD_ELMER_TEST(circuits_harmonic_massive_schur LABELS circuits harmonic 3D mgdyn whitney)
