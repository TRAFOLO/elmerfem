IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  SET(BASE ${CMAKE_CURRENT_SOURCE_DIR}/../circuits_transient_stranded_homogenization)
  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY ${BASE}/sif/2241-circuits.definitions sif/2241.sif.in sif/N2.ladder sif/N4.ladder sif/N2bdf2.ladder
    sif/edge.ladder sif/startup.ladder sif/N2.params sif/N4.params sif/N2bdf2.params sif/edge.params
    sif/startup.params sif/mixed.sif.in sif/mixed.params sif/mixed-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY ${BASE}/2241/mesh.boundary ${BASE}/2241/mesh.header ${BASE}/2241/mesh.nodes ${BASE}/2241/mesh.elements
    ${BASE}/2241/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/2241")
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/2241/dat")

  ADD_ELMER_TEST(circuits_transient_stranded_foster_ladder LABELS 3D circuits transient homogenization stranded mgdyn whitney)
ENDIF()
