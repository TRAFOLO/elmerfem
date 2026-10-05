IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  SET(DC ${CMAKE_CURRENT_SOURCE_DIR}/../circuits_transient_stranded_skinladder_dc)
  SET(RL ${CMAKE_CURRENT_SOURCE_DIR}/../circuits_transient_stranded_skinladder_rl)
  CONFIGURE_FILE(${RL}/sif/2241.sif sif/2241.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY ${DC}/sif/2241-circuits.definitions sif/vcomp.check sif/on_h2e-9.params sif/off_h2e-9.params
    sif/on_h1e-9.params sif/off_h1e-9.params DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY ${DC}/2241/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/2241/dat")
  file(COPY ${DC}/2241/mesh.boundary ${DC}/2241/mesh.header ${DC}/2241/mesh.nodes
    ${DC}/2241/mesh.elements ${DC}/2241/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/2241")

  ADD_ELMER_TEST(circuits_transient_stranded_skinladder_energy LABELS 3D circuits transient homogenization stranded mgdyn whitney)
ENDIF()
