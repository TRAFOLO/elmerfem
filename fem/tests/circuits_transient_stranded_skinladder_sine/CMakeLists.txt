IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/2241.sif sif/2241.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY sif/2241-circuits.definitions sif/N12.params sif/N24.params sif/N48.params sif/N12_iter2.params DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY 2241/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/2241/dat")
  file(COPY 2241/mesh.boundary 2241/mesh.header 2241/mesh.nodes 2241/mesh.elements 2241/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/2241")

  ADD_ELMER_TEST(circuits_transient_stranded_skinladder_sine LABELS 3D circuits transient homogenization stranded mgdyn whitney)
ENDIF()
