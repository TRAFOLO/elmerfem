IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/square.sif sif/square.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY sif/square-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY square/mesh.boundary square/mesh.header square/mesh.nodes square/mesh.elements square/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/square")
  file(COPY square/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/square/dat")

  ADD_ELMER_TEST(circuits_harmonic_foilsheet_square NPROCS 1 4 LABELS circuits harmonic 3D mgdyn whitney)
ENDIF()
