IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/touching_height.sif sif/touching_height.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY sif/touching_height-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY rings/mesh.boundary rings/mesh.header rings/mesh.nodes rings/mesh.elements rings/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/rings")
  file(COPY rings/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/rings/dat")

  ADD_ELMER_TEST(circuits_harmonic_foilsheet_touching_height NPROCS 1 4 LABELS circuits harmonic 3D mgdyn whitney)
ENDIF()
