IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/6459r.sif sif/6459r.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY sif/6459r-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY 6459r/mesh.boundary 6459r/mesh.header 6459r/mesh.nodes 6459r/mesh.elements 6459r/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/6459r")
  file(COPY 6459r/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/6459r/dat")

  ADD_ELMER_TEST(circuits_harmonic_massive_closed NPROCS 1 4 LABELS circuits harmonic 3D mgdyn whitney)
ENDIF()
