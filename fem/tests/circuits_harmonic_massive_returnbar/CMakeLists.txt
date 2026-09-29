IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/returnbar.sif sif/returnbar.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY sif/returnbar-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY returnbar/mesh.boundary returnbar/mesh.header returnbar/mesh.nodes returnbar/mesh.elements returnbar/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/returnbar")
  file(COPY returnbar/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/returnbar/dat")

  ADD_ELMER_TEST(circuits_harmonic_massive_returnbar NPROCS 1 4 LABELS circuits harmonic 3D mgdyn whitney)
ENDIF()
