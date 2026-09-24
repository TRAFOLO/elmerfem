IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/helix2.sif sif/helix2.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY sif/helix2-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY helix2/mesh.boundary helix2/mesh.header helix2/mesh.nodes helix2/mesh.elements helix2/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/helix2")
  file(COPY helix2/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/helix2/dat")

  ADD_ELMER_TEST(circuits_harmonic_massive_multiturn NPROCS 1 4 LABELS circuits harmonic 3D mgdyn whitney)
ENDIF()
