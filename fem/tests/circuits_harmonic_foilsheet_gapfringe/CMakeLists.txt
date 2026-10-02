IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/gapfringe.sif sif/gapfringe.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY sif/gapfringe-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY wedge/mesh.boundary wedge/mesh.header wedge/mesh.nodes wedge/mesh.elements wedge/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/wedge")
  file(COPY wedge/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/wedge/dat")

  ADD_ELMER_TEST(circuits_harmonic_foilsheet_gapfringe NPROCS 1 4 LABELS circuits harmonic 3D mgdyn whitney)
ENDIF()
