IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/switch_diode_freewheel.sif sif/switch_diode_freewheel.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY sif/switch_diode_freewheel_circuit.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY 1381/mesh.boundary 1381/mesh.header 1381/mesh.nodes 1381/mesh.elements 1381/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/1381")
  file(COPY 1381/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/1381/dat")

  ADD_ELMER_TEST(circuits2D_transient_switch_diode_freewheel LABELS 2D circuits transient mgdyn switch diode)
ENDIF()
