IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/layerrule.sif sif/layerrule.sif COPYONLY)

  SET(GAPFRINGE ${CMAKE_CURRENT_SOURCE_DIR}/../circuits_harmonic_foilsheet_gapfringe)
  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY ${GAPFRINGE}/sif/gapfringe-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY sif/t0450.params sif/t0551.params sif/t1450.params sif/t1550.params sif/t1550_explicit3.params
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY ${GAPFRINGE}/wedge/mesh.boundary ${GAPFRINGE}/wedge/mesh.header ${GAPFRINGE}/wedge/mesh.nodes
    ${GAPFRINGE}/wedge/mesh.elements ${GAPFRINGE}/wedge/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/wedge")
  file(COPY ${GAPFRINGE}/wedge/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/wedge/dat")

  ADD_ELMER_TEST(circuits_harmonic_foilsheet_layer_rule NPROCS 1 4 LABELS circuits harmonic 3D mgdyn whitney)
ENDIF()
