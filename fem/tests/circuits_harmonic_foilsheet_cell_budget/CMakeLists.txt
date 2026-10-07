IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/budget.sif sif/budget.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  SET(RINGS ${CMAKE_CURRENT_SOURCE_DIR}/../circuits_harmonic_foilsheet_touching)
  file(COPY ${RINGS}/sif/touching-circuits.definitions sif/fit.params sif/overflow.params
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY ${RINGS}/rings/mesh.boundary ${RINGS}/rings/mesh.header ${RINGS}/rings/mesh.nodes
    ${RINGS}/rings/mesh.elements ${RINGS}/rings/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/rings")
  file(COPY ${RINGS}/rings/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/rings/dat")

  ADD_ELMER_TEST(circuits_harmonic_foilsheet_cell_budget NPROCS 1 4 LABELS circuits harmonic 3D mgdyn whitney)
ENDIF()
