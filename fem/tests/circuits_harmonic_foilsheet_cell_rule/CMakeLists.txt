IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/cellrule.sif sif/cellrule.sif COPYONLY)

  SET(GAPFRINGE ${CMAKE_CURRENT_SOURCE_DIR}/../circuits_harmonic_foilsheet_gapfringe)
  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY ${GAPFRINGE}/sif/gapfringe-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY sif/f50_s20.params sif/f50_cap.params sif/f50_cells15.params sif/dc.params sif/scan_dc_50.params sif/f50_n1200.params sif/scan_n1200.params
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY ${GAPFRINGE}/wedge/mesh.boundary ${GAPFRINGE}/wedge/mesh.header ${GAPFRINGE}/wedge/mesh.nodes
    ${GAPFRINGE}/wedge/mesh.elements ${GAPFRINGE}/wedge/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/wedge")
  file(COPY ${GAPFRINGE}/wedge/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/wedge/dat")

  ADD_ELMER_TEST(circuits_harmonic_foilsheet_cell_rule NPROCS 1 4 LABELS circuits harmonic 3D mgdyn whitney)
ENDIF()
