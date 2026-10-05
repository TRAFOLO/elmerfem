IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  SET(DC ${CMAKE_CURRENT_SOURCE_DIR}/../circuits_transient_stranded_skinladder_dc)
  file(COPY ${DC}/sif/2241-circuits.definitions sif/model.inc sif/y0_negative.sif
    sif/alpha_negative.sif sif/resistance.sif sif/resistance_alpha0.sif sif/restart_write.sif
    sif/restart_no_time.sif sif/restart_time0.sif DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY ${DC}/2241/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/2241/dat")
  file(COPY ${DC}/2241/mesh.boundary ${DC}/2241/mesh.header ${DC}/2241/mesh.nodes
    ${DC}/2241/mesh.elements ${DC}/2241/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/2241")

  ADD_ELMER_TEST(circuits_transient_stranded_skinladder_fatal LABELS 3D circuits transient homogenization stranded mgdyn whitney)
ENDIF()
