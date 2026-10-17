IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  CONFIGURE_FILE(sif/2241.sif sif/2241.sif COPYONLY)

  SET(DC ${CMAKE_CURRENT_SOURCE_DIR}/../circuits_transient_stranded_skinladder_dc)
  file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY ${DC}/sif/2241-circuits.definitions sif/vcomp.check sif/dc_ie_dt2e-7.params sif/dc_ie_dt1e-6.params
    sif/ramp_ie_dt2e-7.params sif/ramp_ie_dt1e-7.params sif/ramp_bdf2_dt2e-7.params
    sif/ramp_bdf2_vardt.params sif/sine_ie_N12.params sif/sine_ie_N24.params
    sif/sweep_eps1e-2.params sif/sweep_eps1e-4.params sif/sweep_eps1e-6.params
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY ${DC}/2241/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/2241/dat")
  file(COPY ${DC}/2241/mesh.boundary ${DC}/2241/mesh.header ${DC}/2241/mesh.nodes
    ${DC}/2241/mesh.elements ${DC}/2241/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/2241")

  ADD_ELMER_TEST(circuits_transient_stranded_skinladder_rl LABELS 3D circuits transient homogenization stranded mgdyn whitney)
ENDIF()
