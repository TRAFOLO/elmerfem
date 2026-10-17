IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)
  ADD_ELMERTEST_MODULE(circuits_transient_stranded_homogenization_energy TestEnergy TestEnergy.F90)

  CONFIGURE_FILE(sif/2241.sif sif/2241.sif COPYONLY)
  SET(base "${CMAKE_CURRENT_SOURCE_DIR}/../circuits_transient_stranded_homogenization")
  FILE(COPY ${base}/ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  FILE(COPY ${base}/sif/2241-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  FILE(COPY ${base}/2241 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(circuits_transient_stranded_homogenization_energy LABELS 3D circuits transient homogenization stranded mgdyn whitney)
ENDIF()
