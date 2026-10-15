IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)
  ADD_ELMERTEST_MODULE(circuits_harmonic_foilsheet_energy TestEnergy
    ../circuits_harmonic_stranded_homogenization_energy/TestEnergy.F90)

  CONFIGURE_FILE(sif/square.sif sif/square.sif COPYONLY)
  SET(base "${CMAKE_CURRENT_SOURCE_DIR}/../circuits_harmonic_foilsheet_square")
  FILE(COPY ${base}/ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  FILE(COPY ${base}/sif/square-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  FILE(COPY ${base}/square DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(circuits_harmonic_foilsheet_energy LABELS circuits harmonic 3D mgdyn whitney)
ENDIF()
