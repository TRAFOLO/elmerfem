IF(HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  # The model of mgdyn_laminate_layered_transient_sheet, steady with an eigen analysis, and the
  # mesh of mgdyn_laminate_layered_harmonic_thick.
  SET(base "${CMAKE_CURRENT_SOURCE_DIR}/..")
  FILE(READ ${base}/mgdyn_laminate_layered_transient_sheet/case.sif sif)
  STRING(REPLACE "Simulation Type = Transient" "Simulation Type = Steady State" sif "${sif}")
  STRING(REPLACE "\"WhitneyAVSolver\"" "\"WhitneyAVSolver\"\n  Eigen Analysis = True\n  Eigen System Values = 2" sif "${sif}")
  FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/case.sif "${sif}")
  FILE(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  FILE(COPY ${base}/mgdyn_laminate_layered_harmonic_thick/mesh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(mgdyn_laminate_layered_transient_sheet_eigen LABELS whitney)
ENDIF()
