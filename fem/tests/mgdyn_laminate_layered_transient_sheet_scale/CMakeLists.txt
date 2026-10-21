IF(HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  # The model of mgdyn_laminate_layered_transient_sheet with "Timestep Scale = 0.5" on the AV
  # solver alone, on the mesh of mgdyn_laminate_layered_harmonic_thick. The AV solver keeps its
  # own timestep only during its solve, so CalcFields would evaluate the sheet loss with the
  # timestep of the simulation: the run stops.
  SET(base "${CMAKE_CURRENT_SOURCE_DIR}/..")
  FILE(READ ${base}/mgdyn_laminate_layered_transient_sheet/case.sif sif)
  STRING(REPLACE "\"WhitneyAVSolver\"" "\"WhitneyAVSolver\"\n  Timestep Scale = Real 0.5" sif "${sif}")
  FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/case.sif "${sif}")
  FILE(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  FILE(COPY ${base}/mgdyn_laminate_layered_harmonic_thick/mesh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(mgdyn_laminate_layered_transient_sheet_scale LABELS whitney transient)
ENDIF()
