IF(HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

  # The model of mgdyn_laminate_layered_harmonic_sheet_heat as one step of a transient
  # simulation (as next to a transient heat solver), on the mesh of
  # mgdyn_laminate_layered_harmonic_thick. CalcFields is called with Transient = True for a
  # harmonic field: the sheet loss, its heat and so the norms stay those of the steady run.
  SET(base "${CMAKE_CURRENT_SOURCE_DIR}/..")
  FILE(READ ${base}/mgdyn_laminate_layered_harmonic_sheet_heat/case.sif sif)
  STRING(REPLACE "Simulation Type = Steady" "Simulation Type = Transient\n  Timestepping Method = BDF\n  BDF Order = 1\n  Timestep Sizes = Real 1.0\n  Timestep Intervals = Integer 1\n  Output Intervals = 0" sif "${sif}")
  FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/case.sif "${sif}")
  FILE(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  FILE(COPY ${base}/mgdyn_laminate_layered_harmonic_thick/mesh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  ADD_ELMER_TEST(mgdyn_laminate_layered_harmonic_sheet_heat_transient LABELS whitney harmonic)
ENDIF()
