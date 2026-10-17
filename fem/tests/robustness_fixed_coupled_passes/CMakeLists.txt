IF(HAVE_UMFPACK OR HAVE_MUMPS)
  INCLUDE(test_macros)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)
  # The inputs of mgdyn2D_scan_homogenization_elementary_solutions: per scan step three
  # coupled passes (min = max = 3) that solve three different excitations.
  SET(_base ${CMAKE_CURRENT_SOURCE_DIR}/../mgdyn2D_scan_homogenization_elementary_solutions)
  file(COPY ${_base}/ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
  file(COPY ${_base}/sif/6801.sif ${_base}/sif/6801-circuits.definitions DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
  file(COPY ${_base}/6801/mesh.boundary ${_base}/6801/mesh.header ${_base}/6801/mesh.nodes ${_base}/6801/mesh.elements ${_base}/6801/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/6801")
  file(COPY ${_base}/6801/dat/gittrack DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/6801/dat")
  ADD_ELMER_TEST(robustness_fixed_coupled_passes LABELS robustness)
ENDIF()
