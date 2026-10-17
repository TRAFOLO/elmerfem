INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(case--hyphen.sif case--hyphen.sif COPYONLY)
file(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
# The unit cube of mgdyn_harmonic_loss
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/../mgdyn_harmonic_loss/mesh3d" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(mgdyn_harmonic_loss_totals LABELS whitney vtu)
