INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)
ADD_ELMERTEST_MODULE(circuits_transient_foilsheet_energy TestEnergy
  ../circuits_transient_stranded_homogenization_energy/TestEnergy.F90)

CONFIGURE_FILE(sif/energy.sif sif/energy.sif COPYONLY)
SET(base "${CMAKE_CURRENT_SOURCE_DIR}/../circuits_transient_foilsheet_skinloss")
FILE(COPY ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")
FILE(COPY sif/ladders.params sif/skin.params sif/beta.params sif/ladders.component
  sif/skin.component sif/beta.component DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
FILE(COPY ${base}/sif/6479-circuits.definitions ${base}/sif/S1Harmonic.dat
  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sif")
FILE(COPY ${base}/6479/mesh.boundary ${base}/6479/mesh.header ${base}/6479/mesh.nodes
  ${base}/6479/mesh.elements ${base}/6479/mesh.names DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/6479")

ADD_ELMER_TEST(circuits_transient_foilsheet_energy NPROCS 1 4 LABELS circuits transient 3D mgdyn whitney)
