include(test_macros)
# The two-coil mesh: the mesh of circuits_transient_foilsheet and a copy of it
# 500 mm along x, bodies +3 and boundaries +12, so the coils do not touch.
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 6479 -clone 2 1 1 -clonesize 500 0 0 -cloneinds -out 6479x2)
# Cut along x: the two coils lie side by side along x, so every partition
# holds the elements of one coil only.
IF(${MPIEXEC_NTASKS} GREATER 1)
  execute_process(COMMAND ${ELMERGRID_BIN} 2 2 6479x2 -partition ${MPIEXEC_NTASKS} 1 1 -nooverwrite)
ENDIF()
RUN_ELMER_TEST()
