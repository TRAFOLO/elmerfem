include(test_macros)
# on 3 ranks the mesh is cut along z: a rank of stack elements, one of both and one of air
# alone, which holds no layered element
IF(MPIEXEC_NTASKS GREATER 1)
  execute_process(COMMAND ${ELMERGRID_BIN} 2 2 mesh -partition 1 1 ${MPIEXEC_NTASKS} -nooverwrite)
ENDIF()
RUN_ELMER_TEST()
