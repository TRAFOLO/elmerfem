include(test_macros)
# on 2 ranks: Metis partitions of the default ElmerGrid run, which keeps the lamination sheets
# inside partitions
IF(MPIEXEC_NTASKS GREATER 1)
  execute_process(COMMAND ${ELMERGRID_BIN} 2 2 mesh -metiskway ${MPIEXEC_NTASKS} -nooverwrite)
ENDIF()
RUN_ELMER_TEST()
