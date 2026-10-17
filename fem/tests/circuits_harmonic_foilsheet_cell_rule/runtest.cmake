include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 wedge -metiskway ${MPIEXEC_NTASKS} -nooverwrite)
# One run per case. RUN_ELMER_TEST checks the pinned norm; each "! expect:" line
# of the case must also appear in the solver output: the final layout and the
# block element count, which must not change with the partitioning.
FOREACH(case f50_s20 f50_cap f50_cells15 dc scan_dc_50 f50_n1200 scan_n1200)
  FILE(READ sif/${case}.params CASE_PARAMS)
  FILE(WRITE sif/case.params "${CASE_PARAMS}")
  RUN_ELMER_TEST()
  STRING(REGEX MATCHALL "! expect: [^\r\n]*" EXPECTS "${CASE_PARAMS}")
  FOREACH(expect IN LISTS EXPECTS)
    STRING(REPLACE "! expect: " "" expect "${expect}")
    STRING(FIND "${TEST_STDOUT_VARIABLE}" "${expect}" POS)
    IF(POS EQUAL -1)
      MESSAGE(FATAL_ERROR "Case ${case}: \"${expect}\" not in the solver output")
    ENDIF()
  ENDFOREACH()
ENDFOREACH()
