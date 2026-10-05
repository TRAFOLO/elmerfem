include(test_macros)
# Input safeguards of the Sigma 33 skin ladder (DEV-1549): each case must stop
# with its Fatal, or run to the end. The restart cases read the result file
# that restart_write leaves in the mesh directory.
SET(FAILED_CASES "")

MACRO(EXPECT_FATAL case text)
  EXECUTE_ELMER_SOLVER(sif/${case}.sif)
  FILE(READ "sif/${case}.sif-stdout.log" CASE_OUT)
  STRING(FIND "${CASE_OUT}" "ERROR:: InitStrandedSkinLadder: " POS_ERR)
  STRING(FIND "${CASE_OUT}" "${text}" POS_TEXT)
  IF(POS_ERR EQUAL -1 OR POS_TEXT EQUAL -1)
    LIST(APPEND FAILED_CASES ${case})
  ENDIF()
ENDMACRO()

MACRO(EXPECT_RUN case)
  EXECUTE_ELMER_SOLVER(sif/${case}.sif)
  FILE(READ "sif/${case}.sif-stdout.log" CASE_OUT)
  STRING(FIND "${CASE_OUT}" "ERROR::" POS_ERR)
  STRING(FIND "${CASE_OUT}" "ALL DONE" POS_DONE)
  IF(NOT POS_ERR EQUAL -1 OR POS_DONE EQUAL -1)
    LIST(APPEND FAILED_CASES ${case})
  ENDIF()
ENDMACRO()

EXPECT_FATAL(y0_negative "is negative. It is the conductivity of the ladder at high frequency")
EXPECT_FATAL(alpha_negative "is negative: the conductivity would rise with frequency")
EXPECT_FATAL(resistance "an explicit \"Resistance\" replaces the skin ladder")
EXPECT_RUN(resistance_alpha0)
EXPECT_RUN(restart_write)
EXPECT_FATAL(restart_no_time "restarted without \"Restart Time\"")
EXPECT_RUN(restart_time0)

IF(FAILED_CASES)
  MESSAGE(FATAL_ERROR "Test failed: ${FAILED_CASES}")
ELSE()
  FILE(WRITE "TEST.PASSED" "1")
ENDIF()
