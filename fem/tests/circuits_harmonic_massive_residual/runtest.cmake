include(test_macros)

RUN_ELMER_TEST()

# The first iterate has a residual, the second (linear case) has none.
IF(NOT TEST_STDOUT_VARIABLE MATCHES "NS \\(ITER=2\\)[^\n]*:: mgdynamics[^a-z]")
  MESSAGE(FATAL_ERROR "The residual measure stopped the loop at the first iteration")
ENDIF()
IF(TEST_STDOUT_VARIABLE MATCHES "NS \\(ITER=3\\)[^\n]*:: mgdynamics[^a-z]")
  MESSAGE(FATAL_ERROR "The residual measure did not converge at the second iteration")
ENDIF()
