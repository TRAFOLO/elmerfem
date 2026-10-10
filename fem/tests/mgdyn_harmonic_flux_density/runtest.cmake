include(test_macros)

RUN_ELMER_TEST()

# Linear case: B is final after the first solve, so the loop must stop at iteration 2.
IF(NOT TEST_STDOUT_VARIABLE MATCHES "NS \\(ITER=2\\)[^\n]*:: mgdynamics[^a-z]")
  MESSAGE(FATAL_ERROR "No second nonlinear iteration of mgdynamics")
ENDIF()
IF(TEST_STDOUT_VARIABLE MATCHES "NS \\(ITER=3\\)[^\n]*:: mgdynamics[^a-z]")
  MESSAGE(FATAL_ERROR "Flux density measure did not stop the linear case at iteration 2")
ENDIF()
