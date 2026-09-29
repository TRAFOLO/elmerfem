include(test_macros)
# The same model with one and with two coupled iterations per timestep,
# both pinned to the same value: a second pass through the solvers in a
# timestep must redo the step, not take another one.
FOREACH(case iter1 iter2)
  FILE(READ sif/${case}.params CASE_PARAMS)
  FILE(WRITE sif/case.params "${CASE_PARAMS}")
  RUN_ELMER_TEST()
ENDFOREACH()
