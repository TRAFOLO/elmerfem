include(test_macros)
# Sine current through a stranded coil whose Sigma 33 skin ladder has a memory,
# at 12, 24 and 48 points per period with implicit Euler. The averaged winding
# heat must converge to the harmonic value to first order. N12_iter2 is N12 with
# two coupled iterations per step, pinned to the same values.
FOREACH(case N12 N24 N48 N12_iter2)
  FILE(READ sif/${case}.params CASE_PARAMS)
  FILE(WRITE sif/case.params "${CASE_PARAMS}")
  RUN_ELMER_TEST()
ENDFOREACH()
