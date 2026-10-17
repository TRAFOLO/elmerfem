include(test_macros)
# A copy of the model 500 mm along x, bodies +3 and boundaries +12. The copy is
# appended after the original, so the foil winding in it is assembled after the
# stranded coil of the original.
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 6459 -clone 2 1 1 -clonesize 500 0 0 -cloneinds -out 6459x2)
RUN_ELMER_TEST()
