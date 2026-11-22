include(test_macros)
# DEV-1667: the compiled wall formulas of the coarse-air conjugate heat transfer give the
# values of the formulas they replace, for every formula, temperature difference and height.
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 block)
RUN_ELMER_EXPECT(SIF wall_law.sif EXIT ZERO MATCH "542 wall law calls, 0 differ from the formulas")
