include(test_macros)
# DEV-1667: the film conductivity of the coarse-air conjugate heat transfer in forced flow:
# full in still air near a wall, gone in moving air, and never below what the cells at a wall
# need to conduct like the wall formula.
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 block RESULT_VARIABLE grid_result)
if(NOT grid_result EQUAL 0)
  message(FATAL_ERROR "ElmerGrid failed: ${grid_result}")
endif()
RUN_ELMER_EXPECT(SIF film.sif EXIT ZERO MATCH "64 film calls, 0 differ from the law")
