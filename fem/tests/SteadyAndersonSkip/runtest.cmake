include(test_macros)
file(READ case.sif case_text)
if(MPIEXEC_NTASKS GREATER 1)
  # Two ranks. Every rank has to keep the same solvers out of the acceleration: the scalar
  # products of the others are taken over all ranks.
  execute_process(COMMAND ${ELMERGRID_BIN} 1 2 channel.grd -partition ${MPIEXEC_NTASKS} 1 1
                  RESULT_VARIABLE grid_result)
  if(NOT grid_result EQUAL 0)
    message(FATAL_ERROR "ElmerGrid failed: ${grid_result}")
  endif()
  string(REPLACE "Linear System Direct Method = UMFPACK" "Linear System Direct Method = MUMPS"
         agree_text "${case_text}")
  file(WRITE agree.sif "${agree_text}")
  RUN_ELMER_EXPECT(SIF agree.sif EXIT ZERO MATCH "accelerated=T")
  # A flag that differs from rank to rank stops the run before a rank takes a path of its own.
  get_filename_component(here "." ABSOLUTE)
  if(EXISTS "${here}/rank_dependent_input")
    string(REPLACE "Steady State Anderson Skip = Logical True"
           "Steady State Anderson Skip = Logical # ELMER_PARALLEL.pe == 0" differ_text "${agree_text}")
    file(WRITE differ.sif "${differ_text}")
    RUN_ELMER_EXPECT(SIF differ.sif EXIT NONZERO
                     MATCH "Steady State Anderson Skip of solver 2 differs between the ranks")
  endif()
else()
  execute_process(COMMAND ${ELMERGRID_BIN} 1 2 channel.grd RESULT_VARIABLE grid_result)
  if(NOT grid_result EQUAL 0)
    message(FATAL_ERROR "ElmerGrid failed: ${grid_result}")
  endif()
  # A solver with a field of its own stays out of the acceleration, which goes on for the flow:
  # the tracer after the last coupled iteration is the one its own solve left there.
  macro(LAST_TWO_RECORDS sweep final)
    file(STRINGS sweeps.dat _sweeps)
    list(GET _sweeps -1 ${sweep})
    file(STRINGS final.dat _final)
    list(GET _final -1 ${final})
    string(STRIP "${${sweep}}" ${sweep})
    string(STRIP "${${final}}" ${final})
  endmacro()
  file(REMOVE sweeps.dat final.dat)
  # the last coupled iteration (the fifth) is one the acceleration accepted
  RUN_ELMER_EXPECT(SIF case.sif EXIT ZERO MATCH "iteration=5 scaled map residual= *[-+.0-9E]+ accelerated=T")
  LAST_TWO_RECORDS(kept_sweep kept_final)
  if(NOT kept_sweep STREQUAL kept_final)
    message(FATAL_ERROR "The tracer kept out was changed after its last solve: '${kept_sweep}' -> '${kept_final}'")
  endif()
  # Control: without the flag the acceleration does change the tracer, so the check above can fail.
  string(REPLACE "  Steady State Anderson Skip = Logical True\n" "" taken_text "${case_text}")
  file(WRITE taken.sif "${taken_text}")
  file(REMOVE sweeps.dat final.dat)
  RUN_ELMER_EXPECT(SIF taken.sif EXIT ZERO MATCH "iteration=5 scaled map residual= *[-+.0-9E]+ accelerated=T")
  LAST_TWO_RECORDS(taken_sweep taken_final)
  if(taken_sweep STREQUAL taken_final)
    message(FATAL_ERROR "The control run shows no change of the tracer by the acceleration: '${taken_final}'")
  endif()
  # The same solver on the first velocity component: that field is a part of the accelerated
  # flow field, which would be extrapolated all the same. The pair is refused.
  string(REPLACE "Variable = Tracer" "Variable = Velocity 1" refused_text "${case_text}")
  file(WRITE refused.sif "${refused_text}")
  RUN_ELMER_EXPECT(SIF refused.sif EXIT NONZERO
                   MATCH "A solver kept out of the acceleration shares its variable with an accelerated one")
endif()
