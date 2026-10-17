include(test_macros)
# The magnetic energy of a transient foil sheet against the circuit's reactive
# power: with both ladders (reluctivity and strand skin), with the strand skin
# ladder alone, and with both ladders stacked along Beta. CheckEnergy fails a
# case beyond 1e-4; each also pins the energy of its last step.
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 6479 -metiskway ${MPIEXEC_NTASKS} -nooverwrite)
FOREACH(case ladders skin beta)
  FILE(READ sif/${case}.params CASE_PARAMS)
  FILE(WRITE sif/case.params "${CASE_PARAMS}")
  FILE(READ sif/${case}.component CASE_COMPONENT)
  FILE(WRITE sif/case.component "${CASE_COMPONENT}")
  RUN_ELMER_TEST()
ENDFOREACH()
