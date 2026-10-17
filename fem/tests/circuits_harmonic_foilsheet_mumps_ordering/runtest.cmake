include(test_macros)
RUN_ELMER_EXPECT(SIF sif/metis.sif EXIT ZERO
  MATCH "ZMumps_SolveSystem: MUMPS ordering used [(]INFOG[(]7[)][)]: 5,"
        "CompareToReferenceSolution: PASSED all")
RUN_ELMER_EXPECT(SIF sif/qamd.sif EXIT ZERO
  MATCH "ZMumps_SolveSystem: MUMPS ordering used [(]INFOG[(]7[)][)]: 6,"
        "CompareToReferenceSolution: PASSED all")
RUN_ELMER_EXPECT(SIF sif/unknown.sif EXIT NONZERO
  MATCH "Unknown \"Mumps Ordering\": nested [(]amd, amf, scotch, pord, metis, qamd or auto[)]")
