include(test_macros)
RUN_ELMER_EXPECT(SIF sif/max2.sif EXIT ZERO MATCH "NOT CONVERGED: nonlinear solver=\"mgdynamics\" iterations=2 change=")
RUN_ELMER_EXPECT(SIF sif/single.sif EXIT ZERO NOMATCH "NOT CONVERGED: nonlinear")
RUN_ELMER_EXPECT(SIF sif/relaxed.sif EXIT ZERO MATCH "NOT CONVERGED: nonlinear solver=\"mgdynamics\" iterations=1 relaxation= 5.000E-01")
