! DEV-1549: the model of the skin ladder safeguard cases (mesh 2241, a stranded
! coil with Transient Homogenization and 1 A DC); each case file adds its
! Simulation and Component sections.
Constants 1
   Permittivity of Vacuum = 8.8542e-12
End
Initial Condition 1
   A {e} = Real 0
End
Solver 1
   Exec Solver = Before all
   Procedure = "DirectionSolver" "DirectionSolver"
   Linear System Solver = Iterative
   Linear System Iterative Method = BiCGStab
   Linear System Max Iterations = 5000
   Linear System Convergence Tolerance = 1.0e-10
   Linear System Abort Not Converged = True
   Linear System Residual Output = 1000
   Equation = Direction Alpha
   Variable =  Alpha
End
Solver 2
   Exec Solver = Before all
   Procedure = "DirectionSolver" "DirectionSolver"
   Linear System Solver = Iterative
   Linear System Iterative Method = BiCGStab
   Linear System Max Iterations = 5000
   Linear System Convergence Tolerance = 1.0e-10
   Linear System Abort Not Converged = True
   Linear System Residual Output = 1000
   Equation = Direction Beta
   Variable =  Beta
End
Solver 3
   Exec Solver = Before All
   Equation = Transformation matrix
   Procedure = "CoordinateTransform" "RotMSolver"
   Variable = -nooutput Dummy
   Optimize Bandwidth = False
   Polar Decomposition Determinant Tolerance = Real 1.0e-9
   Polar Decomposition Max Iterations = Integer 100
   Exported Variable 1 = -nooutput RotM E[RotM E:9]
   Exported Variable 2 = -nooutput Alpha Vector E[Alpha Vector E:3]
   Exported Variable 3 = -nooutput Beta Vector E[Beta Vector E:3]
   Exported Variable 4 = -nooutput Gamma Vector E[Gamma Vector E:3]
   Discontinuous Galerkin = True
End
Solver 4 !-- WpotentialSolver - Wsolve
   Exec Solver = Before All
   Procedure = "WPotentialSolver" "Wsolve"
   Equation = "Wire direction"
   Variable = W
   Linear System Solver = Iterative
   Linear System Iterative Method = CG
   Linear System Max Iterations = 10000
   Linear System Convergence Tolerance = 1.0e-10
   Linear System Abort Not Converged = True
   Linear System Residual Output = 1000
End
Solver 5 !-- CircuitsAndDynamics - CircuitsAndDynamicsHarmonic
   Exec Solver = Always
   Equation = Circuits
   Procedure = "CircuitsAndDynamics" "CircuitsAndDynamics"
End
Solver 6  !---- MagnetoDynamics, WhitneyAVSolver
   Exec Solver = Always
   Equation = "MGDynamics"
   Variable = "A"
   Procedure = "MagnetoDynamics" "WhitneyAVSolver"
   NonLinear System Relaxation Factor=1
   Export Lagrange Multiplier = Logical True
   Exported Variable 1 = -elem Xi Alpha
   Exported Variable 2 = -elem Xi Beta
   Exported Variable 3 = -elem Proximity Loss
   ! GCR, not BiCGStab(l), which breaks down (0/0) in the second step of these
   ! short steps.
   Linear System Solver = Iterative
   Linear System Iterative Method = GCR
   Linear System GCR Restart = 200
   Linear System preconditioning = none
   Linear System Convergence Tolerance = 1.e-11
   Linear System Max Iterations = 3000
   Linear System Residual Output = 1000
   Steady State Convergence Tolerance = 1e-08
End
Solver 7
   Exec Solver = Always
   Equation = "MGDynamicsCalc"
   Procedure = "MagnetoDynamics" "MagnetoDynamicsCalcFields"
   Linear System Symmetric = True
   Potential Variable = String "A"
   Calculate Current Density = Logical True
   Calculate Joule Heating = Logical True
   Calculate Nodal Fields = Logical False
   Loss Estimation = Logical True
   Steady State Convergence Tolerance = 0
   Linear System Solver = "Iterative"
   Linear System Preconditioning = None
   Linear System Residual Output = 1000
   Linear System Max Iterations = 5000
   Linear System Iterative Method = CG
   Linear System Convergence Tolerance = 1.0e-10
End
Solver 8
   Exec Solver = Always
   Equation = Circuits Output
   Procedure = "CircuitsAndDynamics" "CircuitsOutput"
End
Equation 1
   Active Solvers(2) = 6 7
End
Equation 2
   Active Solvers(7) = 1 2 3 4 5 6 7
End
Material 1
   Name = iron
   Electric Conductivity = 0
   Relative Permeability = 2000
   Relative Permittivity = 1e100
End
Material 2
   Name = air
   Electric Conductivity = 0
   Relative Permeability = 1
   Relative Permittivity = 1
End
Material 3
   Name = metal
   Relative Permeability = 1
   Relative Permittivity = 1e100
   Electric Conductivity = 1e7
End
Body 1
   Name = core
   Target Bodies(1) = $ core
   Equation = 1
   Material = 1
   Initial Condition = 1
End
Body 2
   Name = air
   Target Bodies(1) = $ air
   Equation = 1
   Material = 2
   Initial Condition = 1
End
Body 3
   Name = airgap_1_limb_2
   Target Bodies(1) = $ airgap_1_limb_2
   Equation = 1
   Material = 2
   Initial Condition = 1
End
Body 4
   Name = wp1
   Target Bodies(1) = $ wp1
   Equation = 2
   Material = 3
   Initial Condition = 1
   Body Force = 1
   Alpha reference (3) = Real 1 0 0
   Beta reference (3) = Real 0 1 0
End
Boundary Condition 1
   Name = BCn Flux Parallel
   Target Boundaries(2) = $ coreface_xy xy0
   A {e} = Real 0
End
Boundary Condition 2
   Name = Alpha0
   Target Boundaries = $ wp1_alpha0
   Body 4: Alpha = Real 0
End
Boundary Condition 3
   Name = Alpha1
   Target Boundaries = $ wp1_alpha1
   Body 4: Alpha = Real 1
End
Boundary Condition 4
   Name = ground
   Target Boundaries = $ wp1_gamma1
   W = Real 0
   A {e} = Real 0
End
Boundary Condition 5
   Name = current in foil winding
   Target Boundaries = $ wp1_gamma0
   W = Real 1
   A {e} = Real 0
End
Boundary Condition 6
   Name = Beta1
   Target Boundaries = $ wp1_beta1
   Body 4: Beta = Real 1
End
Boundary Condition 7
   Name = Beta0
   Target Boundaries = $ wp1_beta0
   Body 4: Beta = Real 0
End
Body Force 1
   Name = "Circuit"
   ! i_component(1) = -S = 1 A
   S = Real -1.0
End
