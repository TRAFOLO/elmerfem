! Check every local node, including shared nodes, on every MPI rank.
SUBROUTINE CheckBounds(Model,Solver,dt,Transient)
  USE DefUtils
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  TYPE(Variable_t), POINTER :: Temperature
  TYPE(Mesh_t), POINTER :: Mesh
  REAL(KIND=dp) :: NodeTemp
  REAL(KIND=dp), PARAMETER :: BoundsTolerance = 1.0e-8_dp
  INTEGER :: i,p

  Mesh => Solver % Mesh
  Temperature => VariableGet(Mesh % Variables,'Temperature',UnfoundFatal=.TRUE.)
  DO i=1,Mesh % NumberOfNodes
    p = Temperature % Perm(i)
    IF ( p <= 0 ) CALL Fatal('CheckBounds','Missing temperature node')
    NodeTemp = Temperature % Values(p)
    IF ( .NOT. IEEE_IS_FINITE(NodeTemp) ) &
        CALL Fatal('CheckBounds','Nonfinite temperature')
    IF ( NodeTemp < 300.0_dp-BoundsTolerance .OR. NodeTemp > 310.0_dp+BoundsTolerance ) &
        CALL Fatal('CheckBounds','Temperature outside [300,310]')
  END DO
END SUBROUTINE CheckBounds
