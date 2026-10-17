!------------------------------------------------------------------------------
!> mumps_refinement_row_scale: the row scales and the componentwise backward
!> error of the MUMPS refinement (MumpsRowScale, MumpsBackwardError in
!> DirectSolve) on fixed numbers, without a factorization.
!------------------------------------------------------------------------------
SUBROUTINE MumpsRowScaleTest( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE DirectSolve
  USE, INTRINSIC :: IEEE_ARITHMETIC
  IMPLICIT NONE
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  INTEGER, PARAMETER :: n = 7
  REAL(KIND=dp) :: Activity(n), Scale(n), r(n)
  INTEGER :: GBlock(n)

  ! block 1: a row with a real activity, an empty row, a row exactly at the threshold;
  ! block 3: a circuit row with a real activity and an empty one; block 2: all zero
  Activity = [ 1.0_dp, 1.0e-22_dp, 1.0e-12_dp, 5.0_dp, 1.0e-21_dp, 0.0_dp, 0.0_dp ]
  GBlock = [ 1, 1, 1, 3, 3, 2, 2 ]
  IF ( .NOT. MumpsRowScale( Activity, GBlock, Scale ) ) CALL Fatal( 'MumpsRowScaleTest', 'scales not valid' )

  ! rounding noise as large as the activity of the empty rows: measured against the block
  r = 0.0_dp
  r(2) = 1.0e-22_dp
  r(5) = 1.0e-21_dp
  CALL Check( 'noise in empty rows', MumpsBackwardError( r, Scale ), 2.0e-22_dp )

  ! a residual in a row with a real activity is measured against that activity
  r = 0.0_dp
  r(1) = 1.0e-10_dp
  CALL Check( 'residual of an active row', MumpsBackwardError( r, Scale ), 1.0e-10_dp )

  ! a row exactly at the threshold is still measured componentwise
  r = 0.0_dp
  r(3) = 1.0e-23_dp
  CALL Check( 'row at the threshold', MumpsBackwardError( r, Scale ), 1.0e-11_dp )

  ! a zero block keeps zero scales: a zero residual counts nothing, any other is infinite
  r = 0.0_dp
  CALL Check( 'zero residual', MumpsBackwardError( r, Scale ), 0.0_dp )
  r(6) = 1.0e-30_dp
  IF ( MumpsBackwardError( r, Scale ) < HUGE(1.0_dp) ) &
      CALL Fatal( 'MumpsRowScaleTest', 'a residual in a zero row must give HUGE' )

  ! a non-finite activity makes the scales invalid
  Activity(4) = IEEE_VALUE( Activity(4), IEEE_POSITIVE_INF )
  IF ( MumpsRowScale( Activity, GBlock, Scale ) ) &
      CALL Fatal( 'MumpsRowScaleTest', 'a non-finite activity must invalidate the scales' )

  CALL Info( 'MumpsRowScaleTest', 'all row scale checks passed' )

CONTAINS

  SUBROUTINE Check( What, omega, expected )
    CHARACTER(LEN=*) :: What
    REAL(KIND=dp) :: omega, expected

    WRITE( Message, '(A,A,ES12.4,A,ES12.4)' ) What, ': omega', omega, ' expected', expected
    CALL Info( 'MumpsRowScaleTest', Message )
    IF ( ABS( omega - expected ) > 1.0e-6_dp * expected ) CALL Fatal( 'MumpsRowScaleTest', 'wrong omega: '//What )
  END SUBROUTINE Check
!------------------------------------------------------------------------------
END SUBROUTINE MumpsRowScaleTest
!------------------------------------------------------------------------------
