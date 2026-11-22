!------------------------------------------------------------------------------
!> TRAFOLO (DEV-1667): calls WallHeatTransferCoefficient and WallHeatFlux of the
!> HybridWallLaw library on a face of boundary condition 1, for every wall formula
!> and the three kinds of channel face,
!> with the formula parameters set in turn in that boundary condition, at
!> temperature differences from -5 to 200 K and heights below and above the foot
!> of the wall. Each value must equal the formula written out here, which is the
!> one the case writer used as a MATC expression before the library existed.
!------------------------------------------------------------------------------
SUBROUTINE WallLawCheck(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  REAL(KIND=dp), PARAMETER :: TOLERANCE = 1.0e-12_dp
  REAL(KIND=dp), PARAMETER :: T_AIR = 300.0_dp
  REAL(KIND=dp), PARAMETER :: DIFFERENCES(9) = [ -5.0_dp, 0.0_dp, 0.005_dp, 0.02_dp, 0.3_dp, 4.0_dp, &
      30.0_dp, 75.0_dp, 200.0_dp ]
  REAL(KIND=dp), PARAMETER :: HEIGHTS(3) = [ 0.0_dp, 0.0105_dp, 0.06_dp ]
  ! formula code of each case, then A, B, channel e1, forced hf, outside h0, run start, strip exponent;
  ! the last three are a channel face of an open wall (e1 0), a closed slot (e1 < 0) and a slot in forced flow
  INTEGER, PARAMETER :: N_CASES = 10
  INTEGER, PARAMETER :: CODES(N_CASES) = [ 1, 2, 3, 4, 5, 6, 7, 5, 5, 6 ]
  REAL(KIND=dp), PARAMETER :: PARAMS(7, N_CASES) = RESHAPE([ &
      0.146966_dp, 2.20278_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.146966_dp, 2.20278_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.15_dp, &
      0.31_dp, 1.4_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.01_dp, 0.0_dp, &
      0.31_dp, 1.4_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.01_dp, 0.25_dp, &
      0.0_dp, 1.9_dp, 18.9_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.2_dp, 2.1_dp, 0.0_dp, 25.4_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.2_dp, 2.1_dp, 3.5_dp, 25.4_dp, 0.35_dp, 0.0_dp, 0.0_dp, &
      0.1_dp, 1.9_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.1_dp, 1.9_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.2_dp, 2.1_dp, -1.0_dp, 25.4_dp, 0.0_dp, 0.0_dp, 0.0_dp ], [7, N_CASES])
  REAL(KIND=dp), PARAMETER :: STRIP_DT_K = 50.0_dp
  REAL(KIND=dp) :: Arg(3), h, q, hRef, qRef, d, y
  INTEGER(KIND=AddrInt) :: HtcProc, FluxProc
  TYPE(Element_t), POINTER :: Face, Previous
  TYPE(ValueList_t), POINTER :: BC
  INTEGER :: f, c, law, i, j, Calls, Differ

  HtcProc = GetProcAddr('HybridWallLaw WallHeatTransferCoefficient')
  FluxProc = GetProcAddr('HybridWallLaw WallHeatFlux')
  Face => NULL()
  DO f = 1, GetNOFBoundaryElements()
    Face => GetBoundaryElement(f)
    IF (GetBCId(Face) == 1) EXIT
  END DO
  Previous => SetCurrentElement(Face)
  BC => GetBC(Face)

  Calls = 0
  Differ = 0
  DO c = 1, N_CASES
    law = CODES(c)
    CALL ListAddInteger(BC, 'Hybrid Wall Law', law)
    CALL ListAddConstReal(BC, 'Hybrid Wall Law A', PARAMS(1, c))
    CALL ListAddConstReal(BC, 'Hybrid Wall Law B', PARAMS(2, c))
    CALL ListAddConstReal(BC, 'Hybrid Wall Law Channel', PARAMS(3, c))
    CALL ListAddConstReal(BC, 'Hybrid Wall Law Forced', PARAMS(4, c))
    CALL ListAddConstReal(BC, 'Hybrid Wall Law Outside', PARAMS(5, c))
    CALL ListAddConstReal(BC, 'Hybrid Wall Law Run Start', PARAMS(6, c))
    ! the run-strip case leaves the scale out: the library takes 50 K, as the case writer
    IF (law == 4) THEN
      CALL ListRemove(BC, 'Hybrid Wall Law Strip dT')
    ELSE
      CALL ListAddConstReal(BC, 'Hybrid Wall Law Strip dT', STRIP_DT_K)
    END IF
    CALL ListAddConstReal(BC, 'Hybrid Wall Law Strip Exponent', PARAMS(7, c))
    DO i = 1, SIZE(DIFFERENCES)
      DO j = 1, SIZE(HEIGHTS)
        d = DIFFERENCES(i)
        y = HEIGHTS(j)
        Arg = [ T_AIR + d, T_AIR, y ]
        h = ExecRealFunction(HtcProc, CurrentModel, Face % NodeIndexes(1), Arg)
        q = ExecRealFunction(FluxProc, CurrentModel, Face % NodeIndexes(1), Arg)
        CALL Reference(law, d, y, PARAMS(:, c), hRef, qRef)
        Calls = Calls + 2
        IF (.NOT. Same(h, hRef)) CALL Report('coefficient', h, hRef)
        IF (.NOT. Same(q, qRef)) CALL Report('flux', q, qRef)
      END DO
    END DO
  END DO
  ! a node where Elmer found no reference air temperature passes no heat; a channel face with
  ! A = 1 would give h = 1 and q = 0.02 W/m^2 without the guard
  CALL ListAddInteger(BC, 'Hybrid Wall Law', 5)
  CALL ListAddConstReal(BC, 'Hybrid Wall Law A', 1.0_dp)
  CALL ListAddConstReal(BC, 'Hybrid Wall Law B', 0.0_dp)
  CALL ListAddConstReal(BC, 'Hybrid Wall Law Channel', 0.0_dp)
  Arg = [ T_AIR + 30.0_dp, HUGE(1.0_dp), 0.0_dp ]
  h = ExecRealFunction(HtcProc, CurrentModel, Face % NodeIndexes(1), Arg)
  q = ExecRealFunction(FluxProc, CurrentModel, Face % NodeIndexes(1), Arg)
  Calls = Calls + 2
  d = 30.0_dp
  y = 0.0_dp
  IF (h /= 0.0_dp) CALL Report('coefficient without Tref', h, 0.0_dp)
  IF (q /= 0.0_dp) CALL Report('flux without Tref', q, 0.0_dp)
  WRITE(Message, '(I0,A,I0,A)') Calls, ' wall law calls, ', Differ, ' differ from the formulas'
  CALL Info('WallLawCheck', Message, Level=1)
  Previous => SetCurrentElement(Previous)

CONTAINS

  LOGICAL FUNCTION Same(a, b)
    REAL(KIND=dp), INTENT(IN) :: a, b
    Same = ABS(a - b) <= TOLERANCE * MAX(1.0_dp, ABS(b))
  END FUNCTION Same

  SUBROUTINE Report(what, got, expected)
    CHARACTER(LEN=*), INTENT(IN) :: what
    REAL(KIND=dp), INTENT(IN) :: got, expected
    Differ = Differ + 1
    WRITE(Message, '(A,I0,A,A,A,ES12.4,A,ES12.4,A,ES22.14,A,ES22.14)') 'Hybrid Wall Law ', law, ' ', what, &
        ' at dT = ', d, ' K, y = ', y, ' m: ', got, ' against ', expected
    CALL Info('WallLawCheck', Message, Level=1)
  END SUBROUTINE Report

  ! The MATC expressions of the case writer, spelled out.
  SUBROUTINE Reference(law, diff, y, p, hRef, qRef)
    INTEGER, INTENT(IN) :: law
    REAL(KIND=dp), INTENT(IN) :: diff, y, p(7)
    REAL(KIND=dp), INTENT(OUT) :: hRef, qRef
    REAL(KIND=dp) :: dh, dq, x, r

    dh = MAX(diff, 0.02_dp)
    dq = MAX(diff, 0.0_dp)
    x = MAX(y - p(6), 0.003_dp)
    SELECT CASE (law)
    CASE (1)
      hRef = p(1) + p(2) * dh**0.25_dp
      qRef = (p(1) + p(2) * dq**0.25_dp) * dq
    CASE (2)
      hRef = p(1) + p(2) * Strip(dh) * dh**0.25_dp
      qRef = (p(1) + p(2) * Strip(dq) * dq**0.25_dp) * dq
    CASE (3)
      hRef = p(1) + p(2) * dh**0.25_dp / x**0.25_dp
      qRef = (p(1) + p(2) * dq**0.25_dp / x**0.25_dp) * dq
    CASE (4)
      hRef = p(1) + p(2) * Strip(dh) * dh**0.25_dp / x**0.25_dp
      qRef = (p(1) + p(2) * Strip(dq) * dq**0.25_dp / x**0.25_dp) * dq
    CASE (5)
      r = Channel(dh, p(3))
      hRef = (p(1) + p(2) * dh**0.25_dp) * r
      qRef = hRef * dh
    CASE (6, 7)
      r = Channel(dh, p(3))
      hRef = (((p(1) + p(2) * dh**0.25_dp) * r)**3 + p(4)**3)**(1.0_dp / 3.0_dp)
      IF (law == 7) hRef = hRef + p(5)
      qRef = hRef * dq
    END SELECT
  END SUBROUTINE Reference

  REAL(KIND=dp) FUNCTION Channel(diff, e1)
    REAL(KIND=dp), INTENT(IN) :: diff, e1
    Channel = 1.0_dp
    IF (e1 > 0.0_dp) Channel = 1.0_dp / SQRT(1.0_dp + 200.5_dp / (e1 * diff)**1.5_dp)
    IF (e1 < 0.0_dp) Channel = 0.0_dp
  END FUNCTION Channel

  REAL(KIND=dp) FUNCTION Strip(diff)
    REAL(KIND=dp), INTENT(IN) :: diff
    REAL(KIND=dp) :: s
    s = diff / STRIP_DT_K
    IF (s < 0.5_dp) s = 0.5_dp
    IF (s > 2.0_dp) s = 2.0_dp
    Strip = s**PARAMS(7, c)
  END FUNCTION Strip
!------------------------------------------------------------------------------
END SUBROUTINE WallLawCheck
!------------------------------------------------------------------------------
