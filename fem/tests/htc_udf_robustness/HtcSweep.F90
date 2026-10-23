!------------------------------------------------------------------------------
!> TRAFOLO (DEV-1622): calls htc_evaluate of htc_udf
!>  - for one face of each boundary condition at surface temperatures from -HUGE
!>    to +HUGE C, NaN and +-Infinity,
!>  - on the face of boundary condition 1 with every convection parameter in turn
!>    set to NaN, +-Infinity, +-HUGE, -1, 0, 1e-300 or 1e300 (values a SIF cannot
!>    hold), for every orientation code and two fluids.
!> Every heat transfer coefficient must be finite and >= 0, and the physical face
!> of boundary condition 1 must give the correlation's value.
!------------------------------------------------------------------------------
SUBROUTINE HtcSweep(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN, IEEE_POSITIVE_INF, &
                                           IEEE_NEGATIVE_INF, IEEE_IS_FINITE
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  ! Vertical face, Lc_nat 0.01 m, still air at 20 C, surface at 100 C: Churchill-Chu with the
  ! air fits of htc_udf, computed independently (Ra 4641, Nu 4.9185).
  REAL(KIND=dp), PARAMETER :: H_REFERENCE = 14.130427242162584_dp
  REAL(KIND=dp), PARAMETER :: H_REFERENCE_TOLERANCE = 1.0e-10_dp
  REAL(KIND=dp), PARAMETER :: T_REFERENCE = 100.0_dp
  INTEGER, PARAMETER :: MAX_BCS = 16
  CHARACTER(LEN=*), PARAMETER :: REAL_PARAMETERS(7) = [ 'Lc_nat        ', 'Lc_forced     ', &
      'theta_deg     ', 'velocity_mag  ', 'channel_gap   ', 'channel_length', 'T_ambient     ' ]
  ! A valid value of each real parameter; the sweep spoils one at a time.
  REAL(KIND=dp), PARAMETER :: VALID_VALUES(7) = [ 0.01_dp, 0.05_dp, 30.0_dp, 2.0_dp, 0.005_dp, 0.05_dp, 20.0_dp ]
  INTEGER, PARAMETER :: NAT_ORIENTS(7) = [ 0, 1, 2, 3, 4, 5, 99 ]
  INTEGER, PARAMETER :: FCD_ORIENTS(5) = [ 0, 1, 2, 3, 77 ]
  INTEGER, PARAMETER :: FLUIDS(2) = [ 1, 42 ]
  REAL(KIND=dp) :: Temperatures(17), Spoilt(9), Arg(1), h, hReference, Nan, Inf
  INTEGER(KIND=AddrInt) :: Htc
  TYPE(Element_t), POINTER :: Face, Previous, Face1
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Done(MAX_BCS)
  INTEGER :: f, k, BcId, Calls, Bad, ParamCalls, ParamBad, iNat, iFcd, iFluid, iPar, iVal, iT
!------------------------------------------------------------------------------
  Nan = IEEE_VALUE(1.0_dp, IEEE_QUIET_NAN)
  Inf = IEEE_VALUE(1.0_dp, IEEE_POSITIVE_INF)
  Temperatures = [ -HUGE(1.0_dp), -1.0e300_dp, -1.0e6_dp, -273.15_dp, -100.0_dp, 0.0_dp, 20.0_dp, &
                   T_REFERENCE, 1226.0_dp, 1227.0_dp, 5000.0_dp, 1.0e6_dp, 1.0e300_dp, HUGE(1.0_dp), &
                   Nan, Inf, -Inf ]
  Spoilt = [ Nan, Inf, -Inf, -HUGE(1.0_dp), -1.0_dp, 0.0_dp, 1.0e-300_dp, 1.0e300_dp, HUGE(1.0_dp) ]

  Htc = GetProcAddr('htc_udf htc_evaluate')
  Done = .FALSE.
  Calls = 0
  Bad = 0
  hReference = -1.0_dp
  Face1 => NULL()
  DO f = 1, GetNOFBoundaryElements()
    Face => GetBoundaryElement(f)
    BcId = GetBCId(Face)
    IF (BcId < 1 .OR. BcId > MAX_BCS) CYCLE
    IF (Done(BcId)) CYCLE
    Done(BcId) = .TRUE.
    IF (BcId == 1) Face1 => Face
    Previous => SetCurrentElement(Face)
    DO k = 1, SIZE(Temperatures)
      Arg(1) = Temperatures(k)
      h = ExecRealFunction(Htc, CurrentModel, Face % NodeIndexes(1), Arg)
      Calls = Calls + 1
      IF (.NOT. Valid(h)) THEN
        Bad = Bad + 1
        WRITE(Message, '(A,I0,A,ES12.4E3,A,ES12.4E3)') 'Boundary condition ', BcId, ': T = ', Arg(1), &
            ' C gives h = ', h
        CALL Info('HtcSweep', Message, Level=1)
      END IF
      IF (BcId == 1 .AND. k == 8) hReference = h
    END DO
  END DO
  WRITE(Message, '(I0,A,I0,A)') Calls, ' face calls, ', Bad, ' not finite or negative'
  CALL Info('HtcSweep', Message, Level=1)
  WRITE(Message, '(A,ES23.16,A,ES23.16)') 'Boundary condition 1 at 100 C: h = ', hReference, &
      ', expected ', H_REFERENCE
  CALL Info('HtcSweep', Message, Level=1)
  IF (ABS(hReference - H_REFERENCE) <= H_REFERENCE_TOLERANCE * H_REFERENCE) &
    CALL Info('HtcSweep', 'Reference h matches', Level=1)

  ! The parameters of boundary condition 1, spoilt in memory one at a time.
  ParamCalls = 0
  ParamBad = 0
  Params => CurrentModel % BCs(1) % Values
  Previous => SetCurrentElement(Face1)
  DO iNat = 1, SIZE(NAT_ORIENTS)
    DO iFcd = 1, SIZE(FCD_ORIENTS)
      DO iFluid = 1, SIZE(FLUIDS)
        CALL ListAddInteger(Params, 'NatOrient', NAT_ORIENTS(iNat))
        CALL ListAddInteger(Params, 'FcdOrient', FCD_ORIENTS(iFcd))
        CALL ListAddInteger(Params, 'Fluid Id', FLUIDS(iFluid))
        DO iPar = 1, SIZE(REAL_PARAMETERS)
          DO iVal = 1, SIZE(Spoilt)
            DO k = 1, SIZE(REAL_PARAMETERS)
              CALL ListAddConstReal(Params, TRIM(REAL_PARAMETERS(k)), VALID_VALUES(k))
            END DO
            CALL ListAddConstReal(Params, TRIM(REAL_PARAMETERS(iPar)), Spoilt(iVal))
            DO iT = 1, SIZE(Temperatures)
              Arg(1) = Temperatures(iT)
              h = ExecRealFunction(Htc, CurrentModel, Face1 % NodeIndexes(1), Arg)
              ParamCalls = ParamCalls + 1
              IF (.NOT. Valid(h)) THEN
                ParamBad = ParamBad + 1
                WRITE(Message, '(A,I0,A,I0,A,I0,3A,ES12.4E3,A,ES12.4E3,A,ES12.4E3)') 'NatOrient ', &
                    NAT_ORIENTS(iNat), ' FcdOrient ', FCD_ORIENTS(iFcd), ' Fluid ', FLUIDS(iFluid), ' ', &
                    TRIM(REAL_PARAMETERS(iPar)), ' = ', Spoilt(iVal), ', T = ', Arg(1), ' C gives h = ', h
                CALL Info('HtcSweep', Message, Level=1)
              END IF
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO
  WRITE(Message, '(I0,A,I0,A)') ParamCalls, ' parameter calls, ', ParamBad, ' not finite or negative'
  CALL Info('HtcSweep', Message, Level=1)

CONTAINS

  LOGICAL FUNCTION Valid(x)
    REAL(KIND=dp), INTENT(IN) :: x
    Valid = IEEE_IS_FINITE(x) .AND. x >= 0.0_dp
  END FUNCTION Valid
!------------------------------------------------------------------------------
END SUBROUTINE HtcSweep
!------------------------------------------------------------------------------
