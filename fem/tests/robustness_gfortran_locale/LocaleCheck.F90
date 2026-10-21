!------------------------------------------------------------------------------
!> TRAFOLO: formatted writes with ordinary allocation traffic, checking after
!> each that the C runtime's numeric locale is still "C" (GfortranLocale.c:
!> without the workaround it is the user's locale after a few percent of them).
!------------------------------------------------------------------------------
SUBROUTINE LocaleCheck(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  INTERFACE
    FUNCTION CSetlocale(Category, Locale) RESULT(Name) BIND(C, NAME='setlocale')
      IMPORT :: C_INT, C_PTR
      INTEGER(C_INT), VALUE :: Category
      TYPE(C_PTR), VALUE :: Locale
      TYPE(C_PTR) :: Name
    END FUNCTION CSetlocale
  END INTERFACE
  INTEGER(C_INT), PARAMETER :: LC_NUMERIC = 4 ! of the UCRT
  INTEGER, PARAMETER :: Writes = 200000
  CHARACTER(KIND=C_CHAR, LEN=2), TARGET :: CName = 'C'//C_NULL_CHAR
  CHARACTER(LEN=80) :: Msg
  REAL(KIND=dp), ALLOCATABLE :: a(:)
  TYPE(C_PTR) :: Name
  INTEGER :: i, Changed
!------------------------------------------------------------------------------
  Changed = 0
  DO i = 1, Writes
    ALLOCATE(a(1 + MOD(i, 7) * 3))
    a = 0.0_dp
    WRITE(Msg, '(A,G15.8,G15.8)') 'NS (NRM,RELC): (', 1.0_dp / i, SUM(a)
    DEALLOCATE(a)
    IF (.NOT. IsC(CSetlocale(LC_NUMERIC, C_NULL_PTR))) THEN
      Changed = Changed + 1
      Name = CSetlocale(LC_NUMERIC, C_LOC(CName))
    END IF
  END DO
  CALL Info('LocaleCheck', 'Numeric locale changed after '//I2S(Changed)//' of '// &
      I2S(Writes)//' formatted writes', Level=3)

CONTAINS

  LOGICAL FUNCTION IsC(p)
    TYPE(C_PTR), INTENT(IN) :: p
    CHARACTER(KIND=C_CHAR), POINTER :: s(:)
    IsC = .FALSE.
    IF (.NOT. C_ASSOCIATED(p)) RETURN
    CALL C_F_POINTER(p, s, [2])
    IsC = s(1) == 'C' .AND. s(2) == C_NULL_CHAR
  END FUNCTION IsC
!------------------------------------------------------------------------------
END SUBROUTINE LocaleCheck
!------------------------------------------------------------------------------
