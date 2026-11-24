!------------------------------------------------------------------------------
!> TRAFOLO (DEV-1667): calls FilmConductivity of the HybridWallLaw library at the
!> nodes of one element of each body, after setting the wall distance and the
!> air velocity at those nodes, and compares it with the law written out here.
!> Body 1 has the film of forced flow (it fades with the fastest node of the
!> element and the cells at a wall conduct at least like the wall formula),
!> body 2 the film that does not fade.
!------------------------------------------------------------------------------
SUBROUTINE FilmCheck(Model, Solver, dt, Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  REAL(KIND=dp), PARAMETER :: TOLERANCE = 1.0e-12_dp
  REAL(KIND=dp), PARAMETER :: AIR_K = 0.0262_dp, FILM_K = 7.86_dp, THICKNESS = 0.006_dp, &
      FILM_SPEED = 0.5_dp, WALL_COEFFICIENT = 64.0_dp
  INTEGER, PARAMETER :: N_CASES = 8
  ! body, wall distance of the first four and of the last four nodes [m], speed of the fastest node [m/s]:
  ! still air at a wall, moving air at a wall, half the film speed, beyond the film, where the film fades,
  ! a wall cell thicker than the film; then the film that does not fade, inside and beyond its thickness
  INTEGER, PARAMETER :: BODIES(N_CASES) = [ 1, 1, 1, 1, 1, 1, 2, 2 ]
  REAL(KIND=dp), PARAMETER :: CASES(3, N_CASES) = RESHAPE([ &
      0.0_dp, 0.0025_dp, 0.0_dp, &
      0.0_dp, 0.0025_dp, 1.0_dp, &
      0.0_dp, 0.0025_dp, 0.25_dp, &
      0.01_dp, 0.0125_dp, 0.0_dp, &
      0.0069_dp, 0.0074_dp, 0.0_dp, &
      0.0_dp, 0.02_dp, 2.0_dp, &
      0.003_dp, 0.0069_dp, 5.0_dp, &
      0.02_dp, 0.0225_dp, 5.0_dp ], [3, N_CASES])
  TYPE(Variable_t), POINTER :: Distance, Flow
  TYPE(Element_t), POINTER :: Element, Previous
  INTEGER(KIND=AddrInt) :: FilmProc
  REAL(KIND=dp) :: Arg(1), d(8), k, kRef, Fastest
  INTEGER :: c, e, i, n, node, Calls, Differ

  Distance => VariableGet(Solver % Mesh % Variables, 'Wall Distance', UnfoundFatal=.TRUE.)
  Flow => VariableGet(Solver % Mesh % Variables, 'Flow Solution', UnfoundFatal=.TRUE.)
  FilmProc = GetProcAddr('HybridWallLaw FilmConductivity')

  Calls = 0
  Differ = 0
  DO c = 1, N_CASES
    Element => NULL()
    DO e = 1, GetNOFActive()
      Element => GetActiveElement(e)
      IF (Element % BodyId == BODIES(c)) EXIT
    END DO
    n = GetElementNOFNodes(Element)
    IF (n /= 8) CALL Fatal('FilmCheck', 'Eight-node elements expected.')
    Previous => SetCurrentElement(Element)
    d(1:4) = CASES(1, c)
    d(5:8) = CASES(2, c)
    Fastest = CASES(3, c)
    Flow % Values = 0.0_dp
    DO i = 1, n
      node = Element % NodeIndexes(i)
      Distance % Values(Distance % Perm(node)) = d(i)
    END DO
    node = Element % NodeIndexes(n)
    Flow % Values(Flow % DOFs*(Flow % Perm(node)-1)+2) = Fastest
    DO i = 1, n
      Arg(1) = d(i)
      k = ExecRealFunction(FilmProc, CurrentModel, Element % NodeIndexes(i), Arg)
      kRef = Reference(BODIES(c), d(i), MINVAL(d), MAXVAL(d), Fastest)
      Calls = Calls + 1
      IF (.NOT. (ABS(k - kRef) <= TOLERANCE * MAX(1.0_dp, ABS(kRef)))) THEN
        Differ = Differ + 1
        WRITE(Message, '(A,I0,A,ES12.4,A,ES22.14,A,ES22.14)') 'case ', c, ' at wall distance ', d(i), &
            ' m: ', k, ' against ', kRef
        CALL Info('FilmCheck', Message, Level=1)
      END IF
    END DO
    Previous => SetCurrentElement(Previous)
  END DO
  WRITE(Message, '(I0,A,I0,A)') Calls, ' film calls, ', Differ, ' differ from the law'
  CALL Info('FilmCheck', Message, Level=1)

CONTAINS

  REAL(KIND=dp) FUNCTION Reference(body, dist, nearest, farthest, speed)
    INTEGER, INTENT(IN) :: body
    REAL(KIND=dp), INTENT(IN) :: dist, nearest, farthest, speed
    REAL(KIND=dp) :: still, inside

    still = 1.0_dp
    IF (body == 1) still = MIN(MAX(1.0_dp - speed / FILM_SPEED, 0.0_dp), 1.0_dp)
    inside = MIN(MAX((1.3_dp * THICKNESS - dist) / (0.3_dp * THICKNESS), 0.0_dp), 1.0_dp)
    Reference = AIR_K + (FILM_K - AIR_K) * still * inside
    IF (body == 1 .AND. nearest <= 0.0_dp) Reference = MAX(Reference, WALL_COEFFICIENT * farthest)
  END FUNCTION Reference
!------------------------------------------------------------------------------
END SUBROUTINE FilmCheck
!------------------------------------------------------------------------------
