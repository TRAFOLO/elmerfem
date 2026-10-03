!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Module for computing the directions of a coil
! *
! *  Author: Eelis Takala
! *  Email:   eelis.takala@trafotek.fi
! *  Web:     http://www.trafotek.fi
! *  Address: Trafotek
! *           Kaarinantie 700
! *           20540 Turku
! *
! *  Original Date: December 2015
! *
! *****************************************************************************/
 
!> \ingroup Solvers
!> \{
!------------------------------------------------------------------------------
SUBROUTINE DirectionSolver_Init0(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
  CHARACTER(LEN=MAX_NAME_LEN) :: varname
  LOGICAL :: Found
  INTEGER, POINTER :: Active(:)
  INTEGER :: mysolver,i,j,k,l,n,m,vDOFs, soln
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Solver_t), POINTER :: Solvers(:), PSolver
  INTEGER, SAVE :: visited=0

  visited = visited + 1

  ! This is really using DG so we don't need to make any dirty tricks to create DG fields
  ! as is done in this initialization. 
  IF (GetLogical(GetSolverParams(),'Discontinuous Galerkin',Found)) RETURN

  PSolver => Solver
  DO mysolver=1,Model % NumberOfSolvers
    IF ( ASSOCIATED(PSolver,Model % Solvers(mysolver)) ) EXIT
  END DO

  varname = GetString(GetSolverParams(), 'Variable', Found)

  n = Model % NumberOfSolvers
  DO i=1,Model % NumberOFEquations
    Active => ListGetIntegerArray(Model % Equations(i) % Values, &
                'Active Solvers', Found)
    m = SIZE(Active)
    IF ( ANY(Active==mysolver) ) &
      CALL ListAddIntegerArray( Model % Equations(i) % Values,  &
           'Active Solvers', m+1, [Active, n+1] )
  END DO

  ! Create DG solver structures on-the-fly without actually solving the matrix
  ! equations. It is assumed that the DG field within each element is independent
  ! and hence no coupling between elemental fields is needed. 
  ALLOCATE(Solvers(n+1))
  Solvers(1:n) = Model % Solvers
  Solvers(n+1) % Values => ListAllocate()
  SolverParams => Solvers(n+1) % Values
  CALL ListAddLogical( SolverParams, 'Discontinuous Galerkin', .TRUE. )
  Solvers(n+1) % DG = .TRUE.
  Solvers(n+1) % PROCEDURE = 0
  Solvers(n+1) % ActiveElements => NULL()
  CALL ListAddString( SolverParams, 'Exec Solver', 'never' )
  CALL ListAddLogical( SolverParams, 'No Matrix',.TRUE.)
  CALL ListAddLogical( SolverParams, 'Optimize Bandwidth',.FALSE.)
  CALL ListAddString( SolverParams, 'Equation', &
  'elementaladd'//TRIM(varname) )
  CALL ListAddString( SolverParams, 'Procedure', &
              'DirectionSolver DirectionSolver_Dummy',.FALSE. )
  CALL ListAddString( SolverParams, 'Variable', '-nooutput '//TRIM(varname)//'_dummy' )


!  pname = ListGetString( Model % Solvers(soln) % Values, 'Mesh', Found )
!  IF(Found) THEN
!    CALL ListAddString( SolverParams, 'Mesh', pname )
!  END IF

  i = 1
  DO WHILE(.TRUE.)
    IF(ListCheckPresent(SolverParams, "Exported Variable "//i2s(i))) THEN
      i=i+1
    ELSE
      EXIT
    END IF
  END DO

  CALL ListAddString( SolverParams, "Exported Variable "//i2s(i), &
        TRIM(varname)//" Direction" )

  DEALLOCATE(Model % Solvers)
  Model % Solvers => Solvers
  Model % NumberOfSolvers = n+1
!------------------------------------------------------------------------------
END SUBROUTINE DirectionSolver_Init0
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE DirectionSolver_Dummy(Model,Solver,dt,Transient)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
END SUBROUTINE DirectionSolver_Dummy
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE DirectionSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the Poisson equation!
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: AllocationsDone = .FALSE., Found, PosEl, NegEl
  TYPE(Element_t), POINTER :: Element

  CHARACTER(LEN=MAX_NAME_LEN) :: varname, Namespace, VNWithNS, DirMethod
  REAL(KIND=dp) :: Norm
  INTEGER :: n, nb, nd, t, istat, active, NofNameSpaces, ns_iter
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: BodyForce, BC
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:)
  LOGICAL, ALLOCATABLE :: ActiveBody(:)

  ! Boundary values closer than this are the same value, and a wall thinner than
  ! this times the bounding box diagonal has no direction.
  REAL(KIND=dp), PARAMETER :: BC_VALUE_TOL = 1.0e-12_dp
  REAL(KIND=dp), PARAMETER :: REL_THICKNESS_TOL = 1.0e-14_dp
  ! Largest difference, relative to the widest body span, between the boundary
  ! values two touching bodies give a node on faces of both. In a chain they are
  ! equal, so a larger difference contradicts the numbering.
  REAL(KIND=dp), PARAMETER :: REL_SHARED_VALUE_TOL = 1.0e-2_dp

  SAVE STIFF, LOAD, FORCE, AllocationsDone
!------------------------------------------------------------------------------

  varname = GetString(GetSolverParams(), 'Variable', Found)
  IF (.NOT.ASSOCIATED(Solver % Matrix)) RETURN

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()

  DirMethod = GetString(GetSolverParams(), 'Direction Method', Found)
  IF (.NOT. Found) DirMethod = 'laplace'

  SELECT CASE(DirMethod)
  CASE('laplace')
    CONTINUE
  CASE('distance')
    CALL DirectionByDistance()
    RETURN
  CASE DEFAULT
    CALL Fatal('DirectionSolver','Unknown Direction Method: '//TRIM(DirMethod))
  END SELECT

  IF ( .NOT. AllocationsDone ) THEN
     N = Solver % Mesh % MaxElementDOFs  ! just big enough for elemental arrays
     ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'DirectionSolver', 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
  END IF
  NofNameSpaces = Model%numberofbodies 
  DO ns_iter=1,NofNameSpaces 
   Namespace='body '//i2s(ns_iter)//':'
   VNWithNS = TRIM(Namespace)//' '//TRIM(varname)
   IF (ListCheckPresentAnyBC(Model,VNWithNS)) CALL ListSetNameSpace(Namespace)
   !System assembly:
   !----------------
   Active = GetNOFActive()
   CALL DefaultInitialize()
   DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()

      LOAD = 0.0d0
      BodyForce => GetBodyForce()
      IF ( ASSOCIATED(BodyForce) ) &
         Load(1:n) = GetReal( BodyForce, 'Source', Found )

      !Get element local matrix and rhs vector:
      !----------------------------------------
      CALL LocalMatrix(  STIFF, FORCE, LOAD, Element, n, nd+nb )
      CALL CondensateP( nd, nb, STIFF, FORCE )
      CALL DefaultUpdateEquations( STIFF, FORCE )
   END DO

   DO t=1, Solver % Mesh % NumberOfBoundaryElements
     ! get element and BC info
     ! -----------------------
     Element => GetBoundaryElement(t)
     IF ( .NOT.ActiveBoundaryElement() ) CYCLE
     n = GetElementNOFNodes()
     ! no evaluation of Neumann BC’s on points
     IF ( GetElementFamily() == 1 ) CYCLE
     BC => GetBC()
     FORCE = 0.0d00
     STIFF = 0.0d00

     PosEl = .False.
     NegEl = .False.
     PosEl = GetLogical(BC,'Positive Electrode', Found)
     NegEl = GetLogical(BC,'Negative Electrode', Found)

     LOAD = 0._dp
     IF (PosEl) LOAD = 1._dp
     IF (NegEl) LOAD = -1._dp
     
     IF (Solver % Variable % name == 'w') CALL BoundaryCondition(LOAD, FORCE, Element, n)

     CALL DefaultUpdateEquations( STIFF, FORCE )
   END DO
   CALL DefaultFinishAssembly()
   CALL DefaultDirichletBCs()
   ! And finally, solve:
   !--------------------
   Norm = DefaultSolve()
   CALL SaveSolutionWithBodyMethod(ns_iter)
   CALL ListSetNameSpace('')
!------------------------------------------------------------------------------
  END DO ! namespaces
!------------------------------------------------------------------------------
CONTAINS

!------------------------------------------------------------------------------
!> Set the direction variable to the relative distance between the two Dirichlet
!> faces of its body: value = v_lo + (v_hi-v_lo)*d_lo/(d_lo+d_hi). This gives
!> |grad var| = 1/T for a wall of constant thickness T of any shape, whereas the
!> Laplace solution of the default method behaves like ln(r) in a round coil.
!> Faces and nodes are grouped by body, so a model with several windings measures
!> each node against the faces of its own winding only, and each body takes its
!> own pair of values from the faces it is a parent of.
!> Windings that touch share a face and its nodes. Numbered as a chain, winding k
!> spanning k-1..k (any integer shift), the shared face carries one value that
!> ends the pairs of both windings. A shared node gets the value of each of its
!> bodies: the elemental copy of the field (the 'Direction' variable) holds each
!> body's own value, the nodal field that of the first body. On a face with a
!> boundary value of both bodies the values must agree, it is fatal otherwise;
!> across a face free for this field, as between windings of a different height
!> or build, each body keeps its own. The 'body N:' namespace trick of the
!> Laplace path is not supported here: the boundary values are read without a
!> namespace. The Laplace path solves one field over all the bodies, so there
!> touching windings also couple across a shared face that is not a Dirichlet
!> face of the field.
!------------------------------------------------------------------------------
  SUBROUTINE DirectionByDistance()
!------------------------------------------------------------------------------
    IMPLICIT NONE
!------------------------------------------------------------------------------
    INTEGER, PARAMETER :: TriCorners(3,2) = RESHAPE([1,2,3,1,3,4],[3,2])
    TYPE(Element_t), POINTER :: Element
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Variable_t), POINTER :: Var, DirVar
    INTEGER, POINTER :: Perm(:), Indexes(:)
    INTEGER :: i, j, k, l, t, n, b, pass, nbody, nface, ntri, nlo, nhi, nparents, &
        nmore, nshared, ndiff, kcur
    INTEGER :: Parents(2)
    INTEGER, ALLOCATABLE :: TriBodyLo(:), TriBodyHi(:), OffLo(:), OffHi(:), NodeBody(:), &
        NodeMore(:), MoreBody(:), MoreNext(:)
    LOGICAL :: IsLo
    REAL(KIND=dp) :: BVals(Mesh % MaxElementNodes)
    REAL(KIND=dp), ALLOCATABLE :: TriLo(:,:), TriHi(:,:), CenLo(:,:), CenHi(:,:), &
        RadLo(:), RadHi(:), LowBound(:), BodyVlo(:), BodyVhi(:), VMin(:), MoreValue(:), &
        FirstValue(:), AllMin(:), AllMax(:), FaceMin(:), FaceMax(:), BIdMin(:), BIdMax(:)
    REAL(KIND=dp) :: v, dlo, dhi, dsum, dsmin, dsmax, bbox, span, SharedTol, p(3), t0, &
        maxdiff
!------------------------------------------------------------------------------
    t0 = RealTime()

    Var => Solver % Variable
    IF (Var % DOFs /= 1) CALL Fatal('DirectionSolver', &
        'Direction Method = distance needs a scalar variable')
    Perm => Var % Perm
    nbody = Model % NumberOfBodies

    ! The bodies each active node belongs to. A node of touching windings has
    ! several: the first is kept in NodeBody, the others in a list from NodeMore.
    !------------------------------------------------------------------------------
    ALLOCATE(ActiveBody(nbody), NodeBody(Mesh % NumberOfNodes), &
        NodeMore(Mesh % NumberOfNodes), MoreBody(64), MoreNext(64), MoreValue(64))
    ActiveBody = .FALSE.
    NodeBody = 0
    NodeMore = 0
    nmore = 0
    DO t=1,GetNOFActive()
      Element => GetActiveElement(t)
      b = Element % BodyId
      ActiveBody(b) = .TRUE.
      n = GetElementNOFNodes()
      DO i=1,n
        j = Element % NodeIndexes(i)
        IF (NodeBody(j) == 0) THEN
          NodeBody(j) = b
        ELSE IF (NodeBody(j) /= b) THEN
          k = NodeMore(j)
          DO WHILE (k > 0)
            IF (MoreBody(k) == b) EXIT
            k = MoreNext(k)
          END DO
          IF (k == 0) THEN
            IF (nmore == SIZE(MoreBody)) CALL GrowNodeBodyList(MoreBody, MoreNext, MoreValue)
            nmore = nmore + 1
            MoreBody(nmore) = b
            MoreNext(nmore) = NodeMore(j)
            NodeMore(j) = nmore
          END IF
        END IF
      END DO
    END DO

    ! The pair of boundary values of each body, from the Dirichlet faces it is
    ! a parent of. A face between two bodies is a face of both.
    !------------------------------------------------------------------------------
    ALLOCATE(BodyVlo(nbody), BodyVhi(nbody))
    BodyVlo = HUGE(v)
    BodyVhi = -HUGE(v)
    nface = 0
    DO t=1,Mesh % NumberOfBoundaryElements
      Element => GetBoundaryElement(t)
      IF (.NOT. DirichletFace(Element, BC)) CYCLE
      nparents = ActiveParents(Element, Parents)
      IF (nparents == 0) CYCLE

      n = GetElementNOFNodes()
      BVals(1:n) = GetReal(BC, varname, Found)
      IF (.NOT. Found) CYCLE
      IF (MAXVAL(BVals(1:n)) - MINVAL(BVals(1:n)) > BC_VALUE_TOL) THEN
        WRITE(Message,'(A,ES15.8,A,ES15.8)') 'Direction Method = distance needs constant &
            &boundary values, but '//TRIM(varname)//' varies from ', &
            MINVAL(BVals(1:n)),' to ',MAXVAL(BVals(1:n))
        CALL Fatal('DirectionSolver', Message)
      END IF

      nface = nface + 1
      DO k=1,nparents
        b = Parents(k)
        BodyVlo(b) = MIN(BodyVlo(b), BVals(1))
        BodyVhi(b) = MAX(BodyVhi(b), BVals(1))
      END DO
    END DO

    IF (ParallelReduction(nface) == 0) CALL Fatal('DirectionSolver', &
        'Direction Method = distance found no boundary conditions for '//TRIM(varname))

    IF (ParEnv % PEs > 1) THEN
      DO b=1,nbody
        BodyVlo(b) = ParallelReduction(BodyVlo(b),1)
        BodyVhi(b) = ParallelReduction(BodyVhi(b),2)
      END DO
    END IF

    ! A body without faces is reported with the face counts further down.
    span = 0.0_dp
    DO b=1,nbody
      IF (BodyVhi(b) < BodyVlo(b)) CYCLE
      span = MAX(span, BodyVhi(b) - BodyVlo(b))
      IF (.NOT. ActiveBody(b) .OR. BodyVhi(b) - BodyVlo(b) > BC_VALUE_TOL) CYCLE
      WRITE(Message,'(A,I0,A,ES15.8)') 'Direction Method = distance needs two distinct &
          &boundary values on every body, but body ',b,' has only ', BodyVlo(b)
      CALL Fatal('DirectionSolver', Message)
    END DO
    SharedTol = REL_SHARED_VALUE_TOL * span

    ! Collect the faces of both ends of every body as triangles tagged by body,
    ! counting them first. A face is classified per body: in a chain the shared
    ! face is the high end of one winding and the low end of the next.
    !------------------------------------------------------------------------------
    DO pass=1,2
      nlo = 0
      nhi = 0
      DO t=1,Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(t)
        IF (.NOT. DirichletFace(Element, BC)) CYCLE
        nparents = ActiveParents(Element, Parents)
        IF (nparents == 0) CYCLE

        n = GetElementNOFNodes()
        BVals(1:n) = GetReal(BC, varname, Found)
        IF (.NOT. Found) CYCLE
        v = BVals(1)

        ntri = 1
        IF (GetElementFamily() == 4) ntri = 2
        Indexes => Element % NodeIndexes

        DO k=1,nparents
          b = Parents(k)
          IF (ABS(v-BodyVlo(b)) <= BC_VALUE_TOL) THEN
            IsLo = .TRUE.
          ELSE IF (ABS(v-BodyVhi(b)) <= BC_VALUE_TOL) THEN
            IsLo = .FALSE.
          ELSE
            WRITE(Message,'(A,I0,A,ES15.8,A,ES15.8,A,ES15.8)') 'Direction Method = distance &
                &needs exactly two boundary values on every body, but body ',b,' has ', &
                BodyVlo(b),' and ',BodyVhi(b),' and ',v
            CALL Fatal('DirectionSolver', Message)
          END IF

          DO i=1,ntri
            IF (IsLo) THEN
              nlo = nlo + 1
              IF (pass == 2) THEN
                CALL SetTriangle(TriLo(:,nlo), Indexes(TriCorners(:,i)))
                TriBodyLo(nlo) = b
              END IF
            ELSE
              nhi = nhi + 1
              IF (pass == 2) THEN
                CALL SetTriangle(TriHi(:,nhi), Indexes(TriCorners(:,i)))
                TriBodyHi(nhi) = b
              END IF
            END IF
          END DO
        END DO
      END DO

      IF (pass == 1) THEN
        ALLOCATE(TriLo(9,MAX(nlo,1)), TriHi(9,MAX(nhi,1)), &
            TriBodyLo(MAX(nlo,1)), TriBodyHi(MAX(nhi,1)))
      END IF
    END DO

    ! Every partition needs the complete surfaces; duplicated faces are harmless
    ! for a minimum distance. A face between bodies of two partitions is saved
    ! in both, each with its own parent, so it reaches both bodies here.
    !------------------------------------------------------------------------------
    IF (ParEnv % PEs > 1) THEN
      CALL GatherTriangles(TriLo, TriBodyLo, nlo)
      CALL GatherTriangles(TriHi, TriBodyHi, nhi)
    END IF

    ALLOCATE(OffLo(nbody+1), OffHi(nbody+1))
    CALL GroupByBody(TriLo, TriBodyLo, nlo, OffLo)
    CALL GroupByBody(TriHi, TriBodyHi, nhi, OffHi)

    DO b=1,nbody
      IF (.NOT. ActiveBody(b)) CYCLE
      IF (OffLo(b+1) > OffLo(b) .AND. OffHi(b+1) > OffHi(b)) CYCLE
      WRITE(Message,'(A,I0,A,I0,A,I0,A)') 'Direction Method = distance needs both boundary &
          &values on every body, but body ',b,' has ',OffLo(b+1)-OffLo(b),' and ', &
          OffHi(b+1)-OffHi(b),' faces'
      CALL Fatal('DirectionSolver', Message)
    END DO

    ALLOCATE(CenLo(3,MAX(nlo,1)), RadLo(MAX(nlo,1)), CenHi(3,MAX(nhi,1)), RadHi(MAX(nhi,1)), &
        LowBound(MAX(nlo,nhi,1)))
    CALL BoundingSpheres(TriLo, nlo, CenLo, RadLo)
    CALL BoundingSpheres(TriHi, nhi, CenHi, RadHi)

    n = Mesh % NumberOfNodes
    bbox = (MAXVAL(Mesh % Nodes % x(1:n)) - MINVAL(Mesh % Nodes % x(1:n)))**2 + &
           (MAXVAL(Mesh % Nodes % y(1:n)) - MINVAL(Mesh % Nodes % y(1:n)))**2 + &
           (MAXVAL(Mesh % Nodes % z(1:n)) - MINVAL(Mesh % Nodes % z(1:n)))**2
    bbox = SQRT(bbox)
    IF (ParEnv % PEs > 1) bbox = ParallelReduction(bbox,2)

    ! Each node is measured against the faces of its first body, and a shared
    ! node against those of every further body too. The value of each body is
    ! kept, and the ends of the values at the node: of all its bodies, and of the
    ! bodies that have the node on a face with a boundary value.
    !------------------------------------------------------------------------------
    ALLOCATE(FirstValue(SIZE(Var % Values)), AllMin(SIZE(Var % Values)), &
        AllMax(SIZE(Var % Values)), FaceMin(SIZE(Var % Values)), &
        FaceMax(SIZE(Var % Values)), BIdMin(SIZE(Var % Values)), BIdMax(SIZE(Var % Values)))
    FirstValue = 0.0_dp
    AllMin = HUGE(v)
    AllMax = -HUGE(v)
    FaceMin = HUGE(v)
    FaceMax = -HUGE(v)
    BIdMin = HUGE(v)
    BIdMax = -HUGE(v)
    dsmin = HUGE(dsmin)
    dsmax = 0.0_dp
    DO i=1,n
      j = Perm(i)
      IF (j == 0) CYCLE
      b = NodeBody(i)
      IF (b == 0) CALL Fatal('DirectionSolver', &
          'Direction Method = distance found a dof outside the active bodies')
      p = [Mesh % Nodes % x(i), Mesh % Nodes % y(i), Mesh % Nodes % z(i)]
      kcur = 0
      k = NodeMore(i)
      DO
        dlo = MinTriangleDistance(p, TriLo, CenLo, RadLo, OffLo(b)+1, OffLo(b+1), LowBound)
        dhi = MinTriangleDistance(p, TriHi, CenHi, RadHi, OffHi(b)+1, OffHi(b+1), LowBound)
        dsum = dlo + dhi
        IF (dsum < REL_THICKNESS_TOL * bbox) THEN
          v = 0.5_dp * (BodyVlo(b) + BodyVhi(b))
        ELSE
          v = BodyVlo(b) + (BodyVhi(b) - BodyVlo(b)) * dlo / dsum
        END IF
        dsmin = MIN(dsmin, dsum)
        dsmax = MAX(dsmax, dsum)

        IF (kcur == 0) THEN
          Var % Values(j) = v
          FirstValue(j) = v
        ELSE
          MoreValue(kcur) = v
        END IF
        AllMin(j) = MIN(AllMin(j), v)
        AllMax(j) = MAX(AllMax(j), v)
        BIdMin(j) = MIN(BIdMin(j), REAL(b,dp))
        BIdMax(j) = MAX(BIdMax(j), REAL(b,dp))
        ! A corner node of a face is a vertex of its triangles, at distance zero.
        IF (dlo == 0.0_dp .OR. dhi == 0.0_dp) THEN
          FaceMin(j) = MIN(FaceMin(j), v)
          FaceMax(j) = MAX(FaceMax(j), v)
        END IF

        IF (k == 0) EXIT
        b = MoreBody(k)
        kcur = k
        k = MoreNext(k)
      END DO
    END DO

    ! A node on a partition boundary may have its bodies in different partitions,
    ! so the ends are taken over the partitions too. The nodal field keeps the
    ! smallest first-body value, the same in every partition.
    !------------------------------------------------------------------------------
    IF (ParEnv % PEs > 1) THEN
      ALLOCATE(VMin(SIZE(Var % Values)))
      VMin = Var % Values
      CALL ParallelSumNodalVector(Mesh, VMin, Perm, Op=1)
      Var % Values = VMin
      DEALLOCATE(VMin)
      CALL ParallelSumNodalVector(Mesh, AllMin, Perm, Op=1)
      CALL ParallelSumNodalVector(Mesh, AllMax, Perm, Op=2)
      CALL ParallelSumNodalVector(Mesh, FaceMin, Perm, Op=1)
      CALL ParallelSumNodalVector(Mesh, FaceMax, Perm, Op=2)
      CALL ParallelSumNodalVector(Mesh, BIdMin, Perm, Op=1)
      CALL ParallelSumNodalVector(Mesh, BIdMax, Perm, Op=2)
    END IF

    nshared = 0
    ndiff = 0
    maxdiff = 0.0_dp
    DO i=1,n
      j = Perm(i)
      IF (j == 0) CYCLE
      IF (FaceMax(j) > FaceMin(j) + SharedTol) THEN
        WRITE(Message,'(A,I0,A,ES13.6,A,ES13.6,A)') 'Direction Method = distance: node ', &
            GlobalNode(i),' is on faces with a boundary value of two touching bodies, '// &
            TRIM(varname)//' = ',FaceMin(j),' and ',FaceMax(j),'. Touching windings need '// &
            'one value on the face they share.'
        CALL Fatal('DirectionSolver', Message)
      END IF
      IF (ParEnv % PEs > 1) THEN
        IF (Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) /= ParEnv % MyPE) CYCLE
      END IF
      IF (BIdMax(j) > BIdMin(j)) nshared = nshared + 1
      IF (AllMax(j) - AllMin(j) > BC_VALUE_TOL) THEN
        ndiff = ndiff + 1
        maxdiff = MAX(maxdiff, AllMax(j) - AllMin(j))
      END IF
    END DO

    ! The elemental copy gives each element the values of its own body.
    !------------------------------------------------------------------------------
    DirVar => VariableGet(Mesh % Variables, TRIM(varname)//' Direction')
    IF (.NOT. ASSOCIATED(DirVar) .AND. GetNOFActive() > 0) CALL Fatal('DirectionSolver', &
        'Direction variable not found')
    DO t=1,GetNOFActive()
      Element => GetActiveElement(t)
      b = Element % BodyId
      DO l=1,GetElementNOFNodes()
        i = Element % NodeIndexes(l)
        IF (NodeBody(i) == b) THEN
          v = FirstValue(Perm(i))
        ELSE
          k = NodeMore(i)
          DO WHILE (MoreBody(k) /= b)
            k = MoreNext(k)
          END DO
          v = MoreValue(k)
        END IF
        j = DirVar % Perm(Element % DGIndexes(l))
        DirVar % Values(DirVar % DOFs*(j-1)+1:DirVar % DOFs*j) = v
      END DO
    END DO

    Var % Norm = ComputeNorm(Solver, SIZE(Var % Values), Var % Values)

    IF (ParEnv % PEs > 1) THEN
      dsmin = ParallelReduction(dsmin,1)
      dsmax = ParallelReduction(dsmax,2)
    END IF
    nshared = ParallelReduction(nshared)
    ndiff = ParallelReduction(ndiff)
    IF (ParEnv % PEs > 1) maxdiff = ParallelReduction(maxdiff,2)

    WRITE(Message,'(A,I0,A,I0,A,I0,A)') 'Distance direction for '//TRIM(varname)//' from ', &
        nlo,' and ',nhi,' boundary triangles on ',COUNT(ActiveBody),' bodies'
    CALL Info('DirectionSolver', Message, Level=5)
    IF (nshared > 0) THEN
      WRITE(Message,'(A,I0,A,I0,A,ES12.5)') 'Nodes shared by touching bodies: ',nshared, &
          ', with a different value in each body: ',ndiff,', largest difference ',maxdiff
      CALL Info('DirectionSolver', Message, Level=5)
    END IF
    WRITE(Message,'(A,ES12.5,A,ES12.5)') 'Wall thickness between the faces ranges from ', &
        dsmin,' to ',dsmax
    CALL Info('DirectionSolver', Message, Level=5)
    WRITE(Message,'(A,F8.3,A)') 'Distance direction computed in ',RealTime()-t0,' s'
    CALL Info('DirectionSolver', Message, Level=5)

    DEALLOCATE(TriLo, TriHi, TriBodyLo, TriBodyHi, OffLo, OffHi, CenLo, RadLo, CenHi, RadHi, &
        LowBound, ActiveBody, NodeBody, NodeMore, MoreBody, MoreNext, MoreValue, BodyVlo, &
        BodyVhi, FirstValue, AllMin, AllMax, FaceMin, FaceMax, BIdMin, BIdMax)
!------------------------------------------------------------------------------
  END SUBROUTINE DirectionByDistance
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Double the room of the list of further bodies of the shared nodes.
!------------------------------------------------------------------------------
  SUBROUTINE GrowNodeBodyList(MoreBody, MoreNext, MoreValue)
!------------------------------------------------------------------------------
    INTEGER, ALLOCATABLE :: MoreBody(:), MoreNext(:)
    REAL(KIND=dp), ALLOCATABLE :: MoreValue(:)
!------------------------------------------------------------------------------
    INTEGER, ALLOCATABLE :: Tmp(:)
    REAL(KIND=dp), ALLOCATABLE :: TmpR(:)
    INTEGER :: m
!------------------------------------------------------------------------------
    m = SIZE(MoreBody)
    ALLOCATE(Tmp(2*m))
    Tmp(1:m) = MoreBody
    CALL MOVE_ALLOC(Tmp, MoreBody)
    ALLOCATE(Tmp(2*m))
    Tmp(1:m) = MoreNext
    CALL MOVE_ALLOC(Tmp, MoreNext)
    ALLOCATE(TmpR(2*m))
    TmpR(1:m) = MoreValue
    CALL MOVE_ALLOC(TmpR, MoreValue)
!------------------------------------------------------------------------------
  END SUBROUTINE GrowNodeBodyList
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Node number for messages, the global one in parallel.
!------------------------------------------------------------------------------
  FUNCTION GlobalNode(i) RESULT(k)
!------------------------------------------------------------------------------
    INTEGER :: i, k
!------------------------------------------------------------------------------
    k = i
    IF (ParEnv % PEs > 1) k = Mesh % ParallelInfo % GlobalDOFs(i)
!------------------------------------------------------------------------------
  END FUNCTION GlobalNode
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> True for a face that carries a Dirichlet condition of the solver variable.
!------------------------------------------------------------------------------
  FUNCTION DirichletFace(Element, BC) RESULT(Found)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(ValueList_t), POINTER :: BC
    LOGICAL :: Found
!------------------------------------------------------------------------------
    INTEGER :: family
!------------------------------------------------------------------------------
    Found = .FALSE.
    IF (GetElementFamily() == 1) RETURN
    IF (.NOT. ActiveBoundaryElement()) RETURN
    BC => GetBC()
    IF (.NOT. ASSOCIATED(BC)) RETURN
    IF (.NOT. ListCheckPresent(BC, varname)) RETURN

    family = GetElementFamily()
    IF (family == 2) CALL Fatal('DirectionSolver', &
        'Direction Method = distance is 3D only')
    IF (family /= 3 .AND. family /= 4) CALL Fatal('DirectionSolver', &
        'Direction Method = distance cannot triangulate boundary element type '&
        //I2S(Element % TYPE % ElementCode))
    Found = .TRUE.
!------------------------------------------------------------------------------
  END FUNCTION DirichletFace
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> The bodies of the parent elements of a face that are active for this solver.
!> A face between two active bodies belongs to both of them.
!------------------------------------------------------------------------------
  FUNCTION ActiveParents(Element, Parents) RESULT(nparents)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER :: Parents(2), nparents
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Parent
    INTEGER :: i, b
!------------------------------------------------------------------------------
    nparents = 0
    IF (.NOT. ASSOCIATED(Element % BoundaryInfo)) RETURN

    DO i=1,2
      IF (i == 1) THEN
        Parent => Element % BoundaryInfo % Left
      ELSE
        Parent => Element % BoundaryInfo % Right
      END IF
      IF (.NOT. ASSOCIATED(Parent)) CYCLE

      b = Parent % BodyId
      IF (b < 1 .OR. b > SIZE(ActiveBody)) CYCLE
      IF (.NOT. ActiveBody(b)) CYCLE
      IF (nparents == 1) THEN
        IF (Parents(1) == b) CYCLE
      END IF
      nparents = nparents + 1
      Parents(nparents) = b
    END DO
!------------------------------------------------------------------------------
  END FUNCTION ActiveParents
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE SetTriangle(Tri, Indexes)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Tri(9)
    INTEGER :: Indexes(3)
!------------------------------------------------------------------------------
    INTEGER :: i, j
!------------------------------------------------------------------------------
    DO i=1,3
      j = Indexes(i)
      Tri(3*i-2) = Mesh % Nodes % x(j)
      Tri(3*i-1) = Mesh % Nodes % y(j)
      Tri(3*i)   = Mesh % Nodes % z(j)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE SetTriangle
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE BoundingSpheres(Tri, ntri, Cen, Rad)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Tri(:,:), Cen(:,:), Rad(:)
    INTEGER :: ntri
!------------------------------------------------------------------------------
    INTEGER :: i, k
!------------------------------------------------------------------------------
    DO k=1,ntri
      Cen(:,k) = (Tri(1:3,k) + Tri(4:6,k) + Tri(7:9,k)) / 3.0_dp
      Rad(k) = 0.0_dp
      DO i=1,3
        Rad(k) = MAX(Rad(k), SQRT(SUM((Tri(3*i-2:3*i,k) - Cen(:,k))**2)))
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE BoundingSpheres
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Sort the triangles by the body they belong to and return the body offsets, so
!> that body b owns the triangles Offset(b)+1 ... Offset(b+1).
!------------------------------------------------------------------------------
  SUBROUTINE GroupByBody(Tri, BodyId, ntri, Offset)
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE :: Tri(:,:)
    INTEGER, ALLOCATABLE :: BodyId(:)
    INTEGER :: ntri, Offset(:)
!------------------------------------------------------------------------------
    INTEGER :: b, k
    INTEGER, ALLOCATABLE :: Pos(:)
    REAL(KIND=dp), ALLOCATABLE :: Sorted(:,:)
!------------------------------------------------------------------------------
    Offset = 0
    DO k=1,ntri
      Offset(BodyId(k)+1) = Offset(BodyId(k)+1) + 1
    END DO
    DO b=1,SIZE(Offset)-1
      Offset(b+1) = Offset(b+1) + Offset(b)
    END DO

    ALLOCATE(Sorted(9,MAX(ntri,1)), Pos(SIZE(Offset)))
    Pos = Offset
    DO k=1,ntri
      b = BodyId(k)
      Pos(b) = Pos(b) + 1
      Sorted(:,Pos(b)) = Tri(:,k)
    END DO

    DEALLOCATE(Tri, Pos)
    CALL MOVE_ALLOC(Sorted, Tri)
!------------------------------------------------------------------------------
  END SUBROUTINE GroupByBody
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Gather the triangles of all partitions, with the body of each triangle, so that
!> every partition holds the complete surfaces.
!------------------------------------------------------------------------------
  SUBROUTINE GatherTriangles(Tri, BodyId, ntri)
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE :: Tri(:,:)
    INTEGER, ALLOCATABLE :: BodyId(:)
    INTEGER :: ntri
!------------------------------------------------------------------------------
    INTEGER :: i, ierr, comm, pes
    INTEGER, ALLOCATABLE :: Counts(:), Displs(:), Counts9(:), Displs9(:), AllBodyId(:)
    REAL(KIND=dp), ALLOCATABLE :: AllTri(:,:)
!------------------------------------------------------------------------------
    comm = ParEnv % ActiveComm
    CALL MPI_COMM_SIZE(comm, pes, ierr)

    ALLOCATE(Counts(pes), Displs(pes), Counts9(pes), Displs9(pes))
    CALL MPI_ALLGATHER(ntri, 1, MPI_INTEGER, Counts, 1, MPI_INTEGER, comm, ierr)

    Displs(1) = 0
    DO i=2,pes
      Displs(i) = Displs(i-1) + Counts(i-1)
    END DO
    Counts9 = 9 * Counts
    Displs9 = 9 * Displs

    ALLOCATE(AllTri(9,MAX(SUM(Counts),1)), AllBodyId(MAX(SUM(Counts),1)))
    CALL MPI_ALLGATHERV(Tri, 9*ntri, MPI_DOUBLE_PRECISION, AllTri, Counts9, Displs9, &
        MPI_DOUBLE_PRECISION, comm, ierr)
    CALL MPI_ALLGATHERV(BodyId, ntri, MPI_INTEGER, AllBodyId, Counts, Displs, &
        MPI_INTEGER, comm, ierr)

    ntri = SUM(Counts)
    DEALLOCATE(Tri, BodyId, Counts, Displs, Counts9, Displs9)
    CALL MOVE_ALLOC(AllTri, Tri)
    CALL MOVE_ALLOC(AllBodyId, BodyId)
!------------------------------------------------------------------------------
  END SUBROUTINE GatherTriangles
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Smallest distance from a point to the triangles i1 ... i2 of a set. The bounding
!> spheres give a lower bound for each triangle, and only the triangles whose lower
!> bound is below the distance of the most promising one are evaluated exactly.
!------------------------------------------------------------------------------
  FUNCTION MinTriangleDistance(p, Tri, Cen, Rad, i1, i2, LowBound) RESULT(dmin)
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: p(3), Tri(:,:), Cen(:,:), Rad(:), LowBound(:), dmin
    INTEGER :: i1, i2
!------------------------------------------------------------------------------
    INTEGER :: k, kbest
!------------------------------------------------------------------------------
    kbest = i1
    DO k=i1,i2
      LowBound(k) = SQRT(SUM((p - Cen(:,k))**2)) - Rad(k)
      IF (LowBound(k) < LowBound(kbest)) kbest = k
    END DO

    dmin = PointTriangleDistance(p, Tri(1:3,kbest), Tri(4:6,kbest), Tri(7:9,kbest))
    DO k=i1,i2
      IF (LowBound(k) >= dmin) CYCLE
      dmin = MIN(dmin, PointTriangleDistance(p, Tri(1:3,k), Tri(4:6,k), Tri(7:9,k)))
    END DO
!------------------------------------------------------------------------------
  END FUNCTION MinTriangleDistance
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Exact distance from a point to a triangle. The closest point is picked by the
!> barycentric region tests of Ericson, Real-Time Collision Detection, 5.1.5.
!------------------------------------------------------------------------------
  PURE FUNCTION PointTriangleDistance(p, a, b, c) RESULT(dist)
!------------------------------------------------------------------------------
    REAL(KIND=dp), INTENT(IN) :: p(3), a(3), b(3), c(3)
    REAL(KIND=dp) :: dist
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: ab(3), ac(3), ap(3), bp(3), cp(3), q(3)
    REAL(KIND=dp) :: d1, d2, d3, d4, d5, d6, va, vb, vc, denom
!------------------------------------------------------------------------------
    ab = b - a
    ac = c - a
    ap = p - a
    d1 = SUM(ab*ap)
    d2 = SUM(ac*ap)

    IF (d1 <= 0.0_dp .AND. d2 <= 0.0_dp) THEN
      q = a
    ELSE
      bp = p - b
      d3 = SUM(ab*bp)
      d4 = SUM(ac*bp)
      IF (d3 >= 0.0_dp .AND. d4 <= d3) THEN
        q = b
      ELSE
        vc = d1*d4 - d3*d2
        IF (vc <= 0.0_dp .AND. d1 >= 0.0_dp .AND. d3 <= 0.0_dp) THEN
          q = a + (d1/(d1-d3)) * ab
        ELSE
          cp = p - c
          d5 = SUM(ab*cp)
          d6 = SUM(ac*cp)
          IF (d6 >= 0.0_dp .AND. d5 <= d6) THEN
            q = c
          ELSE
            vb = d5*d2 - d1*d6
            IF (vb <= 0.0_dp .AND. d2 >= 0.0_dp .AND. d6 <= 0.0_dp) THEN
              q = a + (d2/(d2-d6)) * ac
            ELSE
              va = d3*d6 - d5*d4
              IF (va <= 0.0_dp .AND. d4-d3 >= 0.0_dp .AND. d5-d6 >= 0.0_dp) THEN
                q = b + ((d4-d3)/((d4-d3)+(d5-d6))) * (c-b)
              ELSE
                denom = va + vb + vc
                IF (denom > 0.0_dp) THEN
                  q = a + ab*(vb/denom) + ac*(vc/denom)
                ELSE
                  q = a  ! zero area face, the vertex regions already cover it
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
    END IF

    dist = SQRT(SUM((p-q)**2))
!------------------------------------------------------------------------------
  END FUNCTION PointTriangleDistance
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE SaveSolutionWithBodyMethod(ns_iter)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  INTEGER :: Active, n, t, nn, ns_iter
  TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
  Active = GetNOFActive()
  DO t=1,Active
     Element => GetActiveElement(t)
     IF (ns_iter /= Element % BodyId) CYCLE
     nn = GetElementNOFNodes()
     CALL SaveElementSolution(Element, nn)
  END DO
!------------------------------------------------------------------------------
  END SUBROUTINE SaveSolutionWithBodyMethod
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE SaveElementSolution(Element, nn)
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  INTEGER :: nn, j, k
  TYPE(Element_t), POINTER :: Element
  TYPE(Valuelist_t), POINTER :: Solverparams
  TYPE(Variable_t), POINTER, SAVE :: directionvar
  LOGICAL, SAVE :: First=.TRUE.
  CHARACTER(LEN=MAX_NAME_LEN), SAVE :: varname
  REAL(KIND=dp) :: direction(nn)

  IF (varname .ne. GetString(GetSolverParams(), 'Variable', Found)) & 
          First =.TRUE.
  varname = GetString(GetSolverParams(), 'Variable', Found)
  IF (First) THEN
    First = .FALSE.
    directionvar => VariableGet( Mesh % Variables, TRIM(varname)//' Direction')
    IF(.NOT. ASSOCIATED(directionvar)) THEN
      CALL Fatal('SaveElementSolution()','Direction variable not found')
    END IF
  END IF
 
  CALL GetLocalSolution(direction, varname)
  DO j=1,nn
    IF (ASSOCIATED(directionvar)) THEN
      DO k=1,directionvar % DOFs
        directionvar % Values( directionvar % DOFs*(directionvar % Perm( &
              Element % DGIndexes(j))-1)+k) = direction(j)
      END DO
    END IF
  END DO 
!------------------------------------------------------------------------------
  END SUBROUTINE SaveElementSolution
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, LOAD, Element, n, nd )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    LOGICAL :: Stat
    INTEGER :: i,t
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
       IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      ! Finally, the elemental matrix & vector:
      !----------------------------------------
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
             MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )
      FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * LoadAtIP * Basis(1:nd)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!----------------------------------------------------------------
  SUBROUTINE BoundaryCondition(LOAD, FORCE, Element, n)
!----------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), DIMENSION(:) :: FORCE, LOAD
    INTEGER :: n
    TYPE(Element_t), POINTER :: Element
!----------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n),dBasisdx(n,3)
    REAL(KIND=dp) :: detJ, LoadAtIP,&
    LocalHeatCapacity, LocalDensity
    LOGICAL :: stat, getSecondDerivatives
    INTEGER :: t,j
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!----------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    FORCE = 0.0d0
    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    !-----------------------------------------------------------------
    ! Loop over Gauss-points (boundary element Integration)
    !-----------------------------------------------------------------
    DO t=1,IP % n
      !Basis function values & derivatives at the integration point:
      !-------------------------------------------------------------
      getSecondDerivatives = .FALSE.
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
      IP % W(t), detJ, Basis, dBasisdx)

      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      DO j=1,n
        FORCE(j) = FORCE(j) + IP % s(t)*DetJ*LoadAtIP*Basis(j)
      END DO
    END DO
  END SUBROUTINE BoundaryCondition

!------------------------------------------------------------------------------
END SUBROUTINE DirectionSolver
!------------------------------------------------------------------------------
