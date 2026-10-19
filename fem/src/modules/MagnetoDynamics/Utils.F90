!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!>  Utilities for the A-V solvers of electromagnetism
!> \ingroup Solvers
!-------------------------------------------------------------------------------
MODULE MagnetoDynamicsUtils

   USE MGDynMaterialUtils
   IMPLICIT NONE

   INTEGER :: JfixPhase
   REAL(KIND=dp), POINTER :: Jfixrhs(:)
   COMPLEX(KIND=dp), POINTER :: JfixRhsC(:)
   INTEGER, POINTER :: JfixSurfacePerm(:) 
   REAL(KIND=dp), ALLOCATABLE, TARGET :: JfixSurfaceVec(:)
   COMPLEX(KIND=dp), ALLOCATABLE :: JfixSurfaceVecC(:)
   
   COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)

   ! A transient 'layered' laminate stack carries the eddy currents of its
   ! laminations as a Foster ladder per element: LaminateXiName holds the states
   ! (three global components per pole), LaminateLadderName the ladder of the
   ! timestep. WhitneyAVSolver_Init declares both on the bodies whose material it
   ! marks with LaminateLadderMask.
   CHARACTER(*), PARAMETER :: LaminateXiName = 'Laminate Xi', &
       LaminateLadderName = 'Laminate Ladder', LaminateLadderMask = 'Laminate Ladder Body'
   INTEGER, PARAMETER :: DefaultLaminateLadderOrder = 6
   ! LaminateLadderName per element: the ladder's relative permeability, Kappa,
   ! the loss density of the last solve, the dynamic reluctivity of the step and
   ! its history field (3), then the residues and the time constants of the
   ! poles (as many of each as the variable has room for).
   INTEGER, PARAMETER :: LamLadMu = 1, LamLadKappa = 2, LamLadLoss = 3, &
       LamLadNuDyn = 4, LamLadHist = 5, LamLadPoles = 8
   ! LaminateXiName per element: the states of the poles, three global
   ! components each, then the element mean of B_par (3) and the stored energy of
   ! the ladder (1) of the last solve, committed with the states for the loss
   INTEGER, PARAMETER :: LamXiExtra = 4

   ! The in-plane energy integral of the static 'layered' law at the knots of
   ! each H-B curve and stacking factor (LaminateInPlaneEnergyIntegral), keyed on
   ! the curve's data: a reloaded model frees its curve entries and makes new ones
   TYPE LaminateEnergyTable_t
     REAL(KIND=dp) :: StackingFactor = -1.0_dp
     LOGICAL :: Cubic = .FALSE.
     REAL(KIND=dp), ALLOCATABLE :: Knots(:), Values(:), CubicCoeff(:), Cumulative(:)
   END TYPE LaminateEnergyTable_t
   INTEGER, PARAMETER :: MaxLaminateEnergyTables = 64
   TYPE(LaminateEnergyTable_t), SAVE :: LamEnergyTables(MaxLaminateEnergyTables)
   INTEGER, SAVE :: NoLamEnergyTables = 0, LastLamEnergyTable = 0
   ! WhitneyAVSolver_Init counts each declaration of the laminate ladder variables:
   ! a solver initialised again (a reloaded model) starts its ladders again
   INTEGER, SAVE :: LaminateInitGeneration = 0

   INTERFACE SetDOFtoValue
     MODULE PROCEDURE SetDOFtoValueR, SetDOFtoValueC
   END INTERFACE

   INTERFACE GetReluctivity
     MODULE PROCEDURE GetReluctivityR, GetReluctivityC, &
         GetReluctivityTensorR, GetReluctivityTensorC
   END INTERFACE

CONTAINS

  RECURSIVE FUNCTION AddConstraintFromBulk(A, M0) RESULT (M)
    USE DefUtils
    TYPE(Matrix_t), INTENT(IN) :: A
    TYPE(Matrix_t), POINTER, INTENT(INOUT) :: M0
    TYPE(Matrix_t), POINTER :: M

    integer :: n, n1

    IF(.NOT. A % FORMAT == MATRIX_CRS) CALL Fatal("AddConstraintFromBulk", "Matrix A is not a CRS matrix")

    IF (.NOT. ASSOCIATED(M0)) THEN
      CALL info("AddConstraintFromBulk", "M0 not associated, creating now.", level=10)
      M => AllocateMatrix()
      ! Rowdata copy
      n = A % NumberOfRows
      ALLOCATE(M % RHS(size(A % RHS)))
      ALLOCATE(M % Rows(size(A%rows)))
      ALLOCATE(M % ConstrainedDof(size(M % Rows))) 
      ALLOCATE(M % Diag(size(A % Diag)))

      M % RHS(:) = 0.0_dp
      M % ConstrainedDOF(:) = .FALSE.

      M % ListMatrix => NULL()
      M % NumberOfRows = A % NumberOfRows

      DO n1 = 1, SIZE(M % Rows)
        M % Rows(n1) = A % Rows(n1)
      END DO

      DO n1 = 1,SIZE(M % diag)
        M % Diag(n1) = A % Diag(n1) 
        IF (M % diag(n1) == 0) THEN
          WRITE (*,*) 'diag', n1, 'is zero'
        end if
      END DO

      ! Columndata copy
      n = size(A % Values)
      ALLOCATE(M % Values(n), M % DValues(n), M % Cols(n))
      M % Cols = A % Cols
      M % DValues = 0.0_dp
      M % Values = 0.0_dp

      !! Permutation copy
      IF(ASSOCIATED(A % Perm)) THEN
        CALL info("AddConstraintFromBulk", "Copying perm.", level=3)
        ALLOCATE(M % perm(size(A%perm)))
        M % perm = A % perm
     END IF

      !! Inverse Permutation copy
      IF(ASSOCIATED(A % InvPerm)) THEN
        CALL info("AddConstraintFromBulk", "Copying inverse perm.", level=3)
        ALLOCATE(M % InvPerm(size(A%Invperm)))
        M % InvPerm = A % InvPerm 
      END IF

      M % FORMAT = A % FORMAT
      M0 => M
    ELSE
      CALL info("AddConstraintFromBulk", "M0 is associated, recursing...", level=10)
      M => AddConstraintFromBulk(A, M0 % ConstraintMatrix)
    END IF
  END FUNCTION AddConstraintFromBulk


!------------------------------------------------------------------------------
  SUBROUTINE SetDOFToValueR(Solver,k,VALUE)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: VALUE,v
    TYPE(Solver_t) :: Solver
    INTEGER :: n,k
    TYPE(Matrix_t), POINTER :: A
    TYPE(Mesh_t), POINTER :: Mesh

    Mesh => GetMesh(Solver)
    n = Solver % Variable % Perm(k+Mesh % NumberOfNodes)
    A => GetMatrix()
    CALL CRS_SetSymmDirichlet(A,A % RHS,n,VALUE)
!------------------------------------------------------------------------------
  END SUBROUTINE SetDOFToValueR
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE SetDOFToValueC(Solver,k,VALUE)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(KIND=dp) :: VALUE
    TYPE(Solver_t) :: Solver
    INTEGER :: n,k
    TYPE(Matrix_t), POINTER :: A
    TYPE(Mesh_t), POINTER :: Mesh

    Mesh => GetMesh(Solver)
    n = Solver % Variable % Perm(k+Mesh % NumberOfNodes)
    A => GetMatrix()

    CALL CRS_SetSymmDirichlet(A,A % RHS,2*(n-1)+1,REAL(VALUE))
    CALL CRS_SetSymmDirichlet(A,A % RHS,2*(n-1)+2,AIMAG(VALUE))
!------------------------------------------------------------------------------
  END SUBROUTINE SetDOFToValueC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetReluctivityR(Material,Acoef,n)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp) :: Acoef(:)
    INTEGER :: n
!------------------------------------------------------------------------------
    LOGICAL :: Found, FirstTime = .TRUE.
    REAL(KIND=dp) :: Avacuum

    SAVE FirstTime, Avacuum 
!------------------------------------------------------------------------------

    IF ( FirstTime ) THEN
      Avacuum = GetConstReal( CurrentModel % Constants, &
              'Permeability of Vacuum', Found )
      IF (.NOT. Found ) Avacuum = PI * 4.0d-7
      FirstTime = .FALSE.
    END IF
  
    Acoef(1:n) = GetReal( Material, 'Relative Permeability', Found )
    IF ( Found ) THEN
      Acoef(1:n) = Avacuum * Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Permeability', Found )
    END IF
    IF ( Found ) THEN
      Acoef(1:n) = 1._dp / Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Reluctivity', Found )
    END IF
    IF( .NOT. Found .AND. &
        .NOT. ListCheckPresent(Material, 'H-B Curve') ) THEN
      CALL Fatal('GetReluctivityR','Give > Relative Permeability < or > Reluctivity <  for material!')
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivityR
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
 SUBROUTINE GetReluctivityC(Material,Acoef,n)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Material
    COMPLEX(KIND=dp) :: Acoef(:)
    INTEGER :: n
!------------------------------------------------------------------------------
    LOGICAL :: L, Found, FirstTime = .TRUE.
    REAL(KIND=dp) :: Avacuum

    SAVE Avacuum, FirstTime
!------------------------------------------------------------------------------

    IF ( FirstTime ) THEN
      Avacuum = GetConstReal( CurrentModel % Constants, &
              'Permeability of Vacuum', Found )
      IF(.NOT. Found ) Avacuum = PI * 4.0d-7
      FirstTime = .FALSE.
    END IF

    Acoef(1:n) = GetReal( Material, 'Relative Permeability', Found )
    IF ( Found ) THEN
      Acoef(1:n) = Avacuum * Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Permeability', Found )
    END IF
    IF ( Found ) THEN
      Acoef(1:n) = 1._dp / Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Reluctivity', Found )
      Acoef(1:n) = CMPLX( REAL(Acoef(1:n)), &
         GetReal( Material, 'Reluctivity im', L ), KIND=dp )
      Found = Found .OR. L
    END IF
    IF( .NOT. Found .AND. &
        .NOT. ListCheckPresent(Material, 'H-B Curve') ) THEN
      CALL Fatal('GetReluctivityC','Give > Relative Permeability < or > Reluctivity <  for material!')
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivityC
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Permeability of vacuum: the model constant if given, else 4 pi 1e-7.
!------------------------------------------------------------------------------
  FUNCTION VacuumPermeability() RESULT(mu0)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: mu0
    LOGICAL :: Found
!------------------------------------------------------------------------------
    mu0 = GetConstReal(CurrentModel % Constants, 'Permeability of Vacuum', Found)
    IF (.NOT. Found) mu0 = PI * 4.0d-7
!------------------------------------------------------------------------------
  END FUNCTION VacuumPermeability
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Conductance factor of the mid-surface sheet that carries the net current of
!> one lamination inside a 'layered' stack body: sinh(k d)/(k d). The stack
!> body already holds mu g over the lamination, which adds i omega mu g d/4 to
!> the net-current impedance between the two faces; with this factor the face
!> impedance is the exact one, x coth(x)/(sigma d). Low frequency limit 1.
!------------------------------------------------------------------------------
  FUNCTION LaminateSheetFactor(Omega, Mu, Sigma, Thickness) RESULT(f)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Omega, Sigma, Thickness
    COMPLEX(KIND=dp), INTENT(IN) :: Mu
    COMPLEX(KIND=dp) :: f
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: kd
    ! below SeriesLimit sinh(kd)/kd is taken from its series; above OverflowLimit it is capped
    REAL(KIND=dp), PARAMETER :: SeriesLimit = 1.0d-3, OverflowLimit = 600.0_dp
!------------------------------------------------------------------------------
    kd = SQRT(im * Omega * Mu * Sigma) * Thickness
    IF (ABS(kd) < SeriesLimit) THEN
      f = 1.0_dp + kd**2 / 6
    ELSE IF (REAL(kd) > OverflowLimit) THEN
      CALL Fatal('LaminateSheetFactor', 'Lamination is too thick for its skin depth (k d > 600)!')
    ELSE
      f = SINH(kd) / kd
    END IF
!------------------------------------------------------------------------------
  END FUNCTION LaminateSheetFactor
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Skin factor g = tanh(x)/x of one insulated lamination of thickness d, with
!> x = k d/2 and k = sqrt(i omega mu sigma) (time dependence exp(i omega t)).
!> mu*g is the mean permeability of the lamination when the tangential field is
!> the same on both faces; sigma*g is the conductivity that links its net sheet
!> current to the mean of the two face electric fields.
!------------------------------------------------------------------------------
  FUNCTION LaminateSkinFactor(Omega, Mu, Sigma, Thickness) RESULT(g)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Omega, Sigma, Thickness
    COMPLEX(KIND=dp), INTENT(IN) :: Mu
    COMPLEX(KIND=dp) :: g
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: x
    ! below SeriesLimit tanh(x)/x is taken from its series, above TanhOneLimit tanh(x) = 1
    REAL(KIND=dp), PARAMETER :: SeriesLimit = 1.0d-3, TanhOneLimit = 20.0_dp
!------------------------------------------------------------------------------
    x = SQRT(im * Omega * Mu * Sigma) * Thickness / 2
    IF (ABS(x) < SeriesLimit) THEN
      g = 1.0_dp - x**2 / 3
    ELSE IF (REAL(x) > TanhOneLimit) THEN
      g = 1.0_dp / x
    ELSE
      g = TANH(x) / x
    END IF
!------------------------------------------------------------------------------
  END FUNCTION LaminateSkinFactor
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Reads the parameters of the 'layered' laminate stack model of a material:
!> lamination thickness, lamination conductivity, stacking factor and the
!> optional lamination normal. Without a normal the stack is treated as bent
!> (a corner) and the in-plane reluctivity is used in every direction.
!------------------------------------------------------------------------------
  SUBROUTINE GetLaminateLayeredParameters(Material, Thickness, Sigma, &
      StackingFactor, Normal, HasNormal, LadderOrder)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp), INTENT(OUT) :: Thickness, Sigma, StackingFactor, Normal(3)
    LOGICAL, INTENT(OUT) :: HasNormal
    ! poles of the transient eddy current ladder
    INTEGER, INTENT(OUT), OPTIONAL :: LadderOrder
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: NormalArray(:,:)
    LOGICAL :: Found
    CHARACTER(*), PARAMETER :: Caller = 'GetLaminateLayeredParameters'
!------------------------------------------------------------------------------
    Thickness = GetConstReal(Material, 'Laminate Thickness', Found)
    IF (.NOT. Found .OR. Thickness <= 0.0_dp) &
        CALL Fatal(Caller, 'Give a positive > Laminate Thickness < for the layered model!')
    Sigma = GetConstReal(Material, 'Laminate Stack Conductivity', Found)
    IF (.NOT. Found .OR. Sigma < 0.0_dp) &
        CALL Fatal(Caller, 'Give > Laminate Stack Conductivity < for the layered model!')
    StackingFactor = GetConstReal(Material, 'Laminate Stacking Factor', Found)
    IF (.NOT. Found) StackingFactor = 1.0_dp
    IF (StackingFactor <= 0.0_dp .OR. StackingFactor > 1.0_dp) &
        CALL Fatal(Caller, '> Laminate Stacking Factor < must be in (0,1]!')

    Normal = 0.0_dp
    NormalArray => ListGetConstRealArray(Material, 'Laminate Normal', HasNormal)
    IF (HasNormal) THEN
      IF (SIZE(NormalArray,1) < 3) CALL Fatal(Caller, '> Laminate Normal < needs three components!')
      Normal = NormalArray(1:3,1)
      IF (SQRT(SUM(Normal**2)) < EPSILON(1.0_dp)) CALL Fatal(Caller, '> Laminate Normal < is a zero vector!')
      Normal = Normal / SQRT(SUM(Normal**2))
    END IF

    IF (PRESENT(LadderOrder)) THEN
      LadderOrder = ListGetInteger(Material, 'Laminate Ladder Order', Found)
      IF (.NOT. Found) LadderOrder = DefaultLaminateLadderOrder
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GetLaminateLayeredParameters
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Reluctivity tensor of a stack of insulated laminations ('layered' laminate
!> stack model): nu_par = 1/(F mu g + (1-F) mu0) in the lamination plane and
!> nu_perp = F nu + (1-F)/mu0 along the normal. NuMaterial is the (complex)
!> reluctivity of the lamination material. Without a normal nu_par is used in
!> every direction.
!------------------------------------------------------------------------------
  FUNCTION LaminateLayeredReluctivity(Omega, NuMaterial, Sigma, Thickness, &
      StackingFactor, Normal, HasNormal) RESULT(Nu)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Omega, Sigma, Thickness, StackingFactor, Normal(3)
    COMPLEX(KIND=dp), INTENT(IN) :: NuMaterial
    LOGICAL, INTENT(IN) :: HasNormal
    COMPLEX(KIND=dp) :: Nu(3,3)
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: mu, g, nu_par, nu_perp
    REAL(KIND=dp) :: mu0
    INTEGER :: i, j
!------------------------------------------------------------------------------
    mu0 = VacuumPermeability()

    mu = 1.0_dp / NuMaterial
    g = LaminateSkinFactor(Omega, mu, Sigma, Thickness)
    nu_par = 1.0_dp / (StackingFactor * mu * g + (1.0_dp - StackingFactor) * mu0)
    nu_perp = StackingFactor * NuMaterial + (1.0_dp - StackingFactor) / mu0

    Nu = 0.0_dp
    DO i = 1, 3
      Nu(i,i) = nu_par
    END DO
    IF (HasNormal) THEN
      DO i = 1, 3
        DO j = 1, 3
          Nu(i,j) = Nu(i,j) + (nu_perp - nu_par) * Normal(i) * Normal(j)
        END DO
      END DO
    END IF
!------------------------------------------------------------------------------
  END FUNCTION LaminateLayeredReluctivity
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Flux density in the laminations of a 'layered' stack whose homogenized field
!> is B: the laminations carry the in-plane flux of the stack, B_par/F, and all
!> of the normal flux B_perp. Without a normal every component is in-plane.
!------------------------------------------------------------------------------
  FUNCTION LaminateMetalFlux(B, StackingFactor, Normal, HasNormal) RESULT(Bs)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(KIND=dp), INTENT(IN) :: B(3)
    REAL(KIND=dp), INTENT(IN) :: StackingFactor, Normal(3)
    LOGICAL, INTENT(IN) :: HasNormal
    COMPLEX(KIND=dp) :: Bs(3)
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: Bperp(3)
!------------------------------------------------------------------------------
    Bperp = 0.0_dp
    IF (HasNormal) Bperp = Normal * SUM(Normal * B)
    Bs = (B - Bperp) / StackingFactor + Bperp
!------------------------------------------------------------------------------
  END FUNCTION LaminateMetalFlux
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Peak flux density |B_s| in the laminations of a 'layered' stack (see
!> LaminateMetalFlux). Weight = B_par/F^2 + B_perp, the map applied twice,
!> gives its change with B, dB_s = conj(Weight).dB/B_s, in the form the Newton
!> term of an H-B curve uses.
!------------------------------------------------------------------------------
  FUNCTION LaminateMetalFluxDensity(B, StackingFactor, Normal, HasNormal, Weight) RESULT(BsAbs)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(KIND=dp), INTENT(IN) :: B(3)
    REAL(KIND=dp), INTENT(IN) :: StackingFactor, Normal(3)
    LOGICAL, INTENT(IN) :: HasNormal
    COMPLEX(KIND=dp), INTENT(OUT), OPTIONAL :: Weight(3)
    REAL(KIND=dp) :: BsAbs
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: Bs(3)
!------------------------------------------------------------------------------
    Bs = LaminateMetalFlux(B, StackingFactor, Normal, HasNormal)
    BsAbs = SQRT(SUM(ABS(Bs)**2))
    IF (PRESENT(Weight)) Weight = LaminateMetalFlux(Bs, StackingFactor, Normal, HasNormal)
!------------------------------------------------------------------------------
  END FUNCTION LaminateMetalFluxDensity
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Derivative of the 'layered' stack tensor with respect to the reluctivity of
!> the lamination material, by a relative step: the material reluctivity also
!> enters through the skin factor g. For the Newton term of an H-B curve.
!------------------------------------------------------------------------------
  FUNCTION LaminateLayeredReluctivityDerivative(Omega, NuMaterial, Sigma, Thickness, &
      StackingFactor, Normal, HasNormal) RESULT(dNu)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Omega, Sigma, Thickness, StackingFactor, Normal(3)
    COMPLEX(KIND=dp), INTENT(IN) :: NuMaterial
    LOGICAL, INTENT(IN) :: HasNormal
    COMPLEX(KIND=dp) :: dNu(3,3)
!------------------------------------------------------------------------------
    REAL(KIND=dp), PARAMETER :: RelativeStep = 1.0e-6_dp
    COMPLEX(KIND=dp) :: Step
!------------------------------------------------------------------------------
    Step = RelativeStep * NuMaterial
    dNu = ( LaminateLayeredReluctivity(Omega, NuMaterial + Step, Sigma, Thickness, &
        StackingFactor, Normal, HasNormal) - LaminateLayeredReluctivity(Omega, &
        NuMaterial, Sigma, Thickness, StackingFactor, Normal, HasNormal) ) / Step
!------------------------------------------------------------------------------
  END FUNCTION LaminateLayeredReluctivityDerivative
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Reluctivity of the lamination steel at a lamination sheet of a 'layered'
!> stack: H(|B_s|)/|B_s| averaged over the volume of the face's parent elements
!> that are layered bodies with an H-B curve, from the potential of PSolver.
!> The sheet conductance then holds for the same permeability as the body.
!> Found is false when no parent has a curve.
!------------------------------------------------------------------------------
  FUNCTION LaminateSheetSteelReluctivity(Face, PSolver, PiolaVersion, EdgeBasisDegree, &
      Found) RESULT(NuSteel)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Face
    TYPE(Solver_t), POINTER :: PSolver
    LOGICAL, INTENT(IN) :: PiolaVersion
    INTEGER, INTENT(IN) :: EdgeBasisDegree
    LOGICAL, INTENT(OUT) :: Found
    REAL(KIND=dp) :: NuSteel
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Parent
    TYPE(ValueList_t), POINTER :: Material
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp), ALLOCATABLE :: Pot(:,:), Basis(:), dBasisdx(:,:), WBasis(:,:), RotWBasis(:,:)
    REAL(KIND=dp) :: Thickness, Sigma, Fill, Normal(3), detJ, BsAbs, NuSum, Volume
    COMPLEX(KIND=dp) :: B(3)
    INTEGER :: side, t, n, nd, np
    LOGICAL :: HasNormal, Stat, GotMaterial
    !$OMP THREADPRIVATE(Nodes)
!------------------------------------------------------------------------------
    Found = .FALSE.
    NuSteel = 0.0_dp
    NuSum = 0.0_dp
    Volume = 0.0_dp
    IF (.NOT. ASSOCIATED(Face % BoundaryInfo)) RETURN
    ! In parallel a face on a partition interface keeps only its local parent
    ! (ElmerGrid drops the other one without a halo): the average would differ
    ! from the serial one, or fall back to the sheet's own permeability.
    IF (ParEnv % PEs > 1) THEN
      IF (.NOT. (ASSOCIATED(Face % BoundaryInfo % Left) .AND. &
          ASSOCIATED(Face % BoundaryInfo % Right))) THEN
        IF (FaceOnPartitionInterface() .AND. AnyLayeredCurve()) CALL Fatal( &
            'LaminateSheetSteelReluctivity', 'A lamination sheet lies on a partition '// &
            'interface in a model with layered H-B bodies, where its parent across the '// &
            'interface is not at hand: partition the mesh so that the sheets lie inside '// &
            'partitions!')
      END IF
    END IF
    DO side = 1, 2
      IF (side == 1) Parent => Face % BoundaryInfo % Left
      IF (side == 2) Parent => Face % BoundaryInfo % Right
      IF (.NOT. ASSOCIATED(Parent)) CYCLE
      Material => GetMaterial(Parent, GotMaterial)
      IF (.NOT. GotMaterial) CYCLE
      IF (.NOT. ListCheckPresent(Material, 'H-B Curve')) CYCLE
      IF (GetString(Material, 'Laminate Stack Model', Stat) /= 'layered') CYCLE
      CALL GetLaminateLayeredParameters(Material, Thickness, Sigma, Fill, Normal, HasNormal)

      n = GetElementNOFNodes(Parent)
      nd = GetElementNOFDOFs(Parent, PSolver)
      np = n * PSolver % Def_Dofs(GetElementFamily(Parent), Parent % BodyId, 1)
      ALLOCATE(Pot(2,nd), Basis(nd), dBasisdx(nd,3), WBasis(nd,3), RotWBasis(nd,3))
      CALL GetElementNodes(Nodes, Parent, PSolver)
      CALL GetLocalSolution(Pot, UElement=Parent, USolver=PSolver)

      IP = GaussPoints(Parent, EdgeBasis=.TRUE., PReferenceElement=PiolaVersion, &
          EdgeBasisDegree=EdgeBasisDegree)
      DO t = 1, IP % n
        Stat = ElementInfo(Parent, Nodes, IP % U(t), IP % V(t), IP % W(t), detJ, Basis, &
            dBasisdx, EdgeBasis=WBasis, RotBasis=RotWBasis, USolver=PSolver)
        B = MATMUL(CMPLX(Pot(1,np+1:nd), Pot(2,np+1:nd), KIND=dp), RotWBasis(1:nd-np,:))
        BsAbs = MAX(LaminateMetalFluxDensity(B, Fill, Normal, HasNormal), 1.0d-8)
        NuSum = NuSum + ListGetFun(Material, 'h-b curve', BsAbs) / BsAbs * detJ * IP % s(t)
        Volume = Volume + detJ * IP % s(t)
      END DO
      DEALLOCATE(Pot, Basis, dBasisdx, WBasis, RotWBasis)
    END DO
    IF (Volume > 0.0_dp) THEN
      Found = .TRUE.
      NuSteel = NuSum / Volume
    END IF

  CONTAINS

    ! all nodes of the face are shared with another partition
    LOGICAL FUNCTION FaceOnPartitionInterface() RESULT(OnInterface)
      TYPE(Mesh_t), POINTER :: Mesh
      OnInterface = .FALSE.
      Mesh => PSolver % Mesh
      IF (.NOT. ASSOCIATED(Mesh % ParallelInfo % GInterface)) RETURN
      OnInterface = ALL(Mesh % ParallelInfo % GInterface(Face % NodeIndexes))
    END FUNCTION FaceOnPartitionInterface

    ! a 'layered' material with an H-B curve somewhere in the model
    LOGICAL FUNCTION AnyLayeredCurve() RESULT(HasCurve)
      INTEGER :: k
      LOGICAL :: GotIt
      HasCurve = .FALSE.
      DO k = 1, CurrentModel % NumberOfMaterials
        IF (ListGetString(CurrentModel % Materials(k) % Values, 'Laminate Stack Model', &
            GotIt) /= 'layered') CYCLE
        IF (ListCheckPresent(CurrentModel % Materials(k) % Values, 'H-B Curve')) HasCurve = .TRUE.
      END DO
    END FUNCTION AnyLayeredCurve
!------------------------------------------------------------------------------
  END FUNCTION LaminateSheetSteelReluctivity
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Incremental reluctivity dH/dB of the H-B curve of a material at |B| = Babs.
!------------------------------------------------------------------------------
  FUNCTION CurveIncrementalReluctivity(Material, Babs) RESULT(dHdB)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp), INTENT(IN) :: Babs
    REAL(KIND=dp) :: dHdB
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: H
!------------------------------------------------------------------------------
    H = ListGetFun(Material, 'h-b curve', Babs, dFdx=dHdB)
    IF (dHdB <= 0.0_dp) THEN
      WRITE(Message,'(A,ES12.4,A)') 'The H-B curve does not rise at |B| =', Babs, ' T!'
      CALL Fatal('CurveIncrementalReluctivity', Message)
    END IF
!------------------------------------------------------------------------------
  END FUNCTION CurveIncrementalReluctivity
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Projector onto the lamination plane of a 'layered' stack, I - n n^T; the
!> identity without a normal, where every direction is in-plane.
!------------------------------------------------------------------------------
  FUNCTION LaminateInPlaneProjector(Normal, HasNormal) RESULT(P)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN) :: Normal(3)
    LOGICAL, INTENT(IN) :: HasNormal
    REAL(KIND=dp) :: P(3,3)
!------------------------------------------------------------------------------
    INTEGER :: i, j
!------------------------------------------------------------------------------
    P = 0.0_dp
    DO i = 1, 3
      P(i,i) = 1.0_dp
    END DO
    IF (HasNormal) THEN
      DO i = 1, 3
        DO j = 1, 3
          P(i,j) = P(i,j) - Normal(i) * Normal(j)
        END DO
      END DO
    END IF
!------------------------------------------------------------------------------
  END FUNCTION LaminateInPlaneProjector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> In-plane energy integral of the static 'layered' law of an H-B curve,
!> I_par(b) = int_0^b beta nu_par(beta) dbeta, nu_par = H/(F beta + (1-F) mu0 H),
!> b >= 0, by five-point Gauss on each interval of the curve. Its values at the
!> knots above zero are tabulated once per curve and stacking factor; a full
!> cache gives its oldest table up for a new one.
!------------------------------------------------------------------------------
  FUNCTION LaminateInPlaneEnergyIntegral(Curve, StackingFactor, b) RESULT(Ipar)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueListEntry_t), POINTER :: Curve
    REAL(KIND=dp), INTENT(IN) :: StackingFactor, b
    REAL(KIND=dp) :: Ipar
!------------------------------------------------------------------------------
    INTEGER :: it, j, n
!------------------------------------------------------------------------------
    n = SIZE(Curve % TValues)
    ! the largest knot in (0, b], if any: the integral up to it is tabulated
    j = n
    DO WHILE (j >= 1)
      IF (Curve % TValues(j) <= b) EXIT
      j = j - 1
    END DO
    ! the table is read where it is found or built, as another call may replace it
    !$OMP CRITICAL(LaminateEnergyTables)
    it = LaminateEnergyTableIndex()
    IF (j < 1) THEN
      Ipar = SegmentIntegral(0.0_dp, b)
    ELSE IF (Curve % TValues(j) <= 0.0_dp) THEN
      Ipar = SegmentIntegral(0.0_dp, b)
    ELSE
      Ipar = LamEnergyTables(it) % Cumulative(j) + SegmentIntegral(Curve % TValues(j), b)
    END IF
    !$OMP END CRITICAL(LaminateEnergyTables)

  CONTAINS

    ! the table of this curve (its knots, values and cubic coefficients: an
    ! ordinary and a monotone cubic differ in these only) and stacking factor,
    ! built at the first call: the cumulative integral from zero at each knot
    ! above zero (zero below). Called inside the critical section.
    FUNCTION LaminateEnergyTableIndex() RESULT(it)
      INTEGER :: it, k
      REAL(KIND=dp) :: lo
      LOGICAL :: Cubic
      Cubic = ASSOCIATED(Curve % CubicCoeff)
      DO it = 1, NoLamEnergyTables
        IF (LamEnergyTables(it) % StackingFactor /= StackingFactor) CYCLE
        IF (LamEnergyTables(it) % Cubic .NEQV. Cubic) CYCLE
        IF (SIZE(LamEnergyTables(it) % Knots) /= n) CYCLE
        IF (ANY(LamEnergyTables(it) % Knots /= Curve % TValues)) CYCLE
        IF (ANY(LamEnergyTables(it) % Values /= Curve % FValues(1,1,:))) CYCLE
        IF (Cubic) THEN
          IF (SIZE(LamEnergyTables(it) % CubicCoeff) /= SIZE(Curve % CubicCoeff)) CYCLE
          IF (ANY(LamEnergyTables(it) % CubicCoeff /= Curve % CubicCoeff)) CYCLE
        END IF
        RETURN
      END DO
      IF (NoLamEnergyTables < MaxLaminateEnergyTables) THEN
        NoLamEnergyTables = NoLamEnergyTables + 1
        it = NoLamEnergyTables
      ELSE
        it = MOD(LastLamEnergyTable, MaxLaminateEnergyTables) + 1
      END IF
      LastLamEnergyTable = it
      IF (ALLOCATED(LamEnergyTables(it) % Knots)) DEALLOCATE(LamEnergyTables(it) % Knots, &
          LamEnergyTables(it) % Values, LamEnergyTables(it) % Cumulative)
      IF (ALLOCATED(LamEnergyTables(it) % CubicCoeff)) DEALLOCATE(LamEnergyTables(it) % CubicCoeff)
      ALLOCATE(LamEnergyTables(it) % Knots(n), LamEnergyTables(it) % Values(n), &
          LamEnergyTables(it) % Cumulative(n))
      LamEnergyTables(it) % Knots = Curve % TValues
      LamEnergyTables(it) % Values = Curve % FValues(1,1,:)
      IF (Cubic) THEN
        ALLOCATE(LamEnergyTables(it) % CubicCoeff(SIZE(Curve % CubicCoeff)))
        LamEnergyTables(it) % CubicCoeff = Curve % CubicCoeff
      END IF
      LamEnergyTables(it) % Cubic = Cubic
      LamEnergyTables(it) % StackingFactor = StackingFactor
      LamEnergyTables(it) % Cumulative = 0.0_dp
      lo = 0.0_dp
      DO k = 1, n
        IF (Curve % TValues(k) <= 0.0_dp) CYCLE
        LamEnergyTables(it) % Cumulative(k) = SegmentIntegral(lo, Curve % TValues(k))
        IF (lo > 0.0_dp) LamEnergyTables(it) % Cumulative(k) = &
            LamEnergyTables(it) % Cumulative(k) + LamEnergyTables(it) % Cumulative(k-1)
        lo = Curve % TValues(k)
      END DO
    END FUNCTION LaminateEnergyTableIndex

    ! int_lo^hi beta nu_par(beta) dbeta by five-point Gauss
    FUNCTION SegmentIntegral(lo, hi) RESULT(s)
      REAL(KIND=dp), INTENT(IN) :: lo, hi
      REAL(KIND=dp) :: s, beta, h, mu0
      REAL(KIND=dp), PARAMETER :: GaussU(5) = [-0.906179845938664_dp, -0.538469310105683_dp, &
          0.0_dp, 0.538469310105683_dp, 0.906179845938664_dp]
      REAL(KIND=dp), PARAMETER :: GaussW(5) = [0.236926885056189_dp, 0.478628670499366_dp, &
          0.568888888888889_dp, 0.478628670499366_dp, 0.236926885056189_dp]
      INTEGER :: k
      s = 0.0_dp
      IF (hi <= lo) RETURN
      mu0 = VacuumPermeability()
      DO k = 1, 5
        beta = 0.5_dp * (lo + hi) + 0.5_dp * (hi - lo) * GaussU(k)
        h = InterpolateCurve(Curve % TValues, Curve % FValues(1,1,:), beta, Curve % CubicCoeff)
        s = s + GaussW(k) * beta * h / (StackingFactor * beta + (1.0_dp - StackingFactor) * mu0 * h)
      END DO
      s = 0.5_dp * (hi - lo) * s
    END FUNCTION SegmentIntegral
!------------------------------------------------------------------------------
  END FUNCTION LaminateInPlaneEnergyIntegral
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Field energy density of the static 'layered' stack law at B (real, transient).
!> The law depends on B through b = |B_s| only, so its energy on the straight
!> path from 0 is |B_par|^2/b^2 I_par(b) + |B_perp|^2/b^2 I_perp(b) with
!> I(b) = int_0^b beta nu(beta) dbeta of the in-plane and the normal stack
!> reluctivity: I_perp = F int H dB + (1-F) b^2/(2 mu0) is the curve integral,
!> I_par is integrated along the curve (LaminateInPlaneEnergyIntegral). With
!> B in the lamination plane or along the normal this is the potential of the
!> law. NuSteel is the steel's reluctivity at b; with HasCurve the H-B curve of
!> the material is integrated, else the steel is linear.
!------------------------------------------------------------------------------
  FUNCTION LaminateLayeredEnergy(Material, B, StackingFactor, Normal, HasNormal, &
      HasCurve, NuSteel) RESULT(W)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp), INTENT(IN) :: B(3), StackingFactor, Normal(3), NuSteel
    LOGICAL, INTENT(IN) :: HasNormal, HasCurve
    REAL(KIND=dp) :: W
!------------------------------------------------------------------------------
    TYPE(ValueListEntry_t), POINTER :: Curve
    COMPLEX(KIND=dp) :: Nu(3,3)
    REAL(KIND=dp) :: Bperp(3), Bpar(3), bs, mu0, CurveEnergy
    LOGICAL :: Found
!------------------------------------------------------------------------------
    IF (.NOT. HasCurve) THEN
      Nu = LaminateLayeredReluctivity(0.0_dp, CMPLX(NuSteel, 0.0_dp, KIND=dp), 0.0_dp, &
          1.0_dp, StackingFactor, Normal, HasNormal)
      W = 0.5_dp * SUM(B * MATMUL(REAL(Nu), B))
      RETURN
    END IF

    W = 0.0_dp
    bs = LaminateMetalFluxDensity(CMPLX(B, 0.0_dp, KIND=dp), StackingFactor, Normal, HasNormal)
    IF (bs <= 0.0_dp) RETURN
    mu0 = VacuumPermeability()
    Bperp = 0.0_dp
    IF (HasNormal) Bperp = Normal * SUM(Normal * B)
    Bpar = B - Bperp

    Curve => ListFind(Material, 'h-b curve')
    CurveEnergy = IntegrateCurve(Curve % TValues, Curve % FValues(1,1,:), Curve % CubicCoeff, &
        0.0_dp, bs, Found=Found)
    W = SUM(Bpar**2) / bs**2 * LaminateInPlaneEnergyIntegral(Curve, StackingFactor, bs) + &
        SUM(Bperp**2) / bs**2 * (StackingFactor * CurveEnergy + &
        (1.0_dp - StackingFactor) * bs**2 / (2 * mu0))
!------------------------------------------------------------------------------
  END FUNCTION LaminateLayeredEnergy
!------------------------------------------------------------------------------

!> Get a real-valued reluctivity tensor. This subroutine seeks values which
!> are strictly given as reluctivity (giving the permeability is not an option
!> here).
!------------------------------------------------------------------------------
  SUBROUTINE GetReluctivityTensorR(Material, Acoef, n, Found)
!-------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER, INTENT(IN) :: Material
    REAL(KIND=dp), POINTER, INTENT(OUT) :: Acoef(:,:,:)
    INTEGER, INTENT(IN) :: n
    LOGICAL, INTENT(OUT) :: Found
!------------------------------------------------------------------------------
    REAL(KIND=dp), SAVE :: nu_vacuum
    LOGICAL, SAVE :: FirstTime = .TRUE.
!------------------------------------------------------------------------------

    IF ( FirstTime ) THEN
      nu_vacuum = GetConstReal( CurrentModel % Constants, &
              'Permeability of Vacuum', Found )
      IF (.NOT. Found ) THEN
        nu_vacuum = 1.0d0/(PI * 4.0d-7)
      ELSE
        nu_vacuum =  1.0d0/nu_vacuum
      END IF
      FirstTime = .FALSE.
    END IF

    CALL GetRealArray( Material, Acoef, 'Reluctivity', Found )

    IF (.NOT. Found) THEN
      CALL GetRealArray( Material, Acoef, 'Relative Reluctivity', Found )
      IF (Found) Acoef = nu_vacuum * Acoef
    END IF
!-------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivityTensorR
!-------------------------------------------------------------------------------

!> Get a complex-valued reluctivity tensor. This subroutine seeks values which
!> are strictly given as reluctivity (giving the permeability is not an option
!> here).
!------------------------------------------------------------------------------
  SUBROUTINE GetReluctivityTensorC(Material, Acoef, n, Found)
!-------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER, INTENT(IN) :: Material
    COMPLEX(KIND=dp), POINTER, INTENT(OUT) :: Acoef(:,:,:)
    INTEGER, INTENT(IN) :: n                                      ! An inactive variable
    LOGICAL, INTENT(OUT) :: Found                                 
!-------------------------------------------------------------------------------
    LOGICAL :: Found_im
    REAL(KIND=dp), POINTER :: work(:,:,:) => NULL()
    INTEGER :: n1, n2, n3
    REAL(KIND=dp), SAVE :: nu_vacuum
    LOGICAL, SAVE :: FirstTime = .TRUE.
!------------------------------------------------------------------------------

    IF ( FirstTime ) THEN
      nu_vacuum = GetConstReal( CurrentModel % Constants, &
              'Permeability of Vacuum', Found )
      IF (.NOT. Found ) THEN
        nu_vacuum = 1.0d0/(PI * 4.0d-7)
      ELSE
        nu_vacuum =  1.0d0/nu_vacuum
      END IF
      FirstTime = .FALSE.
    END IF
    
    IF (ASSOCIATED(Acoef)) DEALLOCATE(Acoef)

    CALL GetRealArray( Material, work, 'Reluctivity', Found )
    IF (.NOT. Found) THEN
      CALL GetRealArray( Material, work, 'Relative Reluctivity', Found )
      IF (Found) work = nu_vacuum * work
    END IF

    IF (Found) THEN
      n1 = SIZE(work,1)
      n2 = SIZE(work,2)
      n3 = SIZE(work,3)
      ALLOCATE(Acoef(n1, n2, n3))
      Acoef(:,:,:) = CMPLX(work(:,:,:), 0.0d0, kind=dp)
    END IF

    CALL GetRealArray( Material, work, 'Reluctivity im', Found_im )
    IF (.NOT. Found_im) THEN
      CALL GetRealArray( Material, work, 'Relative Reluctivity im', Found_im )
      IF (Found_im) work = nu_vacuum * work
    END IF
    
    IF (Found_im) THEN
      n1 = SIZE(work,1)
      n2 = SIZE(work,2)
      n3 = SIZE(work,3)
      IF (.NOT. ASSOCIATED(Acoef)) THEN
        ALLOCATE(Acoef(n1, n2, n3))
        Acoef(:,:,:) = CMPLX(0.0d0, work(:,:,:), kind=dp)
      ELSE
        IF (SIZE(Acoef,1) /= n1 .OR. SIZE(Acoef,2) /= n2 .OR.  SIZE(Acoef,3) /= n3) &
            CALL Fatal('GetReluctivityTensorC', 'Reluctivity and Reluctivity im of different size')
        Acoef(1:n1,1:n2,1:n3) = CMPLX(REAL(Acoef(1:n1,1:n2,1:n3)), work(1:n1,1:n2,1:n3), kind=dp)
      END IF
    END IF
    Found = Found .OR. Found_im

    IF (ASSOCIATED(work)) DEALLOCATE(work)
!-------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivityTensorC
!-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !> Packs rows associated with edge dofs from constraint matrix and
  !> adds ColOffset to column indices of (1,1) block.
  !> If CM % ConstrainedDOF(i) is true then the ith row is treted empty.
  !! @param ColOffset add to column index of that correspond to nodal dofs
  !!                   by this amount 
  !-------------------------------------------------------------------------------
  SUBROUTINE PackEdgeRows(CM, Model, ColOffset)
    !-------------------------------------------------------------------------------
    TYPE(Matrix_t), INTENT(INOUT), POINTER :: CM 
    TYPE(Model_t), INTENT(IN) :: Model
    INTEGER, OPTIONAL :: ColOffset
    !-------------------------------------------------------------------------------
    LOGICAL, ALLOCATABLE :: Emptyrow(:)
    INTEGER :: ColOffsetF
    REAL(KIND=dp), POINTER CONTIG :: V(:), RHS(:)
    INTEGER, POINTER CONTIG :: R(:), C(:), InvPerm(:), Perm(:)
    INTEGER :: i, j, numempty, k, fileind

    ALLOCATE(EmptyRow(CM % NumberOfRows))
    EmptyRow = .TRUE.

    IF (.NOT. PRESENT(ColOffset)) THEN
      ColOffsetF = CM % NumberOfRows
    ELSE
      ColOffsetF = ColOffset
    END IF

    DO j = 1, CM % NumberOfRows
      DO i = CM % Rows(j), CM % Rows(j+1)-1
        IF (CM % Cols(i) <= Model % NumberOfNodes) &
             CM % Cols(i) = CM % Cols(i) + ColOffsetF
      END DO
    END DO

    numempty = 0
    ROW_LOOP: DO i = 1, CM % NumberOfRows
      ! If CM % ConstrainedDOF(i) is true, then this must be a row corresponding to
      ! edge dof so it must be zero. 
      IF ( CM % ConstrainedDOF(i) ) THEN
        Emptyrow(i) = .TRUE.
      ELSE ! Otherwise the row might correspond with dirichlet scalar dof
        ! Such row must not be packed but diagonal entry of the final matrix
        ! must be set to 1. (Done later?)
        COL_LOOP : DO j = CM % Rows(i), CM % Rows(i+1)-1
          IF (abs(CM % values(j)) >= 1e-12_dp) THEN
            EmptyRow(i) = .FALSE.
            EXIT COL_LOOP
          END IF
        END DO COL_LOOP
      END IF

      IF (EmptyRow(i)) THEN
        CM % RHS(i) = 0.0_dp
        IF (CM % Constraineddof(i)) numempty = numempty + 1
      END IF
    END DO ROW_LOOP


    ! Pack empty rows
    ALLOCATE(R(CM % NumberOfRows-numempty+1))
    ALLOCATE(InvPerm(CM % NumberOfRows-numempty))
    ALLOCATE(Perm(CM % NumberOfRows-numempty))
    ALLOCATE(RHS(CM % NumberOfRows-numempty))
    RHS(:) = 0.0_dp
    k = 2
    R(1) = 1
    DO i = 2, CM % NumberOfRows+1
      IF(.not. emptyrow(i-1)) THEN
        IF(ASSOCIATED(CM%InvPerm)) InvPerm(k-1) = CM % InvPerm(i-1)  !< TODO: This might be nonsense
        IF(ASSOCIATED(CM%Perm)) Perm(k-1) = CM % Perm(i-1)
        R(k) = R(k-1) + CM % Rows(i) - CM % rows(i-1)

        k = k + 1
      END IF
    END DO
    k = size(R)
    ALLOCATE(C(R(k)-1), V(R(k)-1))

    k = 1
    DO i = 2, CM % NumberOfRows+1
      IF (.not. EmptyRow(i-1)) THEN
        DO j = 0, R(k+1)-R(k)-1
          C(R(k)+j) = CM % cols(CM % rows(i-1)+j)
          V(R(k)+j) = CM % Values(CM % rows(i-1)+j)
        END DO
        k = k+1
      END IF
    END DO

    deallocate(CM % Rows)
    deallocate(CM % Cols)
    deallocate(CM % Values)
    deallocate(cm % rhs)
    if(ASSOCIATED(CM % InvPerm)) THEN
      DEALLOCATE(CM % InvPerm)
      CM % InvPerm => InvPerm
    END IF
    IF(ASSOCIATED(CM % Perm)) THEN 
      DEALLOCATE(CM % Perm)
      CM % Perm => Perm
    END IF

    CM % rows => R
    CM % cols => C
    CM % values => V
    CM % RHS => RHS
    CM % NumberOfRows = size(R,1)-1

  !-------------------------------------------------------------------------------
  END SUBROUTINE PackEdgeRows
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  SUBROUTINE CommunicateCondReg(Solver, Mesh, CondReg)
  !-------------------------------------------------------------------------------
     TYPE(Mesh_t) :: Mesh
     LOGICAL :: CondReg(:)
     TYPE(Solver_t) :: Solver

     INTEGER, ALLOCATABLE :: s_e(:,:), r_e(:)
     INTEGER :: i,j,k,l,n,ii(ParEnv % PEs), ierr, status(MPI_STATUS_SIZE)

     IF( ParEnv % PEs<=1 ) RETURN

     ALLOCATE(s_e(Mesh % NumberOfNodes, ParEnv % PEs), &
               r_e(Mesh % NUmberOfNodes) )
     ii = 0
     DO i=1,Mesh % NumberOfNodes
       IF(.NOT.CondReg(i) .AND. Mesh % ParallelInfo % GInterface(i) ) THEN
          DO j=1,SIZE(Mesh % ParallelInfo % Neighbourlist(i) % Neighbours)
            k = Mesh % ParallelInfo % Neighbourlist(i) % Neighbours(j)
            IF ( k== ParEnv % MyPE ) CYCLE
            k = k + 1
            ii(k) = ii(k) + 1
            s_e(ii(k),k) = Mesh % ParallelInfo % GlobalDOFs(i)
          END DO
       END IF
     END DO

     DO i=0, ParEnv % PEs-1
       IF(i==ParEnv % myPE) CYCLE
       k = i+1
       CALL MPI_BSEND( ii(k),1,MPI_INTEGER,i,110,Solver % matrix % comm,ierr )
       IF( ii(k) > 0 ) THEN
         CALL MPI_BSEND( s_e(1:ii(k),k),ii(k),MPI_INTEGER,i,111,solver % matrix % comm,ierr )
       END IF
     END DO

     DO i=0, ParENv % PEs-1
       IF(i==ParEnv % myPE) CYCLE
       CALL MPI_RECV( n,1,MPI_INTEGER,i,110,Solver % Matrix % Comm, status,ierr )
       IF ( n>0 ) THEN
         CALL MPI_RECV( r_e,n,MPI_INTEGER,i,111,Solver % Matrix % Comm, status,ierr )
         DO j=1,n
           k = SearchNode( Mesh % ParallelInfo, r_e(j) )
           IF ( k>0 ) CondReg(k) = .FALSE.
         END DO
       END IF
     END DO
     DEALLOCATE(s_e, r_e)
  !-------------------------------------------------------------------------------
  END SUBROUTINE CommunicateCondReg
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE RecvDoneNodesAndEdges(Solver,Mesh,Done,TreeEdges)
  !-------------------------------------------------------------------------------
     TYPE(Solver_t):: Solver
     TYPE(Mesh_t) :: Mesh
     LOGICAL :: Done(:), TreeEdges(:)
  
     INTEGER, ALLOCATABLE :: s_e(:,:), r_e(:), iperm(:)
     INTEGER :: i,j,k,l,n,ii(ParEnv % PEs), ierr, status(MPI_STATUS_SIZE)

     IF(ParEnv % myPE <= 0) RETURN

     ALLOCATE(r_e(Mesh % NumberOfEdges), iperm(SIZE(Solver % Variable % Perm)))
     iperm = 0
     DO i=1,SIZE(Solver % Variable % Perm)
       IF ( Solver % Variable % Perm(i) > 0 ) &
         iperm(Solver % Variable % Perm(i)) = i
     END DO

     DO i=0, ParEnv % MyPE-1
       CALL MPI_RECV( n,1,MPI_INTEGER,i,112,Solver % Matrix % Comm, status,ierr )
       IF ( n>0 ) THEN
         CALL MPI_RECV( r_e,n,MPI_INTEGER,i,113,Solver % Matrix % Comm, status,ierr )
         DO j=1,n
           k = SearchNode( Solver % Matrix % ParallelInfo, r_e(j), Order=Solver % Variable % Perm )
           k = iperm(k) - Mesh % NumberOfNodes
           IF ( k>0 .AND. k<=SIZE(TreeEdges) ) THEN
             TreeEdges(k) = .TRUE.
           END IF
         END DO
       END IF
     END DO

     DO i=0, ParENv % myPE-1
       CALL MPI_RECV( n,1,MPI_INTEGER,i,114,Solver % Matrix % Comm, status,ierr )
       IF ( n>0 ) THEN
         CALL MPI_RECV( r_e,n,MPI_INTEGER,i,115,Solver % Matrix % Comm, status,ierr )
         DO j=1,n
           k = SearchNode( Mesh % ParallelInfo, r_e(j) )
           IF ( k>0 ) Done(k) = .TRUE.
         END DO
       END IF
     END DO
  !-------------------------------------------------------------------------------
  END SUBROUTINE RecvDoneNodesAndEdges
  !-------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  SUBROUTINE SendDoneNodesAndEdges(Solver,Mesh,Done,TreeEdges)
  !-------------------------------------------------------------------------------
     TYPE(Solver_t):: Solver
     TYPE(Mesh_t) :: Mesh
     LOGICAL :: Done(:), TreeEdges(:)
  
     INTEGER, ALLOCATABLE :: s_e(:,:), r_e(:), iperm(:)
     INTEGER :: i,j,k,l,n,ii(ParEnv % PEs), ierr, status(MPI_STATUS_SIZE)

     IF(Parenv % myPE < ParEnv % PEs-1) THEN

      ALLOCATE(s_e(Mesh % NumberOfEdges, ParEnv % PEs) )

      ii = 0
      DO i=1,Mesh % NumberOfedges
        IF ( TreeEdges(i) .AND. Mesh % ParallelInfo % EdgeInterface(i) ) THEN
          DO j=1,SIZE(Mesh % ParallelInfo % EdgeNeighbourlist(i) % Neighbours)
            k = Mesh % ParallelInfo % EdgeNeighbourlist(i) % Neighbours(j)
            IF ( k>ParEnv % myPE ) THEN
               k = k + 1
               ii(k) = ii(k) +1
               s_e(ii(k),k) = Solver % Matrix % ParallelInfo % GlobalDOFs( & 
                    Solver % Variable % Perm(Mesh % NumberOfNodes+i) )
            END IF
          END DO
        END IF
      END DO
 
      DO i=Parenv % mype+1,Parenv % PEs-1
        k = i+1
        CALL MPI_BSEND( ii(k),1,MPI_INTEGER,i,112,Solver % matrix % comm,ierr )
        IF( ii(k) > 0 ) THEN
          CALL MPI_BSEND( s_e(1:ii(k),k),ii(k),MPI_INTEGER,i,113,solver % matrix % comm,ierr )
        END IF
      END DO

      ii = 0
      DO i=1,Mesh % NumberOfNodes
        IF ( Done(i) .AND. Mesh % ParallelInfo % GInterface(i) ) THEN
          DO j=1,SIZE(Mesh % ParallelInfo % Neighbourlist(i) % Neighbours)
            k = Mesh % ParallelInfo % Neighbourlist(i) % Neighbours(j)
            IF ( k>ParEnv % myPE ) THEN
               k = k + 1
               ii(k) = ii(k) +1
               s_e(ii(k),k) = Mesh % ParallelInfo % GlobalDOFs(i)
            END IF
          END DO
        END IF
      END DO

      DO i=Parenv % mype+1,Parenv % PEs-1
        k = i+1
        CALL MPI_BSEND( ii(k),1,MPI_INTEGER,i,114,Solver % matrix % comm,ierr )
        IF( ii(k) > 0 ) THEN
          CALL MPI_BSEND( s_e(1:ii(k),k),ii(k),MPI_INTEGER,i,115,solver % matrix % comm,ierr )
        END IF
      END DO
    END IF

    CALL SparIterBarrier()
  !-------------------------------------------------------------------------------
  END SUBROUTINE SendDoneNodesAndEdges
  !-------------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------------
  ! Mark nodes that are on outer boundary using face elements and node parmutation.
  !-------------------------------------------------------------------------------
  SUBROUTINE MarkOuterNodes(Mesh,Perm,SurfaceNodes,SurfacePerm,EnsureBC) 

    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: Perm(:),SurfacePerm(:)
    INTEGER :: SurfaceNodes
    LOGICAL :: EnsureBC
    
    INTEGER :: snodes0, snodes, i,t,n,ActParents,ParParents
    TYPE(Element_t), POINTER :: Element, P1, P2, P
    LOGICAL, ALLOCATABLE :: BcNode(:)
    
    CALL Info('MarkOuterNodes','Marking outer nodes on outer boundary',Level=8)
    
    SurfaceNodes = 0

    IF( Mesh % NumberOfFaces == 0 ) THEN
      CALL Fatal('MarkOuterNodes','The faces are not created!')
    END IF
    
    n = Mesh % NumberOfNodes
    IF(.NOT. ASSOCIATED( SurfacePerm ) ) THEN
      ALLOCATE( SurfacePerm( n ) )
    END IF
    SurfacePerm = 0
    
       
    DO t=1, Mesh % NumberOfFaces 
      
      Element => Mesh % Faces(t)         
      
      IF( ParEnv % PEs > 1 ) THEN
        ! Don't set BCs on partition interfaces
        IF( Mesh % ParallelInfo % FaceInterface(t) ) CYCLE
      END IF
      
      P1 => Element % BoundaryInfo % Left
      P2 => Element % BoundaryInfo % Right
      
      ActParents = 0
      ParParents = 0

      IF( ASSOCIATED( P1 ) ) THEN
        IF (ALL(Perm(P1 % NodeIndexes)>0)) THEN
          ActParents = ActParents + 1
          IF( P1 % PartIndex == ParEnv % MyPe ) ParParents = ParParents + 1
        ELSE
          NULLIFY( P1 )
        END IF
      END IF
      IF( ASSOCIATED( P2 ) ) THEN
        IF (ALL(Perm(P2 % NodeIndexes)>0)) THEN
          ActParents = ActParents + 1
          IF( P2 % PartIndex == ParEnv % MyPe ) ParParents = ParParents + 1         
        ELSE
          NULLIFY( P2 )
        END IF
      END IF

      ! We have either none or both parents as active.
      ! The BCs will be set only to outer boundaries of the domain. 
      IF( ActParents /= 1 ) CYCLE
      
      ! The one parent is not a true one!
      ! This can happen when we have halo elements. 
      IF( ParEnv % PEs > 0 .AND. ParParents == 0 ) CYCLE
      
      SurfacePerm(Element % NodeIndexes) = 1
    END DO

    IF( EnsureBC ) THEN
      ALLOCATE( BcNode(n) )
      BcNode = .FALSE.

      DO t=Mesh % NumberOfBulkElements +1, &
          Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements      
        Element => Mesh % Elements(t)         
        BcNode(Element % NodeIndexes) = .TRUE.
      END DO

      snodes0 = COUNT( SurfacePerm > 0 )
      snodes0 = ParallelReduction(snodes0) 

      !DO i=1,n
      !  IF( SurfacePerm(i) > 0 .AND. .NOT. BcNode(i) ) THEN
      !    PRINT *,'node:',ParEnv % MyPe, i, Mesh % Nodes % x(i), Mesh % Nodes % y(i), Mesh % Nodes % z(i)
      !  END IF
      !END DO

      WHERE( .NOT. BcNode ) SurfacePerm = 0
      DEALLOCATE( BcNode ) 
    END IF
      
    ! Create numbering for the surface nodes
    snodes = 0
    DO i=1,n
      IF( SurfacePerm(i) > 0 ) THEN
        snodes = snodes + 1
        SurfacePerm(i) = snodes
      END IF
    END DO     
    
    snodes = ParallelReduction(snodes) 
    CALL Info('MarkOuterNodes','Total number of surface nodes: '//I2S(snodes),Level=6)

    IF( EnsureBC ) THEN
      IF( snodes0 > snodes ) THEN
        CALL Info('MarkOuterNodes','Removed number of surface nodes not at BCs: '&
            //I2S(snodes0-snodes),Level=6)
      END IF
    END IF

    SurfaceNodes = snodes
    
  END SUBROUTINE MarkOuterNodes

!------------------------------------------------------------------------------
  SUBROUTINE GaugeTree(Solver,Mesh,TreeEdges,FluxCount,FluxMap,Transient)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Solver_t) :: Solver
    INTEGER :: FluxCount, FluxMap(:)
    LOGICAL :: Transient
    TYPE(Mesh_t) :: Mesh
    LOGICAL, ALLOCATABLE :: TreeEdges(:)

    TYPE(ListMatrixEntry_t), POINTER :: Aentry
    TYPE(ListMatrix_t), POINTER :: Alist(:)
    INTEGER :: i,j,k,l,n,Start
    LOGICAL, ALLOCATABLE :: Done(:), CondReg(:)
    TYPE(ValueList_t), POINTER :: BC
    REAL(KIND=dp) :: Cond1, NodalCond(27)
    TYPE(Element_t), POINTER :: Edge, Boundary, Element

    INTEGER, ALLOCATABLE :: r_e(:), s_e(:,:), iperm(:)
    LOGICAL :: Found
    INTEGER :: ssz, status(MPI_STATUS_SIZE), ierr, ii(ParEnv % PEs)
!------------------------------------------------------------------------------

    IF( ALLOCATED( TreeEdges ) ) THEN
      CALL Info('WhitneyAVSolver','Gauge tree already created',Level=15)
      RETURN
    END IF
      
    ALLOCATE(TreeEdges(Mesh % NumberOfEdges))
    TreeEdges = .FALSE.

    n = Mesh % NumberOfNodes
    ALLOCATE(Done(n)); Done=.FALSE.

    ! Skip Dirichlet BCs in terms of A:
    ! ---------------------------------
    DO i=1,Mesh % NumberOfBoundaryElements
      Boundary => GetBoundaryElement(i)

      SELECT CASE(GetElementFamily())
      CASE(1)
        CYCLE
      CASE(2)
        k = GetBoundaryEdgeIndex(Boundary,1); Element => Mesh % Edges(k)
      CASE(3,4)
        k = GetBoundaryFaceIndex(Boundary)  ; Element => Mesh % Faces(k)
      END SELECT
      IF (.NOT. ActiveBoundaryElement(Element)) CYCLE

      BC => GetBC()
      IF (.NOT.ASSOCIATED(BC)) CYCLE
      IF (.NOT.( ListCheckPresent(BC, 'Mortar BC') .OR. ListCheckPresent( BC, &
                 TRIM(Solver % Variable % Name)//' {e}'))) CYCLE
 
      Done(Element % NodeIndexes) = .TRUE.
    END DO

    IF( Transient ) THEN
      IF ( GetLogical( GetSolverParams(), 'Gauge Tree Skip Conducting Regions', Found) ) THEN
        ! Skip conducting regions:
        ! -------------------------
        ALLOCATE(CondReg(Mesh % NumberOfNodes))
        condReg = .TRUE.
        DO i=1,GetNOFActive()
          Element => GetActiveElement(i)
          n = Element % Type % NumberOfNodes
          NodalCond(1:n) = GetReal(GetMaterial(), 'Electric Conductivity',Found)
          cond1 = MAXVAL(NodalCond(1:n))
          IF (cond1 < EPSILON(cond1) ) condReg(Element % NodeIndexes) = .FALSE.
        END DO

        CALL CommunicateCondReg(Solver,Mesh,CondReg)

        Done = Done.OR.CondReg
        DEALLOCATE(CondReg)
      END IF
    END IF

    ! 
    ! Skip Dirichlet BCs in terms of B:
    ! ---------------------------------
    DO i=1,FluxCount
      j = FluxMap(i)
      IF ( Solver % Variable % Perm(j+n)<=0 ) CYCLE
      Edge => Mesh % Edges(j)
      Done(Edge % NodeIndexes)=.TRUE.
    END DO

    ! 
    ! already set:
    ! ------------

    CALL RecvDoneNodesAndEdges(Solver,Mesh,Done,TreeEdges)

    ! node -> edge list
    ! -----------------
    Alist => NULL()
    n = Mesh % NumberOfNodes
    DO i=1,Mesh % NumberOfEdges
      Edge => Mesh % Edges(i)
      IF ( Solver % Variable % Perm(i+n)<=0 ) CYCLE
      DO j=1,Edge % TYPE % NumberOfNodes
        k=Edge % NodeIndexes(j)
        Aentry=>List_GetMatrixIndex(Alist,k,i)
      END DO
    END DO

    !
    ! generate the tree for all (perhaps disconnected) parts:
    ! -------------------------------------------------------
    DO WHILE(.NOT.ALL(Done))
      DO Start=1,n
        IF (.NOT. Done(Start)) EXIT
      END DO
      CALL DepthFirstSearch(Alist,Done,Start)
    END DO
    CALL List_FreeMatrix(SIZE(Alist),Alist)

    CALL SendDoneNodesAndEdges(Solver,Mesh,Done,TreeEdges)
    DEALLOCATE(Done)

CONTAINS

!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE DepthFirstSearch(Alist,done,i)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ListMatrix_t) :: Alist(:)
    INTEGER :: i
    LOGICAL :: Done(:)
!------------------------------------------------------------------------------
    TYPE(ListMatrixEntry_t), POINTER :: Aentry
    INTEGER :: j,k,l,n
    TYPE(Element_t), POINTER :: Edge
!------------------------------------------------------------------------------

    ! To give better matrix conditioning some directional heuristics
    ! could be added,e.g. select the order of going through the nodes
    ! edge list here:

    Done(i) = .TRUE.

    Aentry => Alist(i) % Head
    DO WHILE(ASSOCIATED(Aentry))
      k = Aentry % Index
      Aentry => Aentry % Next

      Edge => Mesh % Edges(k)
      IF (ALL(Done(Edge % NodeIndexes))) CYCLE

      IF ( .NOT. TreeEdges(k)) CALL SetDOFToValue(Solver,k,0._dp)
      TreeEdges(k)=.TRUE.
      DO l=1,2
        n = Edge % NodeIndexes(l)
        IF (.NOT. Done(n)) CALL DepthFirstSearch(Alist,done,n)
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE DepthFirstSearch
!------------------------------------------------------------------------------


    
!------------------------------------------------------------------------------
  END SUBROUTINE GaugeTree
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE GaugeTreeFluxBC(Solver,Mesh,TreeEdges,BasicCycles,FluxCount,FluxMap)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ListMatrix_t), POINTER :: BasicCycles(:)
    INTEGER :: FluxCount, FluxMap(:)
    TYPE(Solver_t) :: Solver
    TYPE(Mesh_t) :: Mesh
    LOGICAL, ALLOCATABLE :: TreeEdges(:)
!------------------------------------------------------------------------------
!   TYPE(Mesh_t) :: Mesh
!   LOGICAL, ALLOCATABLE :: TreeEdges(:)
!------------------------------------------------------------------------------
    TYPE(ListMatrixEntry_t), POINTER :: Aentry, Ltmp
    TYPE(ValueList_t), POINTER :: BC
    TYPE(ListMatrix_t), POINTER :: Alist(:)
    INTEGER :: i,j,k,l,n,Start,nCount,fixedge
    LOGICAL, ALLOCATABLE :: Done(:)
    INTEGER, ALLOCATABLE :: NodeList(:)
    TYPE(Element_t), POINTER :: Edge, Boundary, Element
!------------------------------------------------------------------------------

    IF( ALLOCATED( TreeEdges ) ) THEN
      CALL Info('WhitneyAVSolver','Boundary Gauge tree already created',Level=15)
      RETURN
    END IF

    ALLOCATE(TreeEdges(Mesh % NumberOfEdges))
    TreeEdges = .FALSE.

    n = Mesh % NumberOfNodes
    ALLOCATE(Done(n)); Done=.FALSE.

    !
    ! list the candidate nodes:
    ! -------------------------
    DO i=1,FluxCount
      j = FluxMap(i)
      Edge => Mesh % Edges(j)
      Done(Edge % NodeIndexes)=.TRUE.
    END DO

    ALLOCATE(NodeList(COUNT(Done)))
    nCount = 0
    DO i=1,n
      IF ( Done(i) ) THEN
        nCount = nCount+1
        NodeList(nCount)=i
      END IF
    END DO

    Done=.FALSE.
    DO i=1,FluxCount
      IF ( TreeEdges(FluxMap(i)) ) THEN
        Edge => Mesh % Edges(FluxMap(i))
        Done(Edge % NodeIndexes)=.TRUE.
      END IF
    END DO

    ! 
    ! Skip Dirichlet BCs in terms of A:
    ! ---------------------------------
    DO i=1,Mesh % NumberOfBoundaryElements
      Boundary => GetBoundaryElement(i)
      SELECT CASE(GetElementFamily())
      CASE(1)
        CYCLE
      CASE(2)
        k = GetBoundaryEdgeIndex(Boundary,1); Element => Mesh % Edges(k)
      CASE(3,4)
        k = GetBoundaryFaceIndex(Boundary)  ; Element => Mesh % Faces(k)
      END SELECT
      IF (.NOT. ActiveBoundaryElement(Element)) CYCLE
      BC => GetBC()
      IF (.NOT.ASSOCIATED(BC)) CYCLE
      IF (.NOT.ListCheckPresent( BC, &
           TRIM(Solver % Variable % Name)//' {e}')) CYCLE
 
      j=1; k=GetBoundaryEdgeIndex(Boundary,j)
      DO WHILE(k>0)
        Edge => Mesh % Edges(k)
        TreeEdges(k) = .TRUE.
        Done(Edge % NodeIndexes) = .TRUE.
        j=j+1; k=GetBoundaryEdgeIndex(Boundary,j)
      END DO
    END DO

    ! node -> edge list
    ! -----------------
    Alist => NULL()
    DO i=1,FluxCount
      j = FluxMap(i)
      IF ( Solver % Variable % Perm(j+n)<=0 ) CYCLE

      Edge => Mesh % Edges(j)
      DO k=1,Edge % TYPE % NumberOfNodes
        l=Edge % NodeIndexes(k)
        Aentry=>List_GetMatrixIndex(Alist,l,j)
      END DO
    END DO

    ALLOCATE(BasicCycles(FluxCount))
    BasicCycles(:) % Degree = 0
    DO i=1,FluxCount
       BasicCycles(i) % Head => NULL()
    END DO

    ! generate the tree for all (perhaps disconnected) parts:
    ! -------------------------------------------------------
    DO WHILE(.NOT.ALL(Done(NodeList)))
      DO i=1,nCount
        Start = NodeList(i)
        IF ( .NOT. Done(Start) ) EXIT
      END DO
      CALL BreadthFirstSearch(Alist,Done,Start,nCount,NodeList)
    END DO
    DEALLOCATE(Done,NodeList)
    CALL List_FreeMatrix(SIZE(Alist),Alist)

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE BreadthFirstSearch(Alist,done,start,nCount,NodeList)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: start,nCount,NodeList(:)
    LOGICAL :: Done(:)
    TYPE(ListMatrix_t) :: Alist(:)
!------------------------------------------------------------------------------
    TYPE(ListMatrixEntry_t), POINTER :: Aentry, Ltmp, Btmp
    INTEGER :: i,j,k,l,n,m,ll,IF,bcycle
    LOGICAL :: FirstTime= .TRUE.
    TYPE(Element_t), POINTER :: Edge,Edge1,Boundary
    LOGICAL, ALLOCATABLE :: DoneL(:)
    INTEGER, ALLOCATABLE :: Fifo(:), Previous(:), FiFo1(:)
!------------------------------------------------------------------------------

   ALLOCATE(DoneL(Mesh % NumberOfEdges)); DoneL=.FALSE.
   ALLOCATE(Fifo(FluxCount),FiFo1(FluxCount))
   ALLOCATE(Previous(Mesh % NumberOfNodes)); Previous=0;

   bcycle = 0
   DO bcycle=0,FluxCount
     IF(.NOT.ASSOCIATED(BasicCycles(bcycle+1) % Head)) EXIT
   END DO

   IF = 0; m=0

   IF( FirstTime ) THEN
     DO i=1,nCount
       j = NodeList(i)
       IF ( Done(j) ) THEN
         m=m+1; fifo1(m)=j
         IF=IF+1; fifo(IF)=j
       END IF
     END DO
     FirstTime = .FALSE.
   END IF

   IF ( IF==0 ) THEN
     Done(Start)=.TRUE.
     m=m+1; fifo1(m)=start
     IF=1; fifo(IF)=start;
   END IF

   IF ( IF>0 ) THEN
     DO WHILE(m>0)
       j = Fifo1(m); m=m-1

       Aentry => Alist(j) % Head
       DO WHILE(ASSOCIATED(Aentry))
         k = Aentry % Index
         Aentry => Aentry % Next

         Edge => Mesh % Edges(k)
         IF (.NOT. TreeEdges(k) .OR. DoneL(k) ) CYCLE
         DoneL(k)=.TRUE.

         l = Edge % NodeIndexes(1)
         IF (l==j) l=Edge % NodeIndexes(2)

         IF=IF+1; Fifo(IF)=l
         m=m+1; Fifo1(m)=l
         Previous(l)=j
       END DO
     END DO
     Start = l
   END IF


   DO WHILE(IF>0)
     j = Fifo(IF); IF=IF-1

     Aentry => Alist(j) % Head
     DO WHILE(ASSOCIATED(Aentry))
       k = Aentry % Index
       Aentry => Aentry % Next

       Edge => Mesh % Edges(k)
       IF ( DoneL(k) ) CYCLE
       DoneL(k)=.TRUE.

       l = Edge % NodeIndexes(1)
       IF (l==j) l=Edge % NodeIndexes(2)

       IF ( Done(l) ) THEN
         ! Generate fundamental cycle
         bcycle = bcycle+1
         CALL AddToCycle(bcycle,k)

         m = j
         DO WHILE(m/=Previous(l))
           Ltmp => Alist(m) % Head
           DO WHILE(ASSOCIATED(Ltmp))
             Edge1 => Mesh % Edges(Ltmp % Index)
             IF ( ANY(Edge1 % NodeIndexes(1:2)==Previous(m)) ) THEN
               CALL AddToCycle(bcycle,Ltmp % Index); EXIT
             END IF
             Ltmp=>Ltmp % Next
           END DO
           IF (ANY(Edge1 % NodeIndexes(1:2) == l) ) EXIT
           m = Previous(m)
           IF(m==0) EXIT
         END DO

         IF (ALL(Edge1 % NodeIndexes(1:2) /= l) ) THEN
           ltmp => Alist(l) % Head
           DO WHILE(ASSOCIATED(ltmp))
             edge1 => Mesh % Edges(Ltmp % Index)
             IF ( ANY(Edge1 % NodeIndexes(1:2)==Previous(l)) ) THEN
               CALL AddToCycle(bcycle,Ltmp % Index); EXIT
             END IF
             ltmp=>ltmp % Next
           END DO
         END IF
       ELSE
         IF (.NOT.TreeEdges(k)) CALL SetDOFToValue(Solver,k,0._dp)
         IF=IF+1; Fifo(IF)=l
         Previous(l)=j
         Done(l)=.TRUE.
         TreeEdges(k) = .TRUE.
       END IF
     END DO
   END DO
   DEALLOCATE(Fifo, Fifo1, DoneL)
!------------------------------------------------------------------------------
  END SUBROUTINE BreadthFirstSearch
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AddToCycle(bcycle,index)
    IMPLICIT NONE
    INTEGER :: bcycle,index
!------------------------------------------------------------------------------
    TYPE(ListMatrixEntry_t), POINTER :: Btmp

    ALLOCATE(Btmp); Btmp % Next => BasicCycles(bcycle) % Head;
    Btmp % Index = index; BasicCycles(bcycle) % Head => Btmp
    BasicCycles(bcycle) % Degree=BasicCycles(bcycle) % Degree+1
!------------------------------------------------------------------------------
  END SUBROUTINE AddToCycle
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE GaugeTreeFluxBC
!------------------------------------------------------------------------------

END MODULE MagnetoDynamicsUtils
