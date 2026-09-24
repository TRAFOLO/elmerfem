!
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
!
!> Block preconditioning for field systems coupled to circuits.
!>
!> The circuit unknowns are few, so their block is assembled densely on every
!> rank, summed with one allreduce and LU-factored with LAPACK; each rank then
!> solves the same small system (no sparse direct solver, no gather to one rank).
!> The field block is preconditioned with ILU (default), nothing, or an
!> auxiliary space cycle: ILU smoothing around a PRESB correction with the real
!> edge operator K + omega*sigma*M, solved by Hypre AMS in the helper solver
!> HarmonicEdgePrecSolver (MagnetoDynamics).
!
!> \ingroup ElmerLib
!> \{

#include "../config.h"

MODULE CircuitAuxPrec

  USE DefUtils
  USE SParIterSolve, ONLY : SParCMatrixVector, SParMatrixVector

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: CircuitBlockCreate, CircuitBlockPrecComplex, CircuitBlockPrecReal

  INTEGER, PARAMETER :: FIELD_ILU = 0, FIELD_NONE = 1, FIELD_AUX = 2
  INTEGER, PARAMETER :: BLOCK_DIAGONAL = 0, BLOCK_SCHUR = 1

  TYPE CircuitBlock_t
    LOGICAL :: IsComplex = .FALSE.
    LOGICAL :: Parallel = .FALSE.
    LOGICAL :: Factorized = .FALSE.
    INTEGER :: nm = 0
    INTEGER :: nc = 0
    INTEGER :: ng = 0
    INTEGER :: nOwnedF = 0
    INTEGER :: nOwnedC = 0
    INTEGER :: FieldMethod = FIELD_ILU
    INTEGER :: BlockMethod = BLOCK_DIAGONAL
    INTEGER :: PreSmooth = 1
    INTEGER :: PostSmooth = 1
    INTEGER :: Comm = 0
    INTEGER, ALLOCATABLE :: FRow(:)
    INTEGER, ALLOCATABLE :: CPos(:), CGlob(:)
    INTEGER, ALLOCATABLE :: Piv(:)
    COMPLEX(KIND=dp), ALLOCATABLE :: LUc(:,:)
    REAL(KIND=dp), ALLOCATABLE :: LUr(:,:)
    ! Schur complement S = C - A_cf P^-1 A_fc, built at the first application:
    ! Craw = C, Zc(:,k) = P^-1 A_fc(:,CoupCol(k)) for the circuit columns coupled to the field.
    LOGICAL :: SchurReady = .FALSE.
    INTEGER :: nCoup = 0
    INTEGER, ALLOCATABLE :: CoupCol(:)
    COMPLEX(KIND=dp), ALLOCATABLE :: Craw(:,:), Zc(:,:)
    REAL(KIND=dp), POINTER :: DiagScaling(:) => NULL()
    ! Edge auxiliary space (HarmonicEdgePrecSolver, Hypre AMS): owned field index -> edge dof
    LOGICAL :: HasEdgeSlave = .FALSE., EdgeMapReady = .FALSE.
    INTEGER, ALLOCATABLE :: EdgeMap(:)
  END TYPE CircuitBlock_t

  TYPE(CircuitBlock_t), SAVE, TARGET :: CB

CONTAINS

!------------------------------------------------------------------------------
!> Called for every linear solve with "Linear System Preconditioning = circuit".
!> A is the (possibly circuit-augmented) matrix after scaling, before gluing.
!------------------------------------------------------------------------------
  SUBROUTINE CircuitBlockCreate(A, Solver)
!------------------------------------------------------------------------------
    TYPE(Matrix_t), TARGET :: A
    TYPE(Solver_t) :: Solver
!------------------------------------------------------------------------------
    INTEGER :: i, j, k, gi, gj, pos, linfo, ierr, step, nOwned
    LOGICAL :: Found
    CHARACTER(:), ALLOCATABLE :: str
    CHARACTER(*), PARAMETER :: Caller = 'CircuitBlockCreate'
!------------------------------------------------------------------------------
    CB % Parallel = ParEnv % PEs > 1
    CB % IsComplex = A % COMPLEX
    CB % nm = A % NumberOfRows - A % ExtraDOFs
    CB % nc = MAX(0, A % ParallelDOFs)
    step = 1
    IF (CB % IsComplex) step = 2
    CB % ng = CB % nc / step
    CB % Comm = A % Comm
    CB % Factorized = .FALSE.

    str = ListGetString(Solver % Values, 'Circuit Prec Field Method', Found)
    IF (.NOT. Found) str = 'ilu'
    SELECT CASE (str)
    CASE ('ilu')
      CB % FieldMethod = FIELD_ILU
    CASE ('none')
      CB % FieldMethod = FIELD_NONE
    CASE ('auxiliary space')
      CB % FieldMethod = FIELD_AUX
      IF (.NOT. CB % IsComplex) CALL Fatal(Caller, &
          '"auxiliary space" field preconditioning is implemented for complex systems only')
    CASE DEFAULT
      CALL Fatal(Caller, 'Unknown "Circuit Prec Field Method": '//TRIM(str))
    END SELECT

    str = ListGetString(Solver % Values, 'Circuit Prec Block Method', Found)
    IF (.NOT. Found) str = 'diagonal'
    SELECT CASE (str)
    CASE ('diagonal')
      CB % BlockMethod = BLOCK_DIAGONAL
    CASE ('schur')
      CB % BlockMethod = BLOCK_SCHUR
      IF (.NOT. CB % IsComplex) CALL Fatal(Caller, &
          '"schur" circuit block is implemented for complex systems only')
    CASE DEFAULT
      CALL Fatal(Caller, 'Unknown "Circuit Prec Block Method": '//TRIM(str))
    END SELECT
    CB % SchurReady = .FALSE.
    CB % nCoup = 0
    CB % PreSmooth = ListGetInteger(Solver % Values, 'Auxiliary Pre Smoothing Iterations', Found)
    IF (.NOT. Found) CB % PreSmooth = 1
    CB % PostSmooth = ListGetInteger(Solver % Values, 'Auxiliary Post Smoothing Iterations', Found)
    IF (.NOT. Found) CB % PostSmooth = 1

    ! Positions of the owned unknowns in the Krylov vectors: owned rows of A in order.
    IF (ALLOCATED(CB % FRow)) DEALLOCATE(CB % FRow)
    IF (ALLOCATED(CB % CPos)) DEALLOCATE(CB % CPos, CB % CGlob)
    ALLOCATE(CB % FRow(CB % nm / step), CB % CPos(MAX(1,CB % ng)), CB % CGlob(MAX(1,CB % ng)))
    pos = 0
    CB % nOwnedF = 0
    CB % nOwnedC = 0
    DO i = 1, A % NumberOfRows, step
      IF (CB % Parallel) THEN
        IF (A % ParallelInfo % NeighbourList(i) % Neighbours(1) /= ParEnv % MyPE) CYCLE
      END IF
      pos = pos + 1
      IF (i <= CB % nm) THEN
        CB % nOwnedF = CB % nOwnedF + 1
        CB % FRow(CB % nOwnedF) = i
      ELSE IF (i <= CB % nm + CB % nc) THEN
        CB % nOwnedC = CB % nOwnedC + 1
        CB % CPos(CB % nOwnedC) = pos
        CB % CGlob(CB % nOwnedC) = (i - CB % nm - 1) / step + 1
      END IF
    END DO
    nOwned = pos

    CB % DiagScaling => NULL()
    IF (A % ScalingMethod == 1 .AND. ASSOCIATED(A % DiagScaling)) CB % DiagScaling => A % DiagScaling

    IF (CB % ng > 0) THEN
      IF (ALLOCATED(CB % Piv)) DEALLOCATE(CB % Piv)
      ALLOCATE(CB % Piv(CB % ng))
      IF (CB % IsComplex) THEN
        IF (ALLOCATED(CB % LUc)) DEALLOCATE(CB % LUc)
        ALLOCATE(CB % LUc(CB % ng, CB % ng))
        CB % LUc = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      ELSE
        IF (ALLOCATED(CB % LUr)) DEALLOCATE(CB % LUr)
        ALLOCATE(CB % LUr(CB % ng, CB % ng))
        CB % LUr = 0.0_dp
      END IF

      ! Each rank holds partial (element-wise) values of the circuit rows: sum them.
      DO i = CB % nm + 1, CB % nm + CB % nc, step
        gi = (i - CB % nm - 1) / step + 1
        DO j = A % Rows(i), A % Rows(i+1) - 1, step
          k = A % Cols(j)
          IF (k <= CB % nm .OR. k > CB % nm + CB % nc) CYCLE
          gj = (k - CB % nm - 1) / step + 1
          IF (CB % IsComplex) THEN
            CB % LUc(gi,gj) = CB % LUc(gi,gj) + CMPLX(A % Values(j), -A % Values(j+1), KIND=dp)
          ELSE
            CB % LUr(gi,gj) = CB % LUr(gi,gj) + A % Values(j)
          END IF
        END DO
      END DO

      IF (CB % Parallel) THEN
        IF (CB % IsComplex) THEN
          CALL ParSumC(CB % LUc, CB % ng**2)
        ELSE
          CALL ParSumR(CB % LUr, CB % ng**2)
        END IF
      END IF

      IF (CB % IsComplex) THEN
        IF (CB % BlockMethod == BLOCK_SCHUR) THEN
          IF (ALLOCATED(CB % Craw)) DEALLOCATE(CB % Craw)
          ALLOCATE(CB % Craw(CB % ng, CB % ng))
          CB % Craw = CB % LUc
        END IF
        CALL ZGETRF(CB % ng, CB % ng, CB % LUc, CB % ng, CB % Piv, linfo)
      ELSE
        CALL DGETRF(CB % ng, CB % ng, CB % LUr, CB % ng, CB % Piv, linfo)
      END IF
      CB % Factorized = (linfo == 0)
      IF (.NOT. CB % Factorized) THEN
        CALL Warn(Caller, 'Circuit block is singular (LAPACK info '//I2S(linfo)// &
            '), circuit unknowns left unpreconditioned')
      END IF
    END IF

    IF (CB % FieldMethod == FIELD_AUX) CALL AuxSetup(Solver)

    CALL Info(Caller, 'Circuit unknowns: '//I2S(CB % ng)//', owned here: '//I2S(CB % nOwnedC)// &
        ', owned field unknowns: '//I2S(CB % nOwnedF)//' of '//I2S(nOwned), Level=8)
!------------------------------------------------------------------------------
  END SUBROUTINE CircuitBlockCreate
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AuxSetup(Solver)
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
!------------------------------------------------------------------------------
    LOGICAL :: Found
    CHARACTER(*), PARAMETER :: Caller = 'CircuitBlockCreate'
!------------------------------------------------------------------------------
    CB % HasEdgeSlave = ListGetInteger(Solver % Values, 'Auxiliary Edge Solver', Found) > 0
    CB % EdgeMapReady = .FALSE.
    IF (.NOT. CB % HasEdgeSlave) CALL Fatal(Caller, &
        '"auxiliary space" needs an "Auxiliary Edge Solver" (HarmonicEdgePrecSolver)')
!------------------------------------------------------------------------------
  END SUBROUTINE AuxSetup
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> u = M^-1 v for complex systems (HUTI convention, owned unknowns only).
!------------------------------------------------------------------------------
  SUBROUTINE CircuitBlockPrecComplex(u, v, ipar)
!------------------------------------------------------------------------------
    INTEGER :: ipar(*)
    COMPLEX(KIND=dp) :: u(*), v(*)
!------------------------------------------------------------------------------
    INTEGER :: ndim, nf
!------------------------------------------------------------------------------
    ndim = ipar(3)
    nf = CB % nOwnedF

    IF (CB % BlockMethod == BLOCK_SCHUR .AND. CB % ng > 0 .AND. .NOT. CB % SchurReady) THEN
      CALL SchurSetup(ipar)
    END IF

    CALL FieldApply(u, v, ipar)
    CALL CircuitSolveComplex(u, v, ipar)
!------------------------------------------------------------------------------
  END SUBROUTINE CircuitBlockPrecComplex
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Field part of the preconditioner: u(1:nOwnedF) from v(1:nOwnedF).
!------------------------------------------------------------------------------
  SUBROUTINE FieldApply(u, v, ipar)
!------------------------------------------------------------------------------
    INTEGER :: ipar(*)
    COMPLEX(KIND=dp) :: u(*), v(*)
!------------------------------------------------------------------------------
    INTEGER :: ndim, nf
!------------------------------------------------------------------------------
    ndim = ipar(3)
    nf = CB % nOwnedF
    SELECT CASE (CB % FieldMethod)
    CASE (FIELD_ILU)
      CALL CRS_ComplexLUPrecondition(u, v, ipar)
    CASE (FIELD_NONE)
      u(1:ndim) = v(1:ndim)
    CASE (FIELD_AUX)
      u(nf+1:ndim) = v(nf+1:ndim)
      CALL AuxFieldCycle(u, v, ipar)
    END SELECT
!------------------------------------------------------------------------------
  END SUBROUTINE FieldApply
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Builds S = C - A_cf P^-1 A_fc for the circuit columns that couple to the
!> field (one field preconditioner application per coupled column) and
!> factorizes it redundantly on all ranks.
!------------------------------------------------------------------------------
  SUBROUTINE SchurSetup(ipar)
!------------------------------------------------------------------------------
    INTEGER :: ipar(*)
!------------------------------------------------------------------------------
    INTEGER :: ndim, nf, ng, j, k, m, linfo, ierr
    LOGICAL, ALLOCATABLE :: Coupled(:)
    REAL(KIND=dp) :: fn
    COMPLEX(KIND=dp), ALLOCATABLE :: t(:), w(:), z(:), colw(:), S(:,:)
    CHARACTER(*), PARAMETER :: Caller = 'CircuitBlockSchur'
!------------------------------------------------------------------------------
    ndim = ipar(3)
    nf = CB % nOwnedF
    ng = CB % ng
    ALLOCATE(t(ndim), w(ndim), z(ndim), colw(ng), S(ng,ng), Coupled(ng))

    ! Field part of each circuit column: w = A e_j restricted to field rows.
    DO j = 1, ng
      CALL UnitColumn(j)
      CALL KrylovMatVecComplex(t, w, ipar)
      fn = SUM(ABS(w(1:nf))**2)
      IF (CB % Parallel) CALL ParSumRS(fn)
      Coupled(j) = fn > 0.0_dp
    END DO
    CB % nCoup = COUNT(Coupled)
    IF (ALLOCATED(CB % CoupCol)) DEALLOCATE(CB % CoupCol)
    IF (ALLOCATED(CB % Zc)) DEALLOCATE(CB % Zc)
    ALLOCATE(CB % CoupCol(MAX(1,CB % nCoup)), CB % Zc(MAX(1,nf), MAX(1,CB % nCoup)))

    S = CB % Craw
    m = 0
    DO j = 1, ng
      IF (.NOT. Coupled(j)) CYCLE
      m = m + 1
      CB % CoupCol(m) = j
      ! z = P^-1 A_fc(:,j)
      CALL UnitColumn(j)
      CALL KrylovMatVecComplex(t, w, ipar)
      t = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      t(1:nf) = w(1:nf)
      z = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      CALL FieldApply(z, t, ipar)
      CB % Zc(1:nf,m) = z(1:nf)
      ! S(:,j) = C(:,j) - A_cf z
      t = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      t(1:nf) = z(1:nf)
      CALL KrylovMatVecComplex(t, w, ipar)
      colw = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO k = 1, CB % nOwnedC
        colw(CB % CGlob(k)) = w(CB % CPos(k))
      END DO
      IF (CB % Parallel) CALL ParSumC(colw, ng)
      S(:,j) = S(:,j) - colw
    END DO

    CB % LUc = S
    CALL ZGETRF(ng, ng, CB % LUc, ng, CB % Piv, linfo)
    CB % Factorized = (linfo == 0)
    IF (.NOT. CB % Factorized) CALL Warn(Caller, 'Schur complement is singular (LAPACK info '//I2S(linfo)//')')
    CB % SchurReady = .TRUE.
    CALL Info(Caller, 'Circuit Schur complement: '//I2S(CB % nCoup)//' of '//I2S(ng)// &
        ' circuit unknowns couple to the field', Level=6)

  CONTAINS

    ! t = e_j in the Krylov vector (set by the owner of circuit unknown j)
    SUBROUTINE UnitColumn(j)
      INTEGER :: j, kk
      t = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO kk = 1, CB % nOwnedC
        IF (CB % CGlob(kk) == j) t(CB % CPos(kk)) = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
      END DO
    END SUBROUTINE UnitColumn
!------------------------------------------------------------------------------
  END SUBROUTINE SchurSetup
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> u = M^-1 v for real systems: ILU (or identity) on the field, dense circuit block.
!------------------------------------------------------------------------------
  SUBROUTINE CircuitBlockPrecReal(u, v, ipar)
!------------------------------------------------------------------------------
    INTEGER :: ipar(*)
    REAL(KIND=dp) :: u(*), v(*)
!------------------------------------------------------------------------------
    INTEGER :: ndim, k, linfo
    REAL(KIND=dp), ALLOCATABLE :: vg(:)
!------------------------------------------------------------------------------
    ndim = ipar(3)
    IF (CB % FieldMethod == FIELD_NONE) THEN
      u(1:ndim) = v(1:ndim)
    ELSE
      CALL CRS_LUPrecondition(u, v, ipar)
    END IF
    IF (CB % ng == 0) RETURN

    ALLOCATE(vg(CB % ng))
    vg = 0.0_dp
    DO k = 1, CB % nOwnedC
      vg(CB % CGlob(k)) = v(CB % CPos(k))
    END DO
    IF (CB % Parallel) THEN
      CALL ParSumR(vg, CB % ng)
    END IF
    IF (CB % Factorized) THEN
      CALL DGETRS('N', CB % ng, 1, CB % LUr, CB % ng, CB % Piv, vg, CB % ng, linfo)
    END IF
    DO k = 1, CB % nOwnedC
      u(CB % CPos(k)) = vg(CB % CGlob(k))
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE CircuitBlockPrecReal
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Overwrites the owned circuit entries of u with the dense circuit block solve.
!------------------------------------------------------------------------------
  SUBROUTINE CircuitSolveComplex(u, v, ipar)
!------------------------------------------------------------------------------
    INTEGER :: ipar(*)
    COMPLEX(KIND=dp) :: u(*), v(*)
!------------------------------------------------------------------------------
    INTEGER :: ndim, k, linfo, ierr
    COMPLEX(KIND=dp), ALLOCATABLE :: vg(:), t(:), w(:)
!------------------------------------------------------------------------------
    IF (CB % ng == 0) RETURN
    ndim = ipar(3)

    ALLOCATE(vg(CB % ng))
    vg = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
    IF (CB % BlockMethod == BLOCK_SCHUR) THEN
      ALLOCATE(t(ndim), w(ndim))
      t = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      t(1:CB % nOwnedF) = u(1:CB % nOwnedF)
      CALL KrylovMatVecComplex(t, w, ipar)
      DO k = 1, CB % nOwnedC
        vg(CB % CGlob(k)) = v(CB % CPos(k)) - w(CB % CPos(k))
      END DO
    ELSE
      DO k = 1, CB % nOwnedC
        vg(CB % CGlob(k)) = v(CB % CPos(k))
      END DO
    END IF
    IF (CB % Parallel) THEN
      CALL ParSumC(vg, CB % ng)
    END IF
    IF (CB % Factorized) THEN
      CALL ZGETRS('N', CB % ng, 1, CB % LUc, CB % ng, CB % Piv, vg, CB % ng, linfo)
    END IF
    DO k = 1, CB % nOwnedC
      u(CB % CPos(k)) = vg(CB % CGlob(k))
    END DO
    ! Back substitution of the block factorization: u_f = y_f - P^-1 A_fc u_c.
    IF (CB % BlockMethod == BLOCK_SCHUR) THEN
      DO k = 1, CB % nCoup
        u(1:CB % nOwnedF) = u(1:CB % nOwnedF) - CB % Zc(1:CB % nOwnedF,k) * vg(CB % CoupCol(k))
      END DO
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE CircuitSolveComplex
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Owned Krylov field index -> dof of the real edge auxiliary solver (0 for
!> nodal unknowns).
!------------------------------------------------------------------------------
  SUBROUTINE EdgeMapSetup(Solver, ESolver)
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver, ESolver
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, ALLOCATABLE :: KIdx(:)
    INTEGER :: e, gi, km, ks, i, nf, cnt
!------------------------------------------------------------------------------
    Mesh => Solver % Mesh
    nf = CB % nOwnedF
    ALLOCATE(KIdx(CB % nm / 2))
    KIdx = 0
    DO i = 1, nf
      KIdx((CB % FRow(i) + 1) / 2) = i
    END DO
    IF (ALLOCATED(CB % EdgeMap)) DEALLOCATE(CB % EdgeMap)
    ALLOCATE(CB % EdgeMap(MAX(1,nf)))
    CB % EdgeMap = 0
    cnt = 0
    DO e = 1, Mesh % NumberOfEdges
      gi = Mesh % MaxNDOFs * Mesh % NumberOfNodes + Mesh % MaxEdgeDOFs * (e-1) + 1
      km = Solver % Variable % Perm(gi)
      ks = ESolver % Variable % Perm(gi)
      IF (km <= 0 .OR. ks <= 0 .OR. km > CB % nm / 2) CYCLE
      i = KIdx(km)
      IF (i == 0) CYCLE
      CB % EdgeMap(i) = ks
      cnt = cnt + 1
    END DO
    CB % EdgeMapReady = .TRUE.
    CALL Info('CircuitBlockEdge', 'Owned field unknowns mapped to the edge space: '//I2S(cnt), Level=6)
    IF (cnt == 0) THEN
      km = 0; ks = 0
      DO e = 1, Mesh % NumberOfEdges
        gi = Mesh % MaxNDOFs * Mesh % NumberOfNodes + Mesh % MaxEdgeDOFs * (e-1) + 1
        IF (Solver % Variable % Perm(gi) > 0) km = km + 1
        IF (ESolver % Variable % Perm(gi) > 0) ks = ks + 1
      END DO
      CALL Warn('CircuitBlockEdge', 'No field unknown maps to the edge space: MaxNDOFs '// &
          I2S(Mesh % MaxNDOFs)//', MaxEdgeDOFs '//I2S(Mesh % MaxEdgeDOFs)//', edges in master '// &
          I2S(km)//', in edge solver '//I2S(ks)//', perm sizes '//I2S(SIZE(Solver % Variable % Perm))// &
          ' '//I2S(SIZE(ESolver % Variable % Perm)))
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE EdgeMapSetup
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Edge-space correction with the real operator H = K + omega*sigma*M (Hypre AMS
!> slave) for r = f + i g in the unscaled field space, applied as PRESB:
!> h = H^-1 (f + g), y = H^-1 (K h - f), dx = (h - y) + i y. The eigenvalues of
!> the preconditioned system lie in [1/2, 1].
!------------------------------------------------------------------------------
  SUBROUTINE AuxEdgeApply(r, dx, nf, ipar)
!------------------------------------------------------------------------------
    INTEGER :: nf, ipar(*)
    COMPLEX(KIND=dp) :: r(nf), dx(nf)
!------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: Solver, ESolver
    TYPE(Variable_t), POINTER :: ResE
    INTEGER :: i, ks, idx, ndim
    REAL(KIND=dp), ALLOCATABLE :: d(:), h(:)
    COMPLEX(KIND=dp), ALLOCATABLE :: t(:), w(:)
    LOGICAL :: Found
!------------------------------------------------------------------------------
    Solver => CurrentModel % Solver
    idx = ListGetInteger(Solver % Values, 'Auxiliary Edge Solver', Found)
    ESolver => CurrentModel % Solvers(idx)
    IF (.NOT. CB % EdgeMapReady) CALL EdgeMapSetup(Solver, ESolver)
    ResE => VariableGet(Solver % Mesh % Variables, 'AuxEdge res', ThisOnly=.TRUE., UnfoundFatal=.TRUE.)

    ALLOCATE(d(nf))
    d = 1.0_dp
    IF (ASSOCIATED(CB % DiagScaling)) THEN
      DO i = 1, nf
        d(i) = CB % DiagScaling(CB % FRow(i))
      END DO
    END IF
    dx = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

    ! h = H^-1 (f + g)
    CALL EdgeSolve()
    ALLOCATE(h(SIZE(ESolver % Variable % Values)))
    h = ESolver % Variable % Values
    ! K h = Re(A h) for real h, with the iterated (scaled) matrix
    ndim = ipar(3)
    ALLOCATE(t(ndim), w(ndim))
    t = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
    DO i = 1, nf
      ks = CB % EdgeMap(i)
      IF (ks > 0) t(i) = CMPLX(h(ks) / d(i), 0.0_dp, KIND=dp)
    END DO
    CALL KrylovMatVecComplex(t, w, ipar)
    ! y = H^-1 (K h - f)
    ResE % Values = 0.0_dp
    DO i = 1, nf
      ks = CB % EdgeMap(i)
      IF (ks > 0) ResE % Values(ks) = (REAL(w(i), KIND=dp) - REAL(r(i), KIND=dp)) / d(i)
    END DO
    CALL DefaultSlaveSolvers(Solver, 'Auxiliary Edge Solver')
    CurrentModel % Solver => Solver
    DO i = 1, nf
      ks = CB % EdgeMap(i)
      IF (ks == 0) CYCLE
      dx(i) = CMPLX(h(ks) - ESolver % Variable % Values(ks), ESolver % Variable % Values(ks), KIND=dp) / d(i)
    END DO

  CONTAINS

    ! Unscaled right-hand side Re r + Im r, then solve.
    SUBROUTINE EdgeSolve()
      INTEGER :: ii, kk
      ResE % Values = 0.0_dp
      DO ii = 1, nf
        kk = CB % EdgeMap(ii)
        IF (kk == 0) CYCLE
        ResE % Values(kk) = (REAL(r(ii), KIND=dp) + AIMAG(r(ii))) / d(ii)
      END DO
      CALL DefaultSlaveSolvers(Solver, 'Auxiliary Edge Solver')
      CurrentModel % Solver => Solver
    END SUBROUTINE EdgeSolve
!------------------------------------------------------------------------------
  END SUBROUTINE AuxEdgeApply
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Field part: smoothing, PRESB edge correction, smoothing, each step with a
!> minimal-residual step length. Only the owned field entries u(1:nOwnedF) are written.
!------------------------------------------------------------------------------
  SUBROUTINE AuxFieldCycle(u, v, ipar)
!------------------------------------------------------------------------------
    INTEGER :: ipar(*)
    COMPLEX(KIND=dp) :: u(*), v(*)
!------------------------------------------------------------------------------
    INTEGER :: nf, ndim, k
    COMPLEX(KIND=dp), ALLOCATABLE :: x(:), r(:), dx(:), t(:), w(:)
!------------------------------------------------------------------------------
    nf = CB % nOwnedF
    ndim = ipar(3)
    ALLOCATE(x(nf), r(nf), dx(nf), t(ndim), w(ndim))
    x = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
    r = v(1:nf)

    DO k = 1, CB % PreSmooth
      CALL SmoothStep()
      CALL ApplyStep()
    END DO

    CALL AuxEdgeApply(r, dx, nf, ipar)
    CALL ApplyStep()

    DO k = 1, CB % PostSmooth
      CALL SmoothStep()
      CALL ApplyStep()
    END DO

    u(1:nf) = x

  CONTAINS

    ! x += alpha dx and r -= alpha (A [dx;0])_f with alpha minimizing the field
    ! residual (robust to badly scaled corrections).
    SUBROUTINE ApplyStep()
      COMPLEX(KIND=dp) :: alpha, num
      REAL(KIND=dp) :: den
      t = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      t(1:nf) = dx
      CALL KrylovMatVecComplex(t, w, ipar)
      num = SUM(CONJG(w(1:nf)) * r)
      den = SUM(ABS(w(1:nf))**2)
      IF (CB % Parallel) THEN
        CALL ParSumCS(num)
        CALL ParSumRS(den)
      END IF
      alpha = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      IF (den > TINY(den)) alpha = num / den
      x = x + alpha * dx
      r = r - alpha * w(1:nf)
    END SUBROUTINE ApplyStep

    ! dx = S r with S = ILU of the iterated matrix if built, else Jacobi.
    SUBROUTINE SmoothStep()
      INTEGER :: i
      IF (ASSOCIATED(GlobalMatrix % CILUValues)) THEN
        t = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
        t(1:nf) = r
        CALL CRS_ComplexLUPrecondition(w, t, ipar)
        dx = w(1:nf)
      ELSE
        DO i = 1, nf
          dx(i) = r(i) / DiagC(i)
        END DO
      END IF
    END SUBROUTINE SmoothStep

    ! Complex diagonal entry of owned field row i of the iterated matrix
    FUNCTION DiagC(i) RESULT(d)
      INTEGER :: i, j
      COMPLEX(KIND=dp) :: d
      j = GlobalMatrix % Diag(2*i-1)
      d = CMPLX(GlobalMatrix % Values(j), -GlobalMatrix % Values(j+1), KIND=dp)
      IF (ABS(d) <= TINY(1.0_dp)) d = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
    END FUNCTION DiagC

!------------------------------------------------------------------------------
  END SUBROUTINE AuxFieldCycle
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> w = A t with the matrix of the running Krylov solve (owned vectors).
!------------------------------------------------------------------------------
  SUBROUTINE KrylovMatVecComplex(t, w, ipar)
!------------------------------------------------------------------------------
    INTEGER :: ipar(*)
    COMPLEX(KIND=dp) :: t(*), w(*)
!------------------------------------------------------------------------------
    IF (CB % Parallel) THEN
      CALL SParCMatrixVector(t, w, ipar)
    ELSE
      CALL CRS_ComplexMatrixVectorProd(t, w, ipar)
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE KrylovMatVecComplex
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Sums over the ranks of the solver communicator. MPI_IN_PLACE is avoided on
!> purpose: with MS-MPI's Fortran header it is not recognized (garbage results).
!------------------------------------------------------------------------------
  SUBROUTINE ParSumC(x, n)
    INTEGER :: n
    COMPLEX(KIND=dp) :: x(n)
    COMPLEX(KIND=dp), ALLOCATABLE :: tmp(:)
    INTEGER :: ierr
    ALLOCATE(tmp(n))
    tmp = x
    CALL MPI_ALLREDUCE(tmp, x, n, MPI_DOUBLE_COMPLEX, MPI_SUM, CB % Comm, ierr)
  END SUBROUTINE ParSumC

  SUBROUTINE ParSumR(x, n)
    INTEGER :: n
    REAL(KIND=dp) :: x(n)
    REAL(KIND=dp), ALLOCATABLE :: tmp(:)
    INTEGER :: ierr
    ALLOCATE(tmp(n))
    tmp = x
    CALL MPI_ALLREDUCE(tmp, x, n, MPI_DOUBLE_PRECISION, MPI_SUM, CB % Comm, ierr)
  END SUBROUTINE ParSumR

  SUBROUTINE ParSumCS(x)
    COMPLEX(KIND=dp) :: x, tmp
    INTEGER :: ierr
    tmp = x
    CALL MPI_ALLREDUCE(tmp, x, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, CB % Comm, ierr)
  END SUBROUTINE ParSumCS

  SUBROUTINE ParSumRS(x)
    REAL(KIND=dp) :: x, tmp
    INTEGER :: ierr
    tmp = x
    CALL MPI_ALLREDUCE(tmp, x, 1, MPI_DOUBLE_PRECISION, MPI_SUM, CB % Comm, ierr)
  END SUBROUTINE ParSumRS
!------------------------------------------------------------------------------

END MODULE CircuitAuxPrec

!> \}
