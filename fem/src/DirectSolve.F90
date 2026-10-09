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
!
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

!> \ingroup ElmerLib 
!> \{

!------------------------------------------------------------------------------
!>  Module containing the direct solvers for linear systems given in CRS format.
!> Included are Lapack band matrix solver, multifrontal Umfpack, MUMPS, SuperLU, 
!> and Pardiso. Note that many of these are linked in with ElmerSolver only 
!> if they are made available at the compilation time. 
!------------------------------------------------------------------------------

#include "../config.h"

MODULE DirectSolve

   USE CRSMatrix
   USE BandMatrix
   USE SParIterSolve

   IMPLICIT NONE

CONTAINS


!------------------------------------------------------------------------------
!> Solver the complex linear system using direct band matrix solver from of Lapack.
!------------------------------------------------------------------------------
   SUBROUTINE ComplexBandSolver( A,x,b, Free_fact )
!------------------------------------------------------------------------------

     LOGICAL, OPTIONAL :: Free_Fact
     TYPE(Matrix_t) :: A
     REAL(KIND=dp) :: x(*),b(*)
!------------------------------------------------------------------------------

   
     INTEGER :: i,j,k,istat,Subband,N
     COMPLEX(KIND=dp), ALLOCATABLE :: BA(:,:)

     REAL(KIND=dp), POINTER CONTIG :: Values(:)
     INTEGER, POINTER CONTIG :: Rows(:), Cols(:), Diag(:)

     SAVE BA
!------------------------------------------------------------------------------

     IF ( PRESENT(Free_Fact) ) THEN
       IF ( Free_Fact ) THEN
         IF ( ALLOCATED(BA) ) DEALLOCATE(BA)
         RETURN
       END IF
     END IF

     Rows => A % Rows
     Cols => A % Cols
     Diag => A % Diag
     Values => A % Values

     n = A % NumberOfRows
     x(1:n) = b(1:n)
     n = n / 2

     IF ( A % Format == MATRIX_CRS .AND. .NOT. A % Symmetric ) THEN
       Subband = 0
       DO i=1,N
         DO j=Rows(2*i-1),Rows(2*i)-1,2
           Subband = MAX(Subband,ABS((Cols(j)+1)/2-i))
         END DO
       END DO

       IF ( .NOT.ALLOCATED( BA ) ) THEN

         ALLOCATE( BA(3*SubBand+1,N),stat=istat )

         IF ( istat /= 0 ) THEN
           CALL Fatal( 'ComplexBandSolver', 'Memory allocation error.' )
         END IF

       ELSE IF ( SIZE(BA,1) /= 3*Subband+1 .OR. SIZE(BA,2) /= N ) THEN

         DEALLOCATE( BA )
         ALLOCATE( BA(3*SubBand+1,N),stat=istat )

         IF ( istat /= 0 ) THEN
           CALL Fatal( 'ComplexBandSolver', 'Memory allocation error.' )
         END IF

       END IF

       BA = 0.0D0
       DO i=1,N
         DO j=Rows(2*i-1),Rows(2*i)-1,2
           k = i - (Cols(j)+1)/2 + 2*Subband + 1
           BA(k,(Cols(j)+1)/2) = CMPLX(Values(j), -Values(j+1), KIND=dp )
         END DO
       END DO

       CALL SolveComplexBandLapack( N,1,BA,x,Subband,3*Subband+1 )

     ELSE IF ( A % Format == MATRIX_CRS ) THEN

       Subband = 0
       DO i=1,N
         DO j=Rows(2*i-1),Diag(2*i-1)
           Subband = MAX(Subband,ABS((Cols(j)+1)/2-i))
         END DO
       END DO

       IF ( .NOT.ALLOCATED( BA ) ) THEN

         ALLOCATE( BA(SubBand+1,N),stat=istat )

         IF ( istat /= 0 ) THEN
           CALL Fatal( 'ComplexBandSolver', 'Memory allocation error.' )
         END IF

       ELSE IF ( SIZE(BA,1) /= Subband+1 .OR. SIZE(BA,2) /= N ) THEN

         DEALLOCATE( BA )
         ALLOCATE( BA(SubBand+1,N),stat=istat )

         IF ( istat /= 0 ) THEN
           CALL Fatal( 'ComplexBandSolver', 'Direct solver memory allocation error.' )
         END IF

       END IF

       BA = 0.0D0
       DO i=1,N
         DO j=Rows(2*i-1),Diag(2*i-1)
           k = i - (Cols(j)+1)/2 + 1
           BA(k,(Cols(j)+1)/2) = CMPLX(Values(j), -Values(j+1), KIND=dp )
         END DO
       END DO

       CALL SolveComplexSBandLapack( N,1,BA,x,Subband,Subband+1 )

     END IF
!------------------------------------------------------------------------------
  END SUBROUTINE ComplexBandSolver 
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solver the real linear system using direct band matrix solver from of Lapack.
!------------------------------------------------------------------------------
   SUBROUTINE BandSolver( A,x,b,Free_Fact )
!------------------------------------------------------------------------------
     LOGICAL, OPTIONAL :: Free_Fact
     TYPE(Matrix_t) :: A
     REAL(KIND=dp) :: x(*),b(*)
!------------------------------------------------------------------------------

     INTEGER :: i,j,k,istat,Subband,N
     REAL(KIND=dp), ALLOCATABLE :: BA(:,:)

     REAL(KIND=dp), POINTER CONTIG :: Values(:)
     INTEGER, POINTER CONTIG :: Rows(:), Cols(:), Diag(:)

     SAVE BA
!------------------------------------------------------------------------------
     IF ( PRESENT(Free_Fact) ) THEN
       IF ( Free_Fact ) THEN
         IF ( ALLOCATED(BA) ) DEALLOCATE(BA)
         RETURN
       END IF
     END IF

     N = A % NumberOfRows

     x(1:n) = b(1:n)

     Rows => A % Rows
     Cols => A % Cols
     Diag => A % Diag
     Values => A % Values

     IF ( A % Format == MATRIX_CRS ) THEN ! .AND. .NOT. A % Symmetric ) THEN
        Subband = 0
        DO i=1,N
          DO j=Rows(i),Rows(i+1)-1
            Subband = MAX(Subband,ABS(Cols(j)-i))
          END DO
        END DO

        IF ( .NOT.ALLOCATED( BA ) ) THEN

          ALLOCATE( BA(3*SubBand+1,N),stat=istat )

          IF ( istat /= 0 ) THEN
            CALL Fatal( 'BandSolver', 'Memory allocation error.' )
          END IF

        ELSE IF ( SIZE(BA,1) /= 3*Subband+1 .OR. SIZE(BA,2) /= N ) THEN

          DEALLOCATE( BA )
          ALLOCATE( BA(3*SubBand+1,N),stat=istat )

          IF ( istat /= 0 ) THEN
            CALL Fatal( 'BandSolver', 'Memory allocation error.' )
          END IF

       END IF

       BA = 0.0D0
       DO i=1,N
         DO j=Rows(i),Rows(i+1)-1
           k = i - Cols(j) + 2*Subband + 1
           BA(k,Cols(j)) = Values(j)
         END DO
       END DO

       CALL SolveBandLapack( N,1,BA,x,Subband,3*Subband+1 )

     ELSE IF ( A % Format == MATRIX_CRS ) THEN

       Subband = 0
       DO i=1,N
         DO j=Rows(i),Diag(i)
           Subband = MAX(Subband,ABS(Cols(j)-i))
         END DO
       END DO

       IF ( .NOT.ALLOCATED( BA ) ) THEN

         ALLOCATE( BA(SubBand+1,N),stat=istat )

         IF ( istat /= 0 ) THEN
           CALL Fatal( 'BandSolver', 'Memory allocation error.' )
         END IF

       ELSE IF ( SIZE(BA,1) /= Subband+1 .OR. SIZE(BA,2) /= N ) THEN

         DEALLOCATE( BA )
         ALLOCATE( BA(SubBand+1,N),stat=istat )

         IF ( istat /= 0 ) THEN
           CALL Fatal( 'BandSolver', 'Memory allocation error.' )
         END IF

       END IF

       BA = 0.0D0
       DO i=1,N
         DO j=Rows(i),Diag(i)
           k = i - Cols(j) + 1
           BA(k,Cols(j)) = Values(j)
         END DO
       END DO

       CALL SolveSBandLapack( N,1,BA,x,Subband,Subband+1 )

     ELSE IF ( A % Format == MATRIX_BAND ) THEN
       CALL SolveBandLapack( N,1,Values,x,Subband,3*Subband+1 )
     ELSE IF ( A % Format == MATRIX_SBAND ) THEN
       CALL SolveSBandLapack( N,1,Values,x,Subband,Subband+1 )
     END IF

!------------------------------------------------------------------------------
  END SUBROUTINE BandSolver 
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solves a linear system using Umfpack multifrontal direct solver courtesy
!> of University of Florida.
!------------------------------------------------------------------------------
  SUBROUTINE UMFPack_SolveSystem( Solver,A,x,b,Free_Fact )
!------------------------------------------------------------------------------
    LOGICAL, OPTIONAL :: Free_Fact
    TYPE(Matrix_t) :: A
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp), TARGET :: x(*), b(*)

    REAL(KIND=dp), POINTER CONTIG :: Values(:)
    INTEGER, POINTER CONTIG :: Rows(:), Cols(:), Diag(:)

#include "../config.h"
#ifdef HAVE_UMFPACK

#ifdef ARCH_32_BITS
#define CAddrInt c_int32_t
#else
#define CAddrInt c_int64_t
#endif
  ! Standard int version
  INTERFACE
    SUBROUTINE umf4def( control ) &
       BIND(C,name='umf4def')
       USE, INTRINSIC :: ISO_C_BINDING
       REAL(C_DOUBLE) :: control(*)
    END SUBROUTINE umf4def

    SUBROUTINE umf4sym( m,n,rows,cols,values,symbolic,control,iinfo ) &
       BIND(C,name='umf4sym')
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(C_INT) :: m,n,rows(*),cols(*)
       INTEGER(CAddrInt) ::  symbolic
       REAL(C_DOUBLE) :: Values(*), control(*),iinfo(*)
    END SUBROUTINE umf4sym

    SUBROUTINE umf4num( rows,cols,values,symbolic,numeric, control,iinfo ) &
       BIND(C,name='umf4num')
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(C_INT) :: rows(*),cols(*)
       INTEGER(CAddrInt) ::  numeric, symbolic
       REAL(C_DOUBLE) :: Values(*), control(*),iinfo(*)
    END SUBROUTINE umf4num

    SUBROUTINE umf4sol( sys, x, b, numeric, control, iinfo ) &
       BIND(C,name='umf4sol')
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(C_INT) :: sys
       INTEGER(CAddrInt) :: numeric
       REAL(C_DOUBLE) :: x(*), b(*), control(*), iinfo(*)
    END SUBROUTINE umf4sol

    SUBROUTINE umf4fsym(symbolic) &
        BIND(C,name='umf4fsym')
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(CAddrInt) :: symbolic
    END SUBROUTINE umf4fsym

    SUBROUTINE umf4fnum(numeric) &
        BIND(C,name='umf4fnum')
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(CAddrInt) :: numeric
    END SUBROUTINE umf4fnum

  END INTERFACE

  ! Long int version
  INTERFACE
    SUBROUTINE umf4_l_def( control ) &
       BIND(C,name='umf4_l_def')
       USE, INTRINSIC :: ISO_C_BINDING
       REAL(C_DOUBLE) :: control(*)
    END SUBROUTINE umf4_l_def

    SUBROUTINE umf4_l_sym( m,n,rows,cols,values,symbolic,control,iinfo ) &
       BIND(C,name='umf4_l_sym')
       USE, INTRINSIC :: ISO_C_BINDING
       !INTEGER(CAddrInt) ::m,n,rows(*),cols(*) 
       INTEGER(UMFPACK_LONG_FORTRAN_TYPE) :: m,n,rows(*),cols(*) !TODO: m,n of are called with AddrInt kind
       INTEGER(CAddrInt) ::  symbolic
       REAL(C_DOUBLE) :: Values(*), control(*),iinfo(*)
    END SUBROUTINE umf4_l_sym

    SUBROUTINE umf4_l_num( rows,cols,values,symbolic,numeric, control,iinfo ) &
       BIND(C,name='umf4_l_num')
       USE, INTRINSIC :: ISO_C_BINDING
       !INTEGER(CAddrInt) :: rows(*),cols(*)
       INTEGER(UMFPACK_LONG_FORTRAN_TYPE) :: rows(*),cols(*)
       INTEGER(CAddrInt) ::  numeric, symbolic
       REAL(C_DOUBLE) :: Values(*), control(*),iinfo(*)
    END SUBROUTINE umf4_l_num

    SUBROUTINE umf4_l_sol( sys, x, b, numeric, control, iinfo ) &
       BIND(C,name='umf4_l_sol')
       USE, INTRINSIC :: ISO_C_BINDING
       !INTEGER(CAddrInt) :: sys
       INTEGER(UMFPACK_LONG_FORTRAN_TYPE) :: sys
       INTEGER(CAddrInt) :: numeric
       REAL(C_DOUBLE) :: x(*), b(*), control(*), iinfo(*)
    END SUBROUTINE umf4_l_sol

    SUBROUTINE umf4_l_fnum(numeric) &
        BIND(C,name='umf4_l_fnum')
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(CAddrInt) :: numeric
    END SUBROUTINE umf4_l_fnum

    SUBROUTINE umf4_l_fsym(symbolic) &
        BIND(C,name='umf4_l_fsym')
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER(CAddrInt) :: symbolic
    END SUBROUTINE umf4_l_fsym
  END INTERFACE

  INTEGER :: i, n, status, sys
  REAL(KIND=dp) :: iInfo(90), Control(20)
  INTEGER(KIND=AddrInt) :: symbolic, zero=0
  INTEGER(KIND=UMFPACK_LONG_FORTRAN_TYPE) :: ln, lsys
  INTEGER(KIND=UMFPACK_LONG_FORTRAN_TYPE), ALLOCATABLE :: LRows(:), LCols(:)

  SAVE iInfo, Control
 
  LOGICAL :: Factorize, FreeFactorize, stat, BigMode

  IF ( PRESENT(Free_Fact) ) THEN
    IF ( Free_Fact ) THEN
      IF ( A % UMFPack_Numeric/=0 ) THEN
        CALL umf4fnum(A % UMFPack_Numeric)
        A % UMFPack_Numeric = 0
      END IF
      RETURN
    END IF
  END IF

  BigMode = ListGetString( Solver % Values, &
      'Linear System Direct Method' ) == 'big umfpack'

  Factorize = ListGetLogical( Solver % Values, &
     'Linear System Refactorize', stat )
  IF ( .NOT. stat ) Factorize = .TRUE.

  n = A % NumberofRows
  Rows => A % Rows
  Cols => A % Cols
  Diag => A % Diag
  Values => A % Values

  IF ( Factorize .OR. A% UmfPack_Numeric==0 ) THEN
    IF ( A % UMFPack_Numeric /= 0 ) THEN
      IF( BigMode ) THEN
        CALL umf4_l_fnum( A % UMFPack_Numeric )
      ELSE
        CALL umf4fnum( A % UMFPack_Numeric )
      END IF
      A % UMFPack_Numeric = 0
    END IF

    IF ( BigMode ) THEN
      ALLOCATE( LRows(SIZE(Rows)), LCols(SIZE(Cols)) )

      DO i=1,n+1
        LRows(i) = Rows(i)-1
      END DO

      DO i=1,SIZE(Cols)
        LCols(i) = Cols(i)-1
      END DO

      ln = n ! TODO: Kludge: ln is AddrInt and n is regular INTEGER
      CALL umf4_l_def( Control )
      CALL umf4_l_sym( ln,ln, LRows, LCols, Values, Symbolic, Control, iInfo )
    ELSE
      Rows = Rows-1
      Cols = Cols-1
      CALL umf4def( Control )
      CALL umf4sym( n,n, Rows, Cols, Values, Symbolic, Control, iInfo )
    END IF

    IF (iInfo(1)<0) THEN
      PRINT *, 'Error occurred in umf4sym: ', iInfo(1)
      STOP EXIT_ERROR
    END IF

    IF ( BigMode ) THEN
      CALL umf4_l_num(LRows, LCols, Values, Symbolic, A % UMFPack_Numeric, Control, iInfo )
    ELSE
      CALL umf4num( Rows, Cols, Values, Symbolic, A % UMFPack_Numeric, Control, iInfo )
    END IF

    IF (iinfo(1)<0) THEN
      PRINT*, 'Error occurred in umf4num: ', iinfo(1)
      STOP EXIT_ERROR
    ENDIF

    IF ( BigMode ) THEN
      DEALLOCATE( LRows, LCols )
      CALL umf4_l_fsym( Symbolic )
    ELSE
      A % Rows = A % Rows+1
      A % Cols = A % Cols+1
      CALL umf4fsym( Symbolic )
    END IF
  END IF

  IF ( BigMode ) THEN
    lsys = 2
    CALL umf4_l_sol( lsys, x, b, A % UMFPack_Numeric, Control, iInfo )
  ELSE
    sys = 2
    CALL umf4sol( sys, x, b, A % UMFPack_Numeric, Control, iInfo )
  END IF

  IF (iinfo(1)<0) THEN
    PRINT*, 'Error occurred in umf4sol: ', iinfo(1)
    STOP EXIT_ERROR
  END IF
 
  FreeFactorize = ListGetLogical( Solver % Values, &
      'Linear System Free Factorization', stat )
  IF ( .NOT. stat ) FreeFactorize = .TRUE.

  IF ( Factorize .AND. FreeFactorize ) THEN
    IF ( BigMode ) THEN
      CALL umf4_l_fnum(A % UMFPack_Numeric)
    ELSE
      CALL umf4fnum(A % UMFPack_Numeric)
    END IF
    A % UMFPack_Numeric = 0
  END IF
#else
   CALL Fatal( 'UMFPack_SolveSystem', 'UMFPACK Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE UMFPack_SolveSystem
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solves a linear system using Cholmod multifrontal direct solver courtesy
!> of University of Florida.
!------------------------------------------------------------------------------
  SUBROUTINE Cholmod_SolveSystem( Solver,A,x,b,Free_fact)
!------------------------------------------------------------------------------
  LOGICAL, OPTIONAL :: Free_Fact
  TYPE(Matrix_t) :: A
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: x(*), b(*)

  INTERFACE
     SUBROUTINE cholmod_ffree(chol) BIND(c,NAME="cholmod_ffree")
       USE Types
       INTEGER(KIND=AddrInt) :: chol
     END SUBROUTINE cholmod_ffree

     FUNCTION cholmod_ffactorize(n,rows,cols,vals,cmplx) RESULT(chol) BIND(c,NAME="cholmod_ffactorize")
        USE Types
        INTEGER :: n, cmplx, Rows(*), Cols(*)
        REAL(KIND=dp) :: Vals(*)
        INTEGER(KIND=dp) :: chol
     END FUNCTION cholmod_ffactorize

     SUBROUTINE cholmod_fsolve(chol, n, x,b) BIND(c,NAME="cholmod_fsolve")
        USE Types
        REAL(KIND=dp) :: x(*), b(*)
        INTEGER :: n
        INTEGER(KIND=dp) :: chol
     END SUBROUTINE cholmod_fsolve
  END INTERFACE

  LOGICAL :: Factorize, FreeFactorize, Found
  INTEGER :: i
  REAL(KIND=dp), POINTER CONTIG :: Vals(:)
  INTEGER, POINTER CONTIG :: Rows(:), Cols(:), Diag(:)

#ifdef HAVE_CHOLMOD
  IF ( PRESENT(Free_Fact) ) THEN
    IF ( Free_Fact ) THEN
      IF ( A % Cholmod/=0 ) THEN
        CALL cholmod_ffree(A % cholmod)
        A % cholmod = 0
      END IF
      RETURN
    END IF
  END IF

  Factorize = ListGetLogical( Solver % Values, &
     'Linear System Refactorize', Found )
  IF ( .NOT. Found ) Factorize = .TRUE.

  IF ( Factorize .OR. A% cholmod==0 ) THEN
    IF ( A % cholmod/=0 ) THEN
      CALL cholmod_ffree(A % cholmod)
      A % cholmod = 0
    END IF

    Rows => A % Rows
    Cols => A % Cols
    Vals => A % Values

    Rows=Rows-1; Cols=Cols-1 ! c numbering
    i=0
    IF(A % Complex) i=1
    A % Cholmod=cholmod_ffactorize(A % NumberOfRows, Rows, Cols, Vals, i)
    Rows=Rows+1; Cols=Cols+1 ! fortran numbering
  END IF

  CALL cholmod_fsolve(A % cholmod, A % NumberOfRows, x, b);

  FreeFactorize = ListGetLogical( Solver % Values, &
      'Linear System Free Factorization', Found )
  IF ( .NOT. Found ) FreeFactorize = .TRUE.

  IF ( Factorize .AND. FreeFactorize ) THEN
    CALL cholmod_ffree(A % cholmod)
    A % cholmod = 0
  END IF
#else
   CALL Fatal( 'Cholmod_SolveSystem', 'Cholmod Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE Cholmod_SolveSystem
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solves a linear system using SuiteSparseQR multifrontal direct solver courtesy
!> of University of Florida.
!------------------------------------------------------------------------------
  SUBROUTINE SPQR_SolveSystem( Solver,A,x,b,Free_fact)
!------------------------------------------------------------------------------
  LOGICAL, OPTIONAL :: Free_Fact
  TYPE(Matrix_t) :: A
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: x(*), b(*)

  INTEGER :: i
  LOGICAL :: Factorize, FreeFactorize, Found

  REAL(KIND=dp), POINTER CONTIG :: Vals(:)
  INTEGER, POINTER CONTIG :: Rows(:), Cols(:), Diag(:)

  INTERFACE
     FUNCTION spqr_ffree(chol) RESULT(stat) BIND(c,NAME="spqr_ffree")
       USE Types
       INTEGER :: stat
       INTEGER(KIND=AddrInt) :: chol
     END FUNCTION spqr_ffree

     FUNCTION spqr_ffactorize(n,rows,cols,vals) RESULT(chol) BIND(c,NAME="spqr_ffactorize")
       USE Types
       INTEGER :: n, rows(*), cols(*)
       REAL(KIND=dp) :: vals(*)
       INTEGER(KIND=AddrInt) :: chol
     END FUNCTION spqr_ffactorize

     SUBROUTINE spqr_fsolve(chol, n, x,b) BIND(c,NAME="spqr_fsolve")
       USE Types
       REAL(KIND=dp) :: x(*), b(*)
       INTEGER :: n
       INTEGER(KIND=AddrInt) :: chol
     END SUBROUTINE spqr_fsolve
  END INTERFACE

#ifdef HAVE_CHOLMOD
  IF ( PRESENT(Free_Fact) ) THEN
    IF ( Free_Fact ) THEN
      IF ( A % Cholmod/=0 ) THEN
        IF(spqr_ffree(A % cholmod)==0) A % Cholmod=0
      END IF
      RETURN
    END IF
  END IF

  Factorize = ListGetLogical( Solver % Values, &
     'Linear System Refactorize', Found )
  IF ( .NOT. Found ) Factorize = .TRUE.

  IF ( Factorize .OR. A% cholmod==0 ) THEN
    IF ( A % cholmod/=0 ) THEN
      i=spqr_ffree(A % cholmod)
      A % cholmod = 0
    END IF

    Rows => A % Rows
    Cols => A % Cols
    Vals => A % Values

    Rows=Rows-1; Cols=Cols-1 ! c numbering
    A % Cholmod=spqr_ffactorize(A % NumberOfRows, Rows, Cols, Vals)
    Rows=Rows+1; Cols=Cols+1 ! fortran numbering
  END IF

  CALL spqr_fsolve(A % cholmod, A % NumberOfRows, x, b);

  FreeFactorize = ListGetLogical( Solver % Values, &
      'Linear System Free Factorization', Found )
  IF ( .NOT. Found ) FreeFactorize = .TRUE.

  IF ( Factorize .AND. FreeFactorize ) THEN
    i=spqr_ffree(A % cholmod)
    A % cholmod = 0
  END IF
#else
   CALL Fatal( 'SPQR_SolveSystem', 'SPQR Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE SPQR_SolveSystem
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Stops the run with a readable message when a MUMPS call failed (status < 0) and reports
!> warnings (status > 0). Pass INFOG(1) and INFOG(2) for a parallel MUMPS instance: they are
!> identical on all ranks of its communicator, so all ranks take the same branch and no rank is
!> left waiting in a collective call. Pass INFO(1) and INFO(2) for an instance on one rank only.
!------------------------------------------------------------------------------
  SUBROUTINE CheckMumpsStatus(Caller, Phase, Label, Status, Detail)
!------------------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(IN) :: Caller, Phase, Label
    INTEGER, INTENT(IN) :: Status, Detail
    CHARACTER(LEN=160) :: Hint
!------------------------------------------------------------------------------
    IF (Status > 0) THEN
      WRITE(Message, '(5A,I0,3A,I0)') 'MUMPS ', Phase, ' finished with a warning: ', Label, '(1) = ', Status, &
          ', ', Label, '(2) = ', Detail
      CALL Info(Caller, Message, Level=4)
    END IF
    IF (Status >= 0) RETURN

    SELECT CASE (Status)
    CASE (-5, -7, -13)
      Hint = 'memory allocation failed (out of memory)'
    CASE (-8, -9, -11, -12, -14, -15, -17, -19, -20)
      Hint = 'internal workspace too small, which a singular matrix can also cause; check the model or ' // &
          'increase "Mumps Percentage Increase Working Space"'
    CASE (-10)
      Hint = 'numerically singular matrix, check the boundary conditions or use "Mumps Null Pivot Detection"'
    CASE (-6)
      Hint = 'structurally singular matrix'
    CASE DEFAULT
      Hint = 'see the MUMPS users'' guide'
    END SELECT
    WRITE(Message, '(5A,I0,3A,I0,2A)') 'MUMPS ', Phase, ' failed: ', Label, '(1) = ', Status, &
        ', ', Label, '(2) = ', Detail, ': ', TRIM(Hint)
    CALL Fatal(Caller, Message)
!------------------------------------------------------------------------------
  END SUBROUTINE CheckMumpsStatus
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Diagnostic, switched on with 'Mumps Residual Check = True'. Computes r = b - A x for the system handed to
!> MUMPS, summed over all ranks in the same way as MUMPS sums the distributed entries. For the rows of the
!> solver's own matrix (field) and for the rows appended to it (constraints, e.g. circuit equations) it
!> reports separately: max|r_i|, max|b_i| (the right-hand side of the current, possibly linearized, system),
!> the largest row activity max_i (|A||x|)_i, the normwise residual max|r| / (max|A||x| + max|b|) and a
!> residual indicator max_i |r_i| / (|A||x| + |b|)_i. The warning compares the normwise residual: on 156
!> TRAFOLO A-V runs it stayed at or below 1.7e-10 for correct solutions, while the indicator reached 0.9.
!> Neither ratio is the exact backward error of the assembled system: on rows shared by ranks, |A| sums the
!> moduli of the rank-local partial entries, which can exceed the modulus of the assembled entry when partial
!> entries cancel, so both ratios can understate the residual; on one rank this effect vanishes. Rows below
!> 1e-12 of the largest activity in their block are measured against that floor (their number is reported).
!> Before the solve, x is the current iterate, so r is the nonlinear residual of that iterate; after the
!> solve, r is the residual of the linear solve. Works on the real form of complex systems.
!------------------------------------------------------------------------------
  SUBROUTINE MumpsResidualCheck(Solver, A, x, b, Comm, What, Level, WarnAbove)
!------------------------------------------------------------------------------
#ifdef HAVE_MUMPS
#  if defined(ELMER_HAVE_MPI_MODULE)
    USE mpi
#  endif
#endif
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp) :: x(*), b(*)
    INTEGER :: Comm
    CHARACTER(LEN=*) :: What
    !> Info level of the report (default 4); normwise residual above which to warn.
    INTEGER, OPTIONAL :: Level
    REAL(KIND=dp), OPTIONAL :: WarnAbove
!------------------------------------------------------------------------------
#ifdef HAVE_MUMPS
#  if defined(ELMER_HAVE_MPIF_HEADER)
    INCLUDE 'mpif.h'
#  endif
    INTEGER :: i, j, k, l, n, nfield, ng, nloc, ierr
    REAL(KIND=dp), ALLOCATABLE :: xg(:), cnt(:), rg(:), ag(:), bg(:), cg(:), buf(:)
    REAL(KIND=dp) :: s, t, den, rmax(2), bmax(2), amax(2), dmax(2), dfloor(2), wmax(2), nrel(2)
    INTEGER :: nfloor(2)
    LOGICAL :: GotName

    n = A % NumberOfRows
    nfield = n
    IF (ASSOCIATED(Solver % Matrix)) nfield = MIN(n, Solver % Matrix % NumberOfRows)
    nloc = MAXVAL(A % Gorder(1:n))
    CALL MPI_ALLREDUCE(nloc, ng, 1, MPI_INTEGER, MPI_MAX, Comm, ierr)
    ALLOCATE(xg(ng), cnt(ng), rg(ng), ag(ng), bg(ng), cg(ng), buf(ng))
    xg = 0.0_dp; cnt = 0.0_dp; rg = 0.0_dp; ag = 0.0_dp; bg = 0.0_dp; cg = 0.0_dp
    DO i = 1, n
      k = A % Gorder(i)
      xg(k) = x(i)
      cnt(k) = 1.0_dp
      IF (i > nfield) cg(k) = 1.0_dp
    END DO
    ! shared rows carry the same solution value on every rank: average the sum
    buf = xg
    CALL MPI_ALLREDUCE(buf, xg, ng, MPI_DOUBLE_PRECISION, MPI_SUM, Comm, ierr)
    buf = cnt
    CALL MPI_ALLREDUCE(buf, cnt, ng, MPI_DOUBLE_PRECISION, MPI_SUM, Comm, ierr)
    buf = cg
    CALL MPI_ALLREDUCE(buf, cg, ng, MPI_DOUBLE_PRECISION, MPI_MAX, Comm, ierr)
    WHERE (cnt > 0.0_dp) xg = xg / cnt
    DO i = 1, n
      s = 0.0_dp
      t = 0.0_dp
      DO j = A % Rows(i), A % Rows(i+1) - 1
        s = s + A % Values(j) * xg(A % Gorder(A % Cols(j)))
        t = t + ABS(A % Values(j)) * ABS(xg(A % Gorder(A % Cols(j))))
      END DO
      k = A % Gorder(i)
      rg(k) = rg(k) + b(i) - s
      ag(k) = ag(k) + t
      bg(k) = bg(k) + b(i)
    END DO
    buf = rg
    CALL MPI_ALLREDUCE(buf, rg, ng, MPI_DOUBLE_PRECISION, MPI_SUM, Comm, ierr)
    buf = ag
    CALL MPI_ALLREDUCE(buf, ag, ng, MPI_DOUBLE_PRECISION, MPI_SUM, Comm, ierr)
    buf = bg
    CALL MPI_ALLREDUCE(buf, bg, ng, MPI_DOUBLE_PRECISION, MPI_SUM, Comm, ierr)
    ! ag holds (|A||x|)_i; on rows shared by ranks it sums the moduli of the partial entries (an upper bound)
    rmax = 0.0_dp; bmax = 0.0_dp; amax = 0.0_dp; dmax = 0.0_dp; wmax = 0.0_dp; nfloor = 0
    DO k = 1, ng
      IF (cnt(k) <= 0.0_dp) CYCLE
      l = 1
      IF (cg(k) > 0.5_dp) l = 2
      rmax(l) = MAX(rmax(l), ABS(rg(k)))
      bmax(l) = MAX(bmax(l), ABS(bg(k)))
      amax(l) = MAX(amax(l), ag(k))
      dmax(l) = MAX(dmax(l), ag(k) + ABS(bg(k)))
    END DO
    dfloor = 1.0e-12_dp * dmax
    DO k = 1, ng
      IF (cnt(k) <= 0.0_dp) CYCLE
      l = 1
      IF (cg(k) > 0.5_dp) l = 2
      den = ag(k) + ABS(bg(k))
      IF (den < dfloor(l)) THEN
        nfloor(l) = nfloor(l) + 1
        den = dfloor(l)
      END IF
      IF (den > 0.0_dp) wmax(l) = MAX(wmax(l), ABS(rg(k)) / den)
    END DO
    nrel = 0.0_dp
    DO l = 1, 2
      IF (amax(l) + bmax(l) > 0.0_dp) nrel(l) = rmax(l) / (amax(l) + bmax(l))
    END DO
    WRITE(Message, '(2A,2(A,5(ES10.3,A),I0,A),A,I0,A)') 'MUMPS residual check ', What, &
        ': field rows: max|r| ', rmax(1), ', max|b| ', bmax(1), ', max|A||x| ', amax(1), &
        ', normwise ', nrel(1), ', indicator ', wmax(1), ' (', nfloor(1), ' floored)', &
        '; constraint rows: max|r| ', rmax(2), ', max|b| ', bmax(2), ', max|A||x| ', amax(2), &
        ', normwise ', nrel(2), ', indicator ', wmax(2), ' (', nfloor(2), ' floored)', &
        '; ', NINT(SUM(cg)), ' constraint rows'
    IF (PRESENT(Level)) THEN
      CALL Info('MumpsResidualCheck', Message, Level=Level)
    ELSE
      CALL Info('MumpsResidualCheck', Message, Level=4)
    END IF
    ! MUMPS reports success on a singular system it factorized without null
    ! pivot detection; only the residual shows that the solution is wrong.
    IF (PRESENT(WarnAbove)) THEN
      IF (MAXVAL(nrel) > WarnAbove) THEN
        WRITE(Message, '(A,ES10.3,A,ES10.3,A)') 'NOT CONVERGED: direct solver="'// &
            ListGetString(Solver % Values,'Equation',GotName)//'" MUMPS normwise residual=', MAXVAL(nrel), &
            ' tolerance=', WarnAbove, ': the solution does not satisfy the linear system. '// &
            'A singular system needs "Mumps Null Pivot Detection = True"; the A-V solvers with '// &
            'circuit coils also need "Use Tree Gauge = False".'
        CALL Warn('MumpsResidualCheck', Message)
      END IF
    END IF
    DEALLOCATE(xg, cnt, rg, ag, bg, cg, buf)
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE MumpsResidualCheck
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> The residual of every MUMPS solution, unless 'Mumps Residual Check = False':
!> reported at info level 4 when the check was asked for, at level 6
!> otherwise, and a warning when the normwise residual exceeds 'Mumps Residual
!> Tolerance' (default 1e-8). Collective over Comm, so every rank calls it.
!------------------------------------------------------------------------------
  SUBROUTINE MumpsCheckSolution(Solver, A, x, b, Comm)
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp) :: x(*), b(*)
    INTEGER :: Comm

    ! Correct TRAFOLO A-V solutions stayed at or below 1.7e-10; failed solves reached 5e-2.
    REAL(KIND=dp), PARAMETER :: DefaultResidualTolerance = 1.0e-8_dp
    REAL(KIND=dp) :: Tol
    LOGICAL :: Check, Asked, GotTol

    Check = ListGetLogical(Solver % Values, 'Mumps Residual Check', Asked)
    IF (.NOT. Asked) Check = .TRUE.
    IF (.NOT. Check) RETURN
    Tol = ListGetConstReal(Solver % Values, 'Mumps Residual Tolerance', GotTol)
    IF (.NOT. GotTol) Tol = DefaultResidualTolerance
    IF (Asked) THEN
      CALL MumpsResidualCheck(Solver, A, x, b, Comm, 'after the solve', Level=4, WarnAbove=Tol)
    ELSE
      CALL MumpsResidualCheck(Solver, A, x, b, Comm, 'after the solve', Level=6, WarnAbove=Tol)
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE MumpsCheckSolution
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE FreeMumpsFactorizations(A)
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A

!------------------------------------------------------------------------------
#ifdef HAVE_MUMPS
    IF(ASSOCIATED(A % SMumpsID)) THEN
      DEALLOCATE( A % SMumpsID % irn_loc, &
         A % SMumpsID % jcn_loc, A % SMumpsID % Rhs,  &
           A % SMumpsID % isol_loc, A % SMumpsID % sol_loc)

      DEALLOCATE( A % SMumpsID % A_loc )
      IF (ASSOCIATED(A % Gorder)) DEALLOCATE(A % Gorder)
      A % Gorder=>Null()

      A % SMumpsID % job = -2
      CALL SMumps(A % SMumpsID)
      DEALLOCATE(A % SMumpsID)
      A % SMumpsId => NULL()
    END IF

    IF(ASSOCIATED(A % CMumpsID)) THEN
      DEALLOCATE( A % CMumpsID % irn_loc, &
         A % CMumpsID % jcn_loc, A % CMumpsID % Rhs,  &
           A % CMumpsID % isol_loc, A % CMumpsID % sol_loc)

      DEALLOCATE( A % CMumpsID % A_loc )
      IF (ASSOCIATED(A % Gorder)) DEALLOCATE(A % Gorder)
      A % Gorder=>Null()

      A % CMumpsID % job = -2
      CALL CMumps(A % CMumpsID)
      DEALLOCATE(A % CMumpsID)
      A % CMumpsId => NULL()
    END IF

    IF(ASSOCIATED(A % MumpsID)) THEN
      DEALLOCATE( A % MumpsID % irn_loc, &
         A % MumpsID % jcn_loc, A % MumpsID % Rhs,  &
           A % MumpsID % isol_loc, A % MumpsID % sol_loc)

      IF(.NOT.ASSOCIATED(A % MumpsID % A_loc, A % Values)) DEALLOCATE( A % MumpsID % A_loc )
      IF (ASSOCIATED(A % Gorder)) DEALLOCATE(A % Gorder)
      A % Gorder=>Null()

      A % MumpsID % job = -2
      CALL DMumps(A % MumpsID)
      DEALLOCATE(A % MumpsID)
      A % MumpsId => NULL()
    END IF

    IF(ASSOCIATED(A % ZMumpsID)) THEN
      DEALLOCATE( A % ZMumpsID % irn_loc, &
         A % ZMumpsID % jcn_loc, A % ZMumpsID % Rhs,  &
           A % ZMumpsID % isol_loc, A % ZMumpsID % sol_loc)

      DEALLOCATE( A % ZMumpsID % A_loc )

      IF (ASSOCIATED(A % Gorder)) DEALLOCATE(A % Gorder)
      A % Gorder=>Null()

      A % ZMumpsID % job = -2
      CALL ZMumps(A % ZMumpsID)
      DEALLOCATE(A % ZMumpsID)
      A % ZMumpsId => NULL()
    END IF
#endif
!------------------------------------------------------------------------------
 END SUBROUTINE FreeMumpsFactorizations
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solves a linear system using MUMPS direct solver. This is a legacy solver
!> with complicated dependencies. Single precision version.
!------------------------------------------------------------------------------
  SUBROUTINE SMumps_SolveSystem( Solver,A,x,b )
!------------------------------------------------------------------------------
#ifdef HAVE_MUMPS
#  if defined(ELMER_HAVE_MPI_MODULE)
  USE mpi
#  endif
#endif

  TYPE(Matrix_t) :: A
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp), TARGET :: x(*), b(*)

#ifdef HAVE_MUMPS
#  if defined(ELMER_HAVE_MPIF_HEADER)
  INCLUDE 'mpif.h'
#  endif

  INTEGER, ALLOCATABLE :: Owner(:)
  INTEGER :: i,j,n,ip,ierr,icntlft,nzloc
  REAL(KIND=dp) :: nullpivtol
  LOGICAL :: Factorize, FreeFactorize, stat, matsym, matspd, scaled

  INTEGER, ALLOCATABLE :: memb(:)
  INTEGER :: Comm_active, Group_active, Group_world

  REAL, ALLOCATABLE :: dbuf(:)

  Factorize = ListGetLogical( Solver % Values, 'Linear System Refactorize', stat )
  IF ( .NOT. stat ) Factorize = .TRUE.

  IF ( Factorize .OR. .NOT.ASSOCIATED(A % SMumpsID) ) THEN
    CALL FreeMumpsFactorizations(A)
    ALLOCATE(A % SMumpsID)

    A % SMumpsID % Comm = A % Comm
    A % SMumpsID % par  =  1
    A % SMumpsID % job  = -1
    A % SMumpsID % Keep =  0

    matsym = ListGetLogical( Solver % Values, 'Linear System Symmetric', stat)
    matspd = ListGetLogical( Solver % Values, 'Linear System Positive Definite', stat)

    ! force unsymmetric mode when "row equilibration" is used
    scaled = ListGetLogical( Solver % Values, 'Linear System Scaling', stat)
    IF(.NOT.stat) scaled = .TRUE.
    IF(scaled) THEN
      IF(ListGetLogical( Solver % Values, 'Linear System Row Equilibration',stat)) matsym=.FALSE.
    END IF

    IF(matsym) THEN
      IF ( matspd) THEN
        A % MumpsID % sym = 1
      ELSE
        A % SMumpsID % sym = 0 ! 2=symmetric, but unsymmetric solver seems faster, at least in a few
                              ! simple cases...  more testing needed...
      END IF
    ELSE
      A % SMumpsID % sym = 0
    END IF

    CALL SMumps(A % SMumpsID)
    CALL CheckMumpsStatus('SMumps_SolveSystem', 'initialization', 'INFOG', &
        A % SMumpsID % infog(1), A % SMumpsID % infog(2))

    IF(ASSOCIATED(A % Gorder)) DEALLOCATE(A % Gorder)

    IF(ASSOCIATED(A % ParallelInfo)) THEN
      n = SIZE(A % ParallelInfo % GlobalDOFs)

      ALLOCATE( A % Gorder(n), Owner(n) )
      CALL ContinuousNumbering( A % ParallelInfo, A % Perm, A % Gorder, Owner )

      CALL MPI_ALLREDUCE( SUM(Owner), A % SMumpsID % n, &
         1, MPI_INTEGER, MPI_SUM, A % SMumpsID % Comm, ierr )
      DEALLOCATE(Owner)
    ELSE
      CALL MPI_ALLREDUCE( A % NumberOfRows, A % SMumpsId % n, &
            1, MPI_INTEGER, MPI_MAX, A % Comm, ierr )

      ALLOCATE(A % Gorder(A % NumberOFrows))
      DO i=1,A % NumberOFRows
        A % Gorder(i) = i
      END DO
    END IF

   ! Set matrix for Mumps (unsymmetric case)
    IF (A % SmumpsID % sym == 0) THEN
      A % SMumpsID % nz_loc = A % Rows(A % NumberOfRows+1)-1

      ALLOCATE( A % SMumpsID % irn_loc(A % SMumpsID % nz_loc) )
      ALLOCATE( A % SMumpsID % a_loc(A % SMumpsId % nz_loc) )
      ALLOCATE( A % SMumpsID % jcn_loc(A % SMumpsId % nz_loc) )

      nzloc = 0
      DO i=1,A % NumberOfRows
        ip = A % Gorder(i)
        DO j=A % Rows(i),A % Rows(i+1)-1
          nzloc = nzloc + 1
          A % SMumpsID % irn_loc(nzloc) = ip
          A % SMumpsID % a_loc(nzloc)   = A % Values(j)
          A % SMumpsID % jcn_loc(nzloc) = A % Gorder(A % Cols(j))
        END DO
      END DO
    ELSE
      ! Set matrix for Mumps (symmetric case)
      nzloc = 0
      DO i=1,A % NumberOfRows
        ! Only output lower triangular part to Mumps
        DO j=A % Rows(i),A % Diag(i)
          nzloc = nzloc + 1
        END DO
      END DO

      A % SMumpsID % nz_loc = nzloc

      ALLOCATE( A % SMumpsID % irn_loc(A % SMumpsID % nz_loc) )
      ALLOCATE( A % SMumpsID % jcn_loc(A % SMumpsId % nz_loc) )
      ALLOCATE( A % SMumpsID % A_loc(A % SMumpsId % nz_loc) )

      nzloc = 0
      DO i=1,A % NumberOfRows
        ! Only output lower triangular part to Mumps
        ip = A % Gorder(i)
        DO j=A % Rows(i),A % Diag(i)
          nzloc = nzloc + 1
          A % SmumpsID % IRN_loc(nzloc) = ip
          A % SmumpsID % A_loc(nzloc) = A % Values(j)
          A % SmumpsID % JCN_loc(nzloc) = A % Gorder(A % Cols(j))
        END DO
      END DO
    END IF


    ALLOCATE(A % SMumpsID % rhs(A % SMumpsId % n))

    ! Tune verbosity of MUMPS.
    i = 0
    IF(InfoActive(20)) i = 1
    A % SMumpsID % icntl(2) = i ! suppress printing of diagnostics and warnings
    A % SMumpsID % icntl(3) = i ! suppress statistics

    A % SMumpsID % icntl(4) = 1 ! the same as the two above, but doesn't seem to work.
    A % SMumpsID % icntl(5) = 0 ! matrix format 'assembled'

    icntlft = ListGetInteger(Solver % Values, &
          'mumps percentage increase working space', stat)
    IF (stat) THEN
       A % SMumpsID % icntl(14) = icntlft
    END IF
    IF (ListGetLogical(Solver % Values, 'mumps null pivot detection', stat)) THEN
       A % SMumpsID % icntl(24) = 1
       CALL Info('SMumps_SolveSystem', 'MUMPS null pivot detection enabled', Level=5)
    END IF
    ! CNTL(3) null-pivot threshold: >0 CNTL(3)*||A||inf, <0 |CNTL(3)|, 0 MUMPS default
    ! (4.10.0: 1e-5*eps*||A||inf; 5.4.0 and later: sqrt(pivots on critical path)*eps*||A||inf). Used: DKEEP(1).
    nullpivtol = ListGetConstReal(Solver % Values, 'mumps null pivot tolerance', stat)
    IF (stat) THEN
       A % SMumpsID % cntl(3) = nullpivtol
    END IF
    IF (ListGetLogical(Solver % Values, 'mumps sequential root', stat)) THEN
       A % SMumpsID % icntl(13) = 1
       CALL Info('SMumps_SolveSystem', 'MUMPS sequential root node (no ScaLAPACK) enabled', Level=5)
    END IF
    A % SMumpsID % icntl(18) = 3 ! 'distributed' matrix 
    A % SMumpsID % icntl(21) = 1 ! 'distributed' solution phase

    A % SMumpsID % job = 4
    CALL SMumps(A % SMumpsID)
    CALL Flush(6)
    CALL CheckMumpsStatus('SMumps_SolveSystem', 'analysis and factorization', 'INFOG', &
        A % SMumpsID % infog(1), A % SMumpsID % infog(2))
    IF (A % SMumpsID % icntl(24) == 1) THEN
       CALL Info('SMumps_SolveSystem', 'MUMPS null pivots detected: '//I2S(A % SMumpsID % infog(28)), Level=5)
    END IF

    A % SMumpsID % lsol_loc = A % Smumpsid % info(23)
    ALLOCATE(A % SMumpsID % sol_loc(A % SMumpsId % lsol_loc))
    ALLOCATE(A % SMumpsID % isol_loc(A % SMumpsId % lsol_loc))
  END IF

 ! sum the rhs from all procs. Could be done for neighbours only (i guess):
 ! ------------------------------------------------------------------------
  A % SMumpsID % RHS = 0
  DO i=1,A % NumberOfRows
    ip = A % Gorder(i)
    A % SMumpsId % RHS(ip) = b(i)
  END DO

  ALLOCATE( dbuf(A % SMumpsID % n) )
  dbuf = A % SMumpsId % RHS
  CALL MPI_ALLREDUCE( dbuf, A % SMumpsID % RHS, &
    A % SMumpsID % n, MPI_REAL, MPI_SUM, A % SMumpsID % Comm, ierr )

 ! Solution:
 ! ---------
  A % SMumpsID % job = 3
  CALL SMumps(A % SMumpsID)
  CALL CheckMumpsStatus('SMumps_SolveSystem', 'solve', 'INFOG', &
      A % SMumpsID % infog(1), A % SMumpsID % infog(2))

 ! Distribute the solution to all:
 ! -------------------------------
  A % SMumpsId % Rhs = 0
  DO i=1,A % SMumpsID % lsol_loc
    A % SMumpsID % RHS(A % SMumpsID % isol_loc(i)) = A % SMumpsID % sol_loc(i)
  END DO
  dbuf = A % SMumpsId % RHS
  CALL MPI_ALLREDUCE( dbuf, A % SMumpsID % RHS, &
    A % SMumpsID % n, MPI_REAL, MPI_SUM,A %  SMumpsID % Comm, ierr )

  DEALLOCATE(dbuf)

 ! Select the values which belong to us:
 ! -------------------------------------
  DO i=1,A % NumberOfRows
    ip = A % Gorder(i)
    x(i) = A % SMumpsId % RHS(ip)
  END DO

  FreeFactorize = ListGetLogical( Solver % Values, 'Linear System Free Factorization', stat )
  IF ( .NOT. stat ) FreeFactorize = .TRUE.

  IF ( Factorize .AND. FreeFactorize ) CALL FreeMumpsFactorizations(A)
#else
   CALL Fatal( 'Mumps_SolveSystem', 'MUMPS Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE SMumps_SolveSystem
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solves a linear system using MUMPS direct solver. This is a legacy solver
!> with complicated dependencies. Single precision complex version.
!------------------------------------------------------------------------------
  SUBROUTINE CMumps_SolveSystem( Solver,A,x,b )
!------------------------------------------------------------------------------
#ifdef HAVE_MUMPS
#  if defined(ELMER_HAVE_MPI_MODULE)
  USE mpi
#  endif
#endif

  TYPE(Matrix_t) :: A
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp), TARGET :: x(*), b(*)

#ifdef HAVE_MUMPS
#  if defined(ELMER_HAVE_MPIF_HEADER)
  INCLUDE 'mpif.h'
#  endif

  INTEGER, ALLOCATABLE :: Owner(:)
  INTEGER :: i,j,n,ip,ierr,icntlft,nzloc
  REAL(KIND=dp) :: nullpivtol
  LOGICAL :: Factorize, FreeFactorize, stat, matsym, matspd, scaled

  INTEGER, ALLOCATABLE :: memb(:)
  INTEGER :: Comm_active, Group_active, Group_world

  COMPLEX, ALLOCATABLE :: dbuf(:)

  Factorize = ListGetLogical( Solver % Values, 'Linear System Refactorize', stat )
  IF ( .NOT. stat ) Factorize = .TRUE.

  IF ( Factorize .OR. .NOT.ASSOCIATED(A % CMumpsID) ) THEN
    CALL FreeMumpsFactorizations(A)
    ALLOCATE(A % CMumpsID)

    A % CMumpsID % Comm = A % Comm
    A % CMumpsID % par  =  1
    A % CMumpsID % job  = -1
    A % CMumpsID % Keep =  0

!   matsym = ListGetLogical( Solver % Values, 'Linear System Symmetric', stat)
!   matspd = ListGetLogical( Solver % Values, 'Linear System Positive Definite', stat)
    matsym = .FALSE.
    matspd = .FALSE.

    ! force unsymmetric mode when "row equilibration" is used
    scaled = ListGetLogical( Solver % Values, 'Linear System Scaling', stat)
    IF(.NOT.stat) scaled = .TRUE.
    IF(scaled) THEN
      IF(ListGetLogical( Solver % Values, 'Linear System Row Equilibration',stat)) matsym=.FALSE.
    END IF

!   IF(matsym) THEN
!     IF ( matspd) THEN
!       A % CMumpsID % sym = 1
!     ELSE
!       A % CMumpsID % sym = 0 ! 2=symmetric, but unsymmetric solver seems faster, at least in a few
!                             ! simple cases...  more testing needed...
!     END IF
!   ELSE
!     A % CMumpsID % sym = 0
!   END IF
     A % CMumpsID % sym = 0

    CALL CMumps(A % CMumpsID)
    CALL CheckMumpsStatus('CMumps_SolveSystem', 'initialization', 'INFOG', &
        A % CMumpsID % infog(1), A % CMumpsID % infog(2))

    IF(ASSOCIATED(A % Gorder)) DEALLOCATE(A % Gorder)

    IF(ASSOCIATED(A % ParallelInfo)) THEN
      n = SIZE(A % ParallelInfo % GlobalDOFs)
      ALLOCATE( A % Gorder(n), Owner(n) )
      CALL ContinuousNumbering( A % ParallelInfo, A % Perm, A % Gorder, Owner )

      CALL MPI_ALLREDUCE( SUM(Owner)/2, A % CMumpsID % n, &
        1, MPI_INTEGER, MPI_SUM, A % CMumpsID % Comm, ierr )
      DEALLOCATE(Owner)
    ELSE
      CALL MPI_ALLREDUCE( A % NumberOfRows/2, A % CMumpsId % n, &
            1, MPI_INTEGER, MPI_MAX, A % Comm, ierr )

      ALLOCATE(A % Gorder(A % NumberOFrows))
      DO i=1,A % NumberOFRows
        A % Gorder(i) = i
      END DO
    END IF

   ! Set matrix for Mumps (unsymmetric case)
    IF (A % CmumpsID % sym == 0) THEN
      A % CMumpsID % nz_loc = (A % Rows(A % NumberOfRows+1)-1)/4

      ALLOCATE( A % CMumpsID % irn_loc(A % CMumpsID % nz_loc) )
      ALLOCATE( A % CMumpsID % a_loc(A % CMumpsId % nz_loc) )
      ALLOCATE( A % CMumpsID % jcn_loc(A % CMumpsId % nz_loc) )

      nzloc = 0
      DO i=1,A % NumberOfRows,2
        ip = (A % Gorder(i)-1)/2+1
        DO j=A % Rows(i),A % Rows(i+1)-1,2
          nzloc = nzloc + 1
          A % CMumpsID % irn_loc(nzloc) = ip
          A % CMumpsID % jcn_loc(nzloc) = (A % Gorder(A % Cols(j))-1)/2+1
          A % CMumpsID % a_loc(nzloc)   = CMPLX( A % Values(j), -A % Values(j+1) )
        END DO
      END DO
    ELSE
      ! Set matrix for Mumps (symmetric case)
      nzloc = 0
      DO i=1,A % NumberOfRows
        ! Only output lower triangular part to Mumps
        DO j=A % Rows(i),A % Diag(i)
          nzloc = nzloc + 1
        END DO
      END DO

      A % CMumpsID % nz_loc = nzloc

      ALLOCATE( A % CMumpsID % irn_loc(A % CMumpsID % nz_loc) )
      ALLOCATE( A % CMumpsID % jcn_loc(A % CMumpsId % nz_loc) )
      ALLOCATE( A % CMumpsID % A_loc(A % CMumpsId % nz_loc) )

      nzloc = 0
      DO i=1,A % NumberOfRows,2
        ! Only output lower triangular part to Mumps
        ip = A % Gorder(i)
        DO j=A % Rows(i),A % Diag(i),2
          nzloc = nzloc + 1
          A % CmumpsID % IRN_loc(nzloc) = ip
          A % CmumpsID % A_loc(nzloc) = A % Values(j)
          A % CmumpsID % JCN_loc(nzloc) = A % Gorder(A % Cols(j))
        END DO
      END DO
    END IF


    ALLOCATE(A % CMumpsID % rhs(A % CMumpsId % n))

    ! Tune verbosity of MUMPS.
    i = 0
    IF(InfoActive(20)) i = 1
    A % CMumpsID % icntl(2) = i ! suppress printing of diagnostics and warnings
    A % CMumpsID % icntl(3) = i ! suppress statistics

    A % CMumpsID % icntl(4) = 1 ! the same as the two above, but doesn't seem to work.
    A % CMumpsID % icntl(5) = 0 ! matrix format 'assembled'

    icntlft = ListGetInteger(Solver % Values, &
          'mumps percentage increase working space', stat)
    IF (stat) THEN
       A % CMumpsID % icntl(14) = icntlft
    END IF
    IF (ListGetLogical(Solver % Values, 'mumps null pivot detection', stat)) THEN
       A % CMumpsID % icntl(24) = 1
       CALL Info('CMumps_SolveSystem', 'MUMPS null pivot detection enabled', Level=5)
    END IF
    ! CNTL(3) null-pivot threshold: >0 CNTL(3)*||A||inf, <0 |CNTL(3)|, 0 MUMPS default
    ! (4.10.0: 1e-5*eps*||A||inf; 5.4.0 and later: sqrt(pivots on critical path)*eps*||A||inf). Used: DKEEP(1).
    nullpivtol = ListGetConstReal(Solver % Values, 'mumps null pivot tolerance', stat)
    IF (stat) THEN
       A % CMumpsID % cntl(3) = nullpivtol
    END IF
    IF (ListGetLogical(Solver % Values, 'mumps sequential root', stat)) THEN
       A % CMumpsID % icntl(13) = 1
       CALL Info('CMumps_SolveSystem', 'MUMPS sequential root node (no ScaLAPACK) enabled', Level=5)
    END IF
    A % CMumpsID % icntl(18) = 3 ! 'distributed' matrix 
    A % CMumpsID % icntl(21) = 1 ! 'distributed' solution phase

    A % CMumpsID % job = 4
    CALL CMumps(A % CMumpsID)
    CALL Flush(6)
    CALL CheckMumpsStatus('CMumps_SolveSystem', 'analysis and factorization', 'INFOG', &
        A % CMumpsID % infog(1), A % CMumpsID % infog(2))
    IF (A % CMumpsID % icntl(24) == 1) THEN
       CALL Info('CMumps_SolveSystem', 'MUMPS null pivots detected: '//I2S(A % CMumpsID % infog(28)), Level=5)
    END IF

    A % CMumpsID % lsol_loc = A % CMumpsid % info(23)
    ALLOCATE(A % CMumpsID % sol_loc(A % CMumpsId % lsol_loc))
    ALLOCATE(A % CMumpsID % isol_loc(A % CMumpsId % lsol_loc))
  END IF

 ! sum the rhs from all procs. Could be done
 ! for neighbours only (i guess):
 ! ------------------------------------------
  A % CMumpsID % RHS = 0
  DO i=1,A % NumberOfRows,2
    ip = (A % Gorder(i)-1)/2+1
    A % CMumpsId % RHS(ip) = CMPLX( b(i), b(i+1) )
  END DO

  ALLOCATE( dbuf(A % CMumpsID % n) )
  dbuf = A % CMumpsId % RHS
  CALL MPI_ALLREDUCE( dbuf, A % CMumpsID % RHS, &
    A % CMumpsID % n, MPI_COMPLEX, MPI_SUM, A % CMumpsID % Comm, ierr )

 ! Solution:
 ! ---------
  A % CMumpsID % job = 3
  CALL CMumps(A % CMumpsID)
  CALL CheckMumpsStatus('CMumps_SolveSystem', 'solve', 'INFOG', &
      A % CMumpsID % infog(1), A % CMumpsID % infog(2))

 ! Distribute the solution to all:
 ! -------------------------------
  A % CMumpsId % Rhs = 0
  DO i=1,A % CMumpsID % lsol_loc
    A % CMumpsID % RHS(A % CMumpsID % isol_loc(i)) = A % CMumpsID % sol_loc(i)
  END DO
  dbuf = A % CMumpsId % RHS
  CALL MPI_ALLREDUCE( dbuf, A % CMumpsID % RHS, &
    A % CMumpsID % N, MPI_COMPLEX, MPI_SUM, A %  CMumpsID % Comm, ierr )

  DEALLOCATE(dbuf)

 ! Select the values which belong to us:
 ! -------------------------------------
  DO i=1,A % NumberOfRows,2
    ip = (A % Gorder(i)-1)/2+1
    x(i)   = REAL( A % CMumpsId % RHS(ip) )
    x(i+1) = AIMAG( A % CMumpsId % RHS(ip) )
  END DO

  FreeFactorize = ListGetLogical( Solver % Values, 'Linear System Free Factorization', stat )
  IF ( .NOT. stat ) FreeFactorize = .TRUE.
  IF ( Factorize .AND. FreeFactorize ) CALL FreeMumpsFactorizations(A)

#else
   CALL Fatal( 'Mumps_SolveSystem', 'MUMPS Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE CMumps_SolveSystem
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Solves a linear system using MUMPS direct solver. This is a legacy solver
!> with complicated dependencies. This is only available in parallel. 
!------------------------------------------------------------------------------
  SUBROUTINE Mumps_SolveSystem( Solver,A,x,b )
!------------------------------------------------------------------------------
#ifdef HAVE_MUMPS
#  if defined(ELMER_HAVE_MPI_MODULE)
  USE mpi
#  endif
#endif

  TYPE(Matrix_t) :: A
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp), TARGET :: x(*), b(*)

#ifdef HAVE_MUMPS
#  if defined(ELMER_HAVE_MPIF_HEADER)
  INCLUDE 'mpif.h'
#  endif

  INTEGER, ALLOCATABLE :: Owner(:)
  INTEGER :: i,j,n,ip,ierr,icntlft,nzloc
  REAL(KIND=dp) :: nullpivtol
  LOGICAL :: Factorize, FreeFactorize, stat, matsym, matspd, scaled

  INTEGER, ALLOCATABLE :: memb(:)
  INTEGER :: Comm_active, Group_active, Group_world

  REAL(KIND=dp), ALLOCATABLE :: dbuf(:)


  Factorize = ListGetLogical( Solver % Values, 'Linear System Refactorize', stat )
  IF ( .NOT. stat ) Factorize = .TRUE.

  IF ( Factorize .OR. .NOT.ASSOCIATED(A % MumpsID) ) THEN
    CALL FreeMumpsFactorizations(A)
    ALLOCATE(A % MumpsID)

    A % MumpsID % Comm = A % Comm
    A % MumpsID % par  =  1
    A % MumpsID % job  = -1
    A % MumpsID % Keep =  0

    matsym = ListGetLogical( Solver % Values, 'Linear System Symmetric', stat)
    matspd = ListGetLogical( Solver % Values, 'Linear System Positive Definite', stat)

    ! force unsymmetric mode when "row equilibration" is used
    scaled = ListGetLogical( Solver % Values, 'Linear System Scaling', stat)
    IF(.NOT.stat) scaled = .TRUE.
    IF(scaled) THEN
      IF(ListGetLogical( Solver % Values, 'Linear System Row Equilibration',stat)) matsym=.FALSE.
    END IF

    IF(matsym) THEN
      IF ( matspd) THEN
        A % MumpsID % sym = 1
      ELSE
        A % MumpsID % sym = 0 ! 2=symmetric, but unsymmetric solver seems faster, at least in a few
                              ! simple cases...  more testing needed...
      END IF
    ELSE
      A % MumpsID % sym = 0
    END IF

    CALL DMumps(A % MumpsID)
    CALL CheckMumpsStatus('Mumps_SolveSystem', 'initialization', 'INFOG', &
        A % MumpsID % infog(1), A % MumpsID % infog(2))

    IF(ASSOCIATED(A % Gorder)) DEALLOCATE(A % Gorder)

    IF(ASSOCIATED(A % ParallelInfo)) THEN
      n = SIZE(A % ParallelInfo % GlobalDOFs)

      ALLOCATE( A % Gorder(n), Owner(n) )
      CALL ContinuousNumbering( A % ParallelInfo, A % Perm, A % Gorder, Owner )

      CALL MPI_ALLREDUCE( SUM(Owner), A % MumpsID % n, &
         1, MPI_INTEGER, MPI_SUM, A % MumpsID % Comm, ierr )
      DEALLOCATE(Owner)
    ELSE
      CALL MPI_ALLREDUCE( A % NumberOfRows, A % MumpsId % n, &
            1, MPI_INTEGER, MPI_MAX, A % Comm, ierr )

      ALLOCATE(A % Gorder(A % NumberOFrows))
      DO i=1,A % NumberOFRows
        A % Gorder(i) = i
      END DO
    END IF

   ! Set matrix for Mumps (unsymmetric case)
    IF (A % mumpsID % sym == 0) THEN
      A % MumpsID % nz_loc = A % Rows(A % NumberOfRows+1)-1

      ALLOCATE( A % MumpsID % irn_loc(A % MumpsID % nz_loc) )
      DO i=1,A % NumberOfRows
        ip = A % Gorder(i)
        DO j=A % Rows(i),A % Rows(i+1)-1
          A % MumpsID % irn_loc(j) = ip
        END DO
      END DO

      ALLOCATE( A % MumpsID % jcn_loc(A % MumpsId % nz_loc) )
      DO i=1,A % MumpsID % nz_loc
        A % MumpsID % jcn_loc(i) = A % Gorder(A % Cols(i))
      END DO
      A % MumpsID % a_loc   => A % values
    ELSE
      ! Set matrix for Mumps (symmetric case)
      nzloc = 0
      DO i=1,A % NumberOfRows
        ! Only output lower triangular part to Mumps
        DO j=A % Rows(i),A % Diag(i)
          nzloc = nzloc + 1
        END DO
      END DO

      A % MumpsID % nz_loc = nzloc

      ALLOCATE( A % MumpsID % irn_loc(A % MumpsID % nz_loc) )
      ALLOCATE( A % MumpsID % jcn_loc(A % MumpsId % nz_loc) )
      ALLOCATE( A % MumpsID % A_loc(A % MumpsId % nz_loc) )

      nzloc = 0
      DO i=1,A % NumberOfRows
        ! Only output lower triangular part to Mumps
        DO j=A % Rows(i),A % Diag(i)
          nzloc = nzloc + 1
          A % mumpsID % IRN_loc(nzloc) = A % Gorder(i)
          A % mumpsID % A_loc(nzloc) = A % Values(j)
          A % mumpsID % JCN_loc(nzloc) = A % Gorder(A % Cols(j))
        END DO
      END DO
    END IF


    ALLOCATE(A % MumpsID % rhs(A % MumpsId % n))

    ! Tune verbosity of MUMPS.
    i = 0
    IF(InfoActive(20)) i = 1
    A % MumpsID % icntl(2) = i ! suppress printing of diagnostics and warnings
    A % MumpsID % icntl(3) = i ! suppress statistics

    A % MumpsID % icntl(4) = 1 ! the same as the two above, but doesn't seem to work.
    A % MumpsID % icntl(5) = 0 ! matrix format 'assembled'

    icntlft = ListGetInteger(Solver % Values, &
          'mumps percentage increase working space', stat)
    IF (stat) THEN
       A % MumpsID % icntl(14) = icntlft
    END IF
    IF (ListGetLogical(Solver % Values, 'mumps null pivot detection', stat)) THEN
       A % MumpsID % icntl(24) = 1
       CALL Info('Mumps_SolveSystem', 'MUMPS null pivot detection enabled', Level=5)
    END IF
    ! CNTL(3) null-pivot threshold: >0 CNTL(3)*||A||inf, <0 |CNTL(3)|, 0 MUMPS default
    ! (4.10.0: 1e-5*eps*||A||inf; 5.4.0 and later: sqrt(pivots on critical path)*eps*||A||inf). Used: DKEEP(1).
    nullpivtol = ListGetConstReal(Solver % Values, 'mumps null pivot tolerance', stat)
    IF (stat) THEN
       A % MumpsID % cntl(3) = nullpivtol
    END IF
    IF (ListGetLogical(Solver % Values, 'mumps sequential root', stat)) THEN
       A % MumpsID % icntl(13) = 1
       CALL Info('Mumps_SolveSystem', 'MUMPS sequential root node (no ScaLAPACK) enabled', Level=5)
    END IF
    A % MumpsID % icntl(18) = 3 ! 'distributed' matrix 
    A % MumpsID % icntl(21) = 1 ! 'distributed' solution phase

    A % MumpsID % job = 4
    CALL DMumps(A % MumpsID)
    CALL Flush(6)
    CALL CheckMumpsStatus('Mumps_SolveSystem', 'analysis and factorization', 'INFOG', &
        A % MumpsID % infog(1), A % MumpsID % infog(2))
    IF (A % MumpsID % icntl(24) == 1) THEN
       CALL Info('Mumps_SolveSystem', 'MUMPS null pivots detected: '//I2S(A % MumpsID % infog(28)), Level=5)
    END IF
    WRITE(Message, '(A,I0,A,I0,A,ES10.3,A,ES10.3)') 'MUMPS ordering used (INFOG(7)): ', A % MumpsID % infog(7), &
        ', entries in factors (INFOG(29)): ', A % MumpsID % infog(29), ', flops (RINFOG(3)): ', A % MumpsID % rinfog(3), &
        ', null-pivot threshold (DKEEP(1)): ', A % MumpsID % dkeep(1)
    CALL Info('Mumps_SolveSystem', Message, Level=5)

    A % MumpsID % lsol_loc = A % mumpsid % info(23)
    ALLOCATE(A % MumpsID % sol_loc(A % MumpsId % lsol_loc))
    ALLOCATE(A % MumpsID % isol_loc(A % MumpsId % lsol_loc))
  END IF

 ! sum the rhs from all procs. Could be done for neighbours only (i guess):
 ! ------------------------------------------------------------------------
  A % MumpsID % RHS = 0
  DO i=1,A % NumberOfRows
    ip = A % Gorder(i)
    A % MumpsId % RHS(ip) = b(i)
  END DO

  ALLOCATE( dbuf(A % MumpsID % n) )
  dbuf = A % MumpsId % RHS
  CALL MPI_ALLREDUCE( dbuf, A % MumpsID % RHS, &
    A % MumpsID % n, MPI_DOUBLE_PRECISION, MPI_SUM, A % MumpsID % Comm, ierr )

 ! Solution:
 ! ---------
  IF (ListGetLogical(Solver % Values, 'Mumps Residual Check', stat)) &
      CALL MumpsResidualCheck(Solver, A, x, b, A % MumpsID % Comm, 'before the solve')
  A % MumpsID % job = 3
  CALL DMumps(A % MumpsID)
  CALL CheckMumpsStatus('Mumps_SolveSystem', 'solve', 'INFOG', &
      A % MumpsID % infog(1), A % MumpsID % infog(2))

 ! Distribute the solution to all:
 ! -------------------------------
  A % MumpsId % Rhs = 0
  DO i=1,A % MumpsID % lsol_loc
    A % MumpsID % RHS(A % MumpsID % isol_loc(i)) = A % MumpsID % sol_loc(i)
  END DO
  dbuf = A % MumpsId % RHS
  CALL MPI_ALLREDUCE( dbuf, A % MumpsID % RHS, &
    A % MumpsID % N, MPI_DOUBLE_PRECISION, MPI_SUM,A %  MumpsID % Comm, ierr )

  DEALLOCATE(dbuf)

 ! Select the values which belong to us:
 ! -------------------------------------
  DO i=1,A % NumberOfRows
    ip = A % Gorder(i)
    x(i) = A % MumpsId % RHS(ip)
  END DO
  CALL MumpsCheckSolution(Solver, A, x, b, A % MumpsID % Comm)

  FreeFactorize = ListGetLogical( Solver % Values, 'Linear System Free Factorization', stat )
  IF ( .NOT. stat ) FreeFactorize = .TRUE.
  IF ( Factorize .AND. FreeFactorize ) CALL FreeMumpsFactorizations(A)

#else
   CALL Fatal( 'Mumps_SolveSystem', 'MUMPS Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE Mumps_SolveSystem
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Solves a complex linear system using MUMPS direct solver. This is a legacy solver
!> with complicated dependencies. This is only available in parallel. 
!------------------------------------------------------------------------------
  SUBROUTINE ZMumps_SolveSystem( Solver,A,x,b )
!------------------------------------------------------------------------------
#ifdef HAVE_MUMPS
#  if defined(ELMER_HAVE_MPI_MODULE)
  USE mpi
#  endif
#endif

  TYPE(Matrix_t) :: A
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp), TARGET :: x(*), b(*)

#ifdef HAVE_MUMPS
#  if defined(ELMER_HAVE_MPIF_HEADER)
  INCLUDE 'mpif.h'
#  endif

  INTEGER, ALLOCATABLE :: Owner(:)
  INTEGER :: i,j,k,l,n,ip,ierr,icntlft,nzloc
  REAL(KIND=dp) :: nullpivtol
  LOGICAL :: Factorize, FreeFactorize, stat, matsym, matspd, scaled

  INTEGER, ALLOCATABLE :: memb(:)
  INTEGER :: Comm_active, Group_active, Group_world

  COMPLEX(KIND=dp), ALLOCATABLE :: dbuf(:)

  Factorize = ListGetLogical( Solver % Values, 'Linear System Refactorize', stat )
  IF ( .NOT. stat ) Factorize = .TRUE.

  IF ( Factorize .OR. .NOT.ASSOCIATED(A % ZMumpsID) ) THEN
    CALL FreeMumpsFactorizations(A)
    ALLOCATE(A % ZMumpsID)

    A % ZMumpsID % Comm = A % Comm
    A % ZMumpsID % par  =  1
    A % ZMumpsID % job  = -1
    A % ZMumpsID % Keep =  0

!   matsym = ListGetLogical( Solver % Values, 'Linear System Symmetric', stat)
!   matspd = ListGetLogical( Solver % Values, 'Linear System Positive Definite', stat)
    matsym = .FALSE.
    matspd = .FALSE.

    ! force unsymmetric mode when "row equilibration" is used
    scaled = ListGetLogical( Solver % Values, 'Linear System Scaling', stat)
    IF(.NOT.stat) scaled = .TRUE.
    IF(scaled) THEN
      IF(ListGetLogical( Solver % Values, 'Linear System Row Equilibration',stat)) matsym=.FALSE.
    END IF

    A % ZMumpsID % sym = 0

!   IF(matsym) THEN
!     IF ( matspd) THEN
!       A % ZMumpsID % sym = 1
!     ELSE
!       A % ZMumpsID % sym = 2 ! 2=symmetric, but unsymmetric solver seems faster, at least in a few
!                             ! simple cases...  more testing needed...
!     END IF
!   ELSE
!     A % ZMumpsID % sym = 0
!   END IF

    CALL ZMumps(A % ZMumpsID)
    CALL CheckMumpsStatus('ZMumps_SolveSystem', 'initialization', 'INFOG', &
        A % ZMumpsID % infog(1), A % ZMumpsID % infog(2))

    IF(ASSOCIATED(A % Gorder)) DEALLOCATE(A % Gorder)

    IF(ASSOCIATED(A % ParallelInfo)) THEN
      n = SIZE(A % ParallelInfo % GlobalDOFs)
      ALLOCATE( A % Gorder(n), Owner(n) )
      CALL ContinuousNumbering( A % ParallelInfo, A % Perm, A % Gorder, Owner )
      CALL MPI_ALLREDUCE( SUM(Owner)/2,A % ZMumpsID % n,1,MPI_INTEGER, MPI_SUM,A % ZMumpsID % Comm,ierr )
      DEALLOCATE(Owner)
    ELSE
      CALL MPI_ALLREDUCE( A % NumberOfRows/2, A % ZMumpsId % n, 1, MPI_INTEGER, MPI_MAX, A % Comm, ierr )

      ALLOCATE(A % Gorder(A % NumberOFrows))
      DO i=1,A % NumberOFRows
        A % Gorder(i) = i
      END DO
    END IF

   ! Set matrix for Mumps (unsymmetric case)
   IF (A % ZmumpsID % sym == 0) THEN ! SYMMETRIC CASE NOT  DONE
      A % ZMumpsID % nz_loc = (A % Rows(A % NumberOfRows+1)-1)/4

      ALLOCATE( A % ZMumpsID % irn_loc(A % ZMumpsID % nz_loc) )
      ALLOCATE( A % ZMumpsID % jcn_loc(A % ZMumpsId % nz_loc) )
      ALLOCATE( A % ZMumpsID % a_loc(A % ZMumpsId % nz_loc) )

      nzloc = 0
      DO i=1,A % NumberOfRows,2
        ip = (A % Gorder(i)-1)/2 + 1
        DO j=A % Rows(i), A % Rows(i+1)-1,2
          nzloc = nzloc + 1
          A % ZMumpsID % irn_loc(nzloc) = ip
          A % ZMumpsID % jcn_loc(nzloc) = (A % Gorder(A % Cols(j))-1)/2+1
          A % ZMumpsID % a_loc(nzloc)   = CMPLX( A % Values(j), -A % Values(j+1), KIND=dp)
        END DO
      END DO
    ELSE
      ! Set matrix for Mumps (symmetric case)
      nzloc = 0
      DO i=1,A % NumberOfRows,2
        ! Only output lower triangular part to Mumps
        DO j=A % Rows(i),A % Diag(i),2
          nzloc = nzloc + 1
        END DO
      END DO

      A % ZMumpsID % nz_loc = nzloc

      ALLOCATE( A % ZMumpsID % irn_loc(A % ZMumpsID % nz_loc) )
      ALLOCATE( A % ZMumpsID % jcn_loc(A % ZMumpsId % nz_loc) )
      ALLOCATE( A % ZMumpsID % A_loc(A % ZMumpsId % nz_loc) )

      nzloc = 0
      DO i=1,A % NumberOfRows,2
        ! Only output lower triangular part to Mumps
        DO j=A % Rows(i),A % Diag(i),2
          nzloc = nzloc + 1
          A % ZmumpsID % IRN_loc(nzloc) = (A % Gorder(i)-1)/2+1
          A % ZmumpsID % JCN_loc(nzloc) = (A % Gorder(A % Cols(j))-1)/2+1
          A % ZmumpsID % A_loc(nzloc) = CMPLX( A % Values(j), -A % Values(j+1), KIND=dp)
        END DO
      END DO
    END IF

    ALLOCATE(A % ZMumpsID % rhs(A % ZMumpsId % n))

    ! Tune verbosity of MUMPS.
    i = 0
    IF(InfoActive(20)) i = 1
    A % ZMumpsID % icntl(2) = i ! suppress printing of diagnostics and warnings
    A % ZMumpsID % icntl(3) = i ! suppress statistics

    A % ZMumpsID % icntl(4) = 1 ! the same as the two above, but doesn't seem to work.
    A % ZMumpsID % icntl(5) = 0 ! matrix format 'assembled'

    icntlft = ListGetInteger(Solver % Values, 'mumps percentage increase working space', stat)
    IF (stat) THEN
       A % ZMumpsID % icntl(14) = icntlft
    END IF
    IF (ListGetLogical(Solver % Values, 'mumps null pivot detection', stat)) THEN
       A % ZMumpsID % icntl(24) = 1
       CALL Info('ZMumps_SolveSystem', 'MUMPS null pivot detection enabled', Level=5)
    END IF
    ! CNTL(3) null-pivot threshold: >0 CNTL(3)*||A||inf, <0 |CNTL(3)|, 0 MUMPS default
    ! (4.10.0: 1e-5*eps*||A||inf; 5.4.0 and later: sqrt(pivots on critical path)*eps*||A||inf). Used: DKEEP(1).
    nullpivtol = ListGetConstReal(Solver % Values, 'mumps null pivot tolerance', stat)
    IF (stat) THEN
       A % ZMumpsID % cntl(3) = nullpivtol
    END IF
    IF (ListGetLogical(Solver % Values, 'mumps sequential root', stat)) THEN
       A % ZMumpsID % icntl(13) = 1
       CALL Info('ZMumps_SolveSystem', 'MUMPS sequential root node (no ScaLAPACK) enabled', Level=5)
    END IF
    A % ZMumpsID % icntl(18) = 3 ! 'distributed' matrix 
    A % ZMumpsID % icntl(21) = 1 ! 'distributed' solution phase

    A % ZMumpsID % job = 4
    CALL ZMumps(A % ZMumpsID)
    CALL Flush(6)
    CALL CheckMumpsStatus('ZMumps_SolveSystem', 'analysis and factorization', 'INFOG', &
        A % ZMumpsID % infog(1), A % ZMumpsID % infog(2))
    IF (A % ZMumpsID % icntl(24) == 1) THEN
       CALL Info('ZMumps_SolveSystem', 'MUMPS null pivots detected: '//I2S(A % ZMumpsID % infog(28)), Level=5)
    END IF
    WRITE(Message, '(A,I0,A,I0,A,ES10.3,A,ES10.3)') 'MUMPS ordering used (INFOG(7)): ', A % ZMumpsID % infog(7), &
        ', entries in factors (INFOG(29)): ', A % ZMumpsID % infog(29), ', flops (RINFOG(3)): ', A % ZMumpsID % rinfog(3), &
        ', null-pivot threshold (DKEEP(1)): ', A % ZMumpsID % dkeep(1)
    CALL Info('ZMumps_SolveSystem', Message, Level=5)

    A % ZMumpsID % lsol_loc = A % Zmumpsid % info(23)
    ALLOCATE(A % ZMumpsID % sol_loc(A % ZMumpsId % lsol_loc))
    ALLOCATE(A % ZMumpsID % isol_loc(A % ZMumpsId % lsol_loc))
  END IF

 ! sum the rhs from all procs. Could be done
 ! for neighbours only (i guess):
 ! ------------------------------------------
  A % ZMumpsID % RHS = 0
  DO i=1,A % NumberOfRows,2
    ip = (A % Gorder(i)-1)/2+1
    A % ZMumpsId % RHS(ip) = CMPLX(b(i), b(i+1), KIND=dp)
  END DO
  ALLOCATE( dbuf(A % ZMumpsID % n) )
  dbuf = A % ZMumpsId % RHS
  CALL MPI_ALLREDUCE( dbuf, A % ZMumpsID % RHS, &
    A % ZMumpsID % n, MPI_DOUBLE_COMPLEX, MPI_SUM, A % ZMumpsID % Comm, ierr )

 ! Solution:
 ! ---------
  IF (ListGetLogical(Solver % Values, 'Mumps Residual Check', stat)) &
      CALL MumpsResidualCheck(Solver, A, x, b, A % ZMumpsID % Comm, 'before the solve')
  A % ZMumpsID % job = 3
  CALL ZMumps(A % ZMumpsID)
  CALL CheckMumpsStatus('ZMumps_SolveSystem', 'solve', 'INFOG', &
      A % ZMumpsID % infog(1), A % ZMumpsID % infog(2))

 ! Distribute the solution to all:
 ! -------------------------------
  A % ZMumpsId % Rhs = 0
  DO i=1,A % ZMumpsID % lsol_loc
    A % ZMumpsID % RHS(A % ZMumpsID % isol_loc(i)) = A % ZMumpsID % sol_loc(i)
  END DO
  dbuf = A % ZMumpsId % RHS
  CALL MPI_ALLREDUCE( dbuf, A % ZMumpsID % RHS, &
    A % ZMumpsID % N, MPI_DOUBLE_COMPLEX, MPI_SUM,A %  ZMumpsID % Comm, ierr )

  DEALLOCATE(dbuf)

 ! Select the values which belong to us:
 ! -------------------------------------
  DO i=1,A % NumberOfRows,2
    ip = (A % Gorder(i)-1)/2+1
    x(i)   = REAL( A % ZMumpsId % RHS(ip) )
    x(i+1) = AIMAG( A % ZMumpsId % RHS(ip) )
  END DO
  CALL MumpsCheckSolution(Solver, A, x, b, A % ZMumpsID % Comm)

  FreeFactorize = ListGetLogical( Solver % Values, 'Linear System Free Factorization', stat )
  IF ( .NOT. stat ) FreeFactorize = .TRUE.
  IF ( Factorize .AND. FreeFactorize ) CALL FreeMumpsFactorizations(A)
#else
   CALL Fatal( 'ZMumps_SolveSystem', 'MUMPS Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE ZMumps_SolveSystem
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solves local linear system using MUMPS direct solver. If the solved system
!> is singular, optionally one possible solution is returned.
!------------------------------------------------------------------------------
  SUBROUTINE MumpsLocal_SolveSystem( Solver, A, x, b, Free_Fact )
!------------------------------------------------------------------------------
     IMPLICIT NONE

     TYPE(Matrix_t) :: A
     TYPE(Solver_t) :: Solver
     REAL(KIND=dp), TARGET :: x(*), b(*)
     LOGICAL, OPTIONAL :: Free_Fact

     INTEGER :: i
     LOGICAL :: Factorize, FreeFactorize, stat

#ifdef HAVE_MUMPS
     ! Free local Mumps instance if requested
     IF (PRESENT(Free_Fact)) THEN
       IF (Free_Fact) THEN
         CALL MumpsLocal_Free(A)
         RETURN
       END IF
     END IF

     ! Refactorize local matrix if needed
     Factorize = ListGetLogical( Solver % Values, &
         'Linear System Refactorize', stat )
     IF (.NOT. stat) Factorize = .TRUE.

     IF (Factorize .OR. .NOT. ASSOCIATED(A % mumpsIDL)) THEN
       CALL MumpsLocal_Factorize(Solver, A)
     END IF
       
     ! Set RHS
     A % mumpsIDL % NRHS = 1
     A % mumpsIDL % LRHS = A % mumpsIDL % n
     DO i=1,A % NumberOfRows
       A % mumpsIDL % RHS(i) = b(i)
     END DO
     ! We could use BLAS here..
     ! CALL DCOPY(A % NumberOfRows, b, 1, A % mumpsIDL % RHS, 1)

     ! SOLUTION PHASE
     A % mumpsIDL % job = 3
     CALL DMumps(A % mumpsIDL)
     CALL CheckMumpsStatus('MumpsLocal_SolveSystem', 'solve', 'INFO', &
         A % mumpsIDL % info(1), A % mumpsIDL % info(2))

     ! TODO: If solution is not local, redistribute the solution vector here

     ! Set local solution
     DO i=1,A % NumberOfRows
       x(i) = A % mumpsIDL % RHS(i)
     END DO
     ! We could use BLAS here..
     ! CALL DCOPY(A % NumberOfRows, A % mumpsIDL % RHS, 1, x, 1)

     FreeFactorize = ListGetLogical( Solver % Values, &
         'Linear System Free Factorization', stat )
     IF (.NOT. stat) FreeFactorize = .TRUE.

     IF (Factorize .AND. FreeFactorize) CALL MumpsLocal_Free(A)
#else
     CALL Fatal( 'MumpsLocal_SolveSystem', 'MUMPS Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE MumpsLocal_SolveSystem
!------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
!> Solves local linear system using MUMPS direct solver. If the solved system
!> is singular, optionally one possible solution is returned.
!------------------------------------------------------------------------------
  SUBROUTINE ZMumpsLocal_SolveSystem( Solver, A, x, b, Free_Fact )
!------------------------------------------------------------------------------
     IMPLICIT NONE

     TYPE(Matrix_t) :: A
     TYPE(Solver_t) :: Solver
     REAL(KIND=dp), TARGET :: x(*), b(*)
     LOGICAL, OPTIONAL :: Free_Fact

     INTEGER :: i,j
     LOGICAL :: Factorize, FreeFactorize, stat

#ifdef HAVE_MUMPS
     ! Free local Mumps instance if requested
     IF (PRESENT(Free_Fact)) THEN
       IF (Free_Fact) THEN
         CALL MumpsLocal_Free(A)
         RETURN
       END IF
     END IF

     ! Refactorize local matrix if needed
     Factorize = ListGetLogical( Solver % Values,'Linear System Refactorize', stat )
     IF (.NOT. stat) Factorize = .TRUE.

     IF (Factorize .OR. .NOT. ASSOCIATED(A % ZmumpsIDL)) THEN
       CALL ZMumpsLocal_Factorize(Solver, A)
     END IF

     ! Set RHS
     A % ZmumpsIDL % NRHS = 1
     A % ZmumpsIDL % LRHS = A % ZmumpsIDL % n
     j = 0
     DO i=1,A % NumberOfRows,2
       j = j + 1
       A % ZmumpsIDL % RHS(j) = CMPLX(b(i), b(i+1), KIND=dp)
     END DO
     ! We could use BLAS here..
     ! CALL DCOPY(A % NumberOfRows, b, 1, A % mumpsIDL % RHS, 1)

     ! SOLUTION PHASE
     A % ZmumpsIDL % job = 3
     CALL ZMumps(A % ZmumpsIDL)
     CALL CheckMumpsStatus('ZMumpsLocal_SolveSystem', 'solve', 'INFO', &
         A % ZmumpsIDL % info(1), A % ZmumpsIDL % info(2))

     ! TODO: If solution is not local, redistribute the solution vector here

     ! Set local solution
     j = 0
     DO i=1,A % NumberOfRows,2
       j = j + 1
       x(i) = REAL(A % ZmumpsIDL % RHS(j) )
       x(i+1) = AIMAG(A % ZmumpsIDL % RHS(j) )
     END DO
     ! We could use BLAS here..
     ! CALL DCOPY(A % NumberOfRows, A % mumpsIDL % RHS, 1, x, 1)

     FreeFactorize = ListGetLogical( Solver % Values, 'Linear System Free Factorization', stat )
     IF (.NOT. stat) FreeFactorize = .TRUE.

     IF (Factorize .AND. FreeFactorize) CALL MumpsLocal_Free(A)
#else
     CALL Fatal( 'ZMumpsLocal_SolveSystem', 'MUMPS Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
   END SUBROUTINE ZMumpsLocal_SolveSystem
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Factorize local matrix with Mumps
!------------------------------------------------------------------------------
  SUBROUTINE MumpsLocal_Factorize(Solver, A)
!------------------------------------------------------------------------------
#ifdef HAVE_MUMPS
#  if defined(ELMER_HAVE_MPI_MODULE)
    USE mpi
#  endif
#endif
    IMPLICIT NONE

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A

#ifdef HAVE_MUMPS
#  if defined(ELMER_HAVE_MPIF_HEADER)
    INCLUDE 'mpif.h'
#  endif

    INTEGER :: i, j, n, nz, allocstat, icntlft, ptype, nzloc
    LOGICAL :: matspd, matsym, nullpiv, stat

    ! INTEGER :: myrank, ierr
    ! CHARACTER(len=32) :: buf

    IF ( ASSOCIATED(A % mumpsIDL) ) CALL MumpsLocal_Free(A)
    ALLOCATE(A % mumpsIDL)

    ! INITIALIZATION PHASE

    ! Initialize local instance of Mumps
    ! TODO (Hybridization): change this if local system needs to be solved
    ! with several cores
    A % mumpsIDL % COMM = MPI_COMM_SELF
    A % mumpsIDL % PAR = 1 ! Host (=self) takes part in factorization

    ! Check if matrix is symmetric or spd
    matsym = ListGetLogical(Solver % Values,'Linear System Symmetric',stat)
    matspd = ListGetLogical(Solver % Values, &
                                   'Linear System Positive Definite', stat)

    A % mumpsIDL % SYM = 0
    IF (matsym) THEN
      IF (matspd) THEN
        ! Matrix is symmetric positive definite
        A % mumpsIDL % SYM = 1
      ELSE
        ! Matrix is symmetric
        A % mumpsIDL % SYM = 2
     END IF
    ELSE
      ! Matrix is unsymmetric
      A % mumpsIDL % SYM = 0
    END IF
    A % mumpsIDL % JOB  = -1 ! Initialize
    CALL DMumps(A % mumpsIDL)
    CALL CheckMumpsStatus('MumpsLocal_Factorize', 'initialization', 'INFO', &
        A % mumpsIDL % INFO(1), A % mumpsIDL % INFO(2))

    ! FACTORIZE PHASE

    ! Set stdio parameters

    A % mumpsIDL % ICNTL(1)  =  6  ! Error messages to stdout
    A % mumpsIDL % ICNTL(2)  = -1 ! No diagnostic and warning messages
    A % mumpsIDL % ICNTL(3)  = -1 ! No statistic messages
    A % mumpsIDL % ICNTL(4)  =  1  ! Print only error messages

    ! Set matrix format
    A % mumpsIDL % ICNTL(5)  = 0  ! Assembled matrix format
    A % mumpsIDL % ICNTL(18) = 0 ! Centralized matrix
    A % mumpsIDL % ICNTL(21) = 0 ! Centralized dense solution phase

    ! Check if solution of singular systems is ok
    A % mumpsIDL % ICNTL(24) = 0
    nullpiv = ListGetLogical(Solver % Values, 'Mumps Solve Singular', stat)
    IF (nullpiv) THEN
      A % mumpsIDL % ICNTL(24) = 1
      A % mumpsIDL % CNTL(1) = 1D-2     ! Pivoting threshold
      A % mumpsIDL % CNTL(3) = 1D-9    ! Null pivot detection threshold
      A % mumpsIDL % CNTL(5) = 1D6      ! Fixation value for null pivots
      A % mumpsIDL % CNTL(13) = 1       ! Do not use ScaLAPACK on the root node
      ! TODO: if needed, here set CNTL(3) and CNTL(5) as parameters for
      ! more accurate null pivot detection
    END IF

    ! Set permutation strategy for Mumps
    ptype = ListGetInteger(Solver % Values, 'Mumps Permutation Type', stat)
    IF (stat) THEN
      A % mumpsIDL % ICNTL(6) = ptype
    END IF

    ! TODO: Change this if system is larger than local.
    ! For larger than local systems define global->local numbering
    n = A % NumberofRows
    nz = A % Rows(A % NumberOfRows+1)-1
    A % mumpsIDL % N  = n
    A % mumpsIDL % NZ = nz
    ! A % mumpsIDL % nz_loc = nz

    ! Allocate rows and columns for MUMPS
    ALLOCATE( A % mumpsIDL % IRN(nz), &
          A % mumpsIDL % JCN(nz), &
          A % mumpsIDL % A(nz), STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('MumpsLocal_Factorize', &
            'Memory allocation for MUMPS row and column indices failed.')
    END IF

    ! Set matrix for Mumps (unsymmetric case)
    IF (A % mumpsIDL % sym == 0) THEN
      DO i=1,A % NumberOfRows
         DO j=A % Rows(i),A % Rows(i+1)-1
            A % mumpsIDL % IRN(j) = i
         END DO
       END DO

       ! Set columns and values
       DO i=1,A % mumpsIDL % nz
         A % mumpsIDL % JCN(i) = A % Cols(i)
       END DO
       DO i=1,A % mumpsIDL % nz
         A % mumpsIDL % A(i) = A % Values(i)
       END DO
    ELSE
      ! Set matrix for Mumps (symmetric case)
      nzloc = 0
      DO i=1,n
        DO j=A % Rows(i),A % Rows(i+1)-1
          ! Only output lower triangular part to Mumps
          IF (i<=A % Cols(j)) THEN
            nzloc = nzloc + 1
            A % mumpsIDL % IRN(nzloc) = i
            A % mumpsIDL % JCN(nzloc) = A % Cols(j)
            A % mumpsIDL % A(nzloc) = A % Values(j)
          END IF
        END DO
      END DO
    END IF

    icntlft = ListGetInteger(Solver % Values, 'mumps percentage increase working space', stat)
    IF (stat) THEN
       A % mumpsIDL % ICNTL(14) = icntlft
    END IF

    A % mumpsIDL % JOB = 1 ! Perform analysis
    CALL DMumps(A % mumpsIDL)
    CALL Flush(6)

    ! Check return status
    IF (A % mumpsIDL % INFO(1)<0) THEN
      CALL Fatal('MumpsLocal_Factorize','Mumps analysis phase failed')
    END IF

    A % mumpsIDL % JOB = 2 ! Perform factorization
    CALL DMumps(A % mumpsIDL)
    CALL Flush(6)

    ! Check return status
    IF (A % mumpsIDL % INFO(1)<0) THEN
      CALL Fatal('MumpsLocal_Factorize','Mumps factorize phase failed')
    END IF

    ! Allocate RHS
    ALLOCATE(A % mumpsIDL % RHS(A % mumpsIDL % N), STAT=allocstat)
    IF (allocstat /= 0) THEN
         CALL Fatal('MumpsLocal_Factorize', &
                 'Memory allocation for MUMPS solution vector and RHS failed.' )
    END IF
#else
   CALL Fatal( 'MumpsLocal_Factorize', 'MUMPS Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE MumpsLocal_Factorize
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Factorize local matrix with Mumps
!------------------------------------------------------------------------------
  SUBROUTINE ZMumpsLocal_Factorize(Solver, A)
!------------------------------------------------------------------------------
#ifdef HAVE_MUMPS
#  if defined(ELMER_HAVE_MPI_MODULE)
    USE mpi
#  endif
#endif
    IMPLICIT NONE

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A

#ifdef HAVE_MUMPS
#  if defined(ELMER_HAVE_MPIF_HEADER)
    INCLUDE 'mpif.h'
#  endif

    INTEGER :: i, j, k, n, nz, allocstat, icntlft, ptype, nzloc
    LOGICAL :: matspd, matsym, nullpiv, stat

    ! INTEGER :: myrank, ierr
    ! CHARACTER(len=32) :: buf

    IF ( ASSOCIATED(A % ZmumpsIDL) ) CALL MumpsLocal_Free(A)

    ! INITIALIZATION PHASE
    ALLOCATE(A % ZmumpsIDL)

    ! Initialize local instance of Mumps
    ! TODO (Hybridization): change this if local system needs to be solved
    ! with several cores
    A % ZmumpsIDL % COMM = MPI_COMM_SELF
    A % ZmumpsIDL % PAR = 1 ! Host (=self) takes part in factorization

    ! Check if matrix is symmetric or spd
    matsym = ListGetLogical(Solver % Values, &
        'Linear System Symmetric', stat)
 
    matspd = ListGetLogical(Solver % Values, &
        'Linear System Positive Definite', stat)
 
    A % ZmumpsIDL % SYM = 0
!   IF (matsym) THEN
!     IF (matspd) THEN
!       ! Matrix is symmetric positive definite
!       A % ZmumpsIDL % SYM = 1
!     ELSE
!       ! Matrix is symmetric
!       A % ZmumpsIDL % SYM = 2
!    END IF
!   ELSE
!     ! Matrix is unsymmetric
!     A % ZmumpsIDL % SYM = 0
!   END IF

    A % ZmumpsIDL % JOB  = -1 ! Initialize
    CALL ZMumps(A % ZmumpsIDL)
    CALL CheckMumpsStatus('ZMumpsLocal_Factorize', 'initialization', 'INFO', &
        A % ZmumpsIDL % INFO(1), A % ZmumpsIDL % INFO(2))

    ! FACTORIZE PHASE

    ! Set stdio parameters
    A % ZmumpsIDL % ICNTL(1)  = 6  ! Error messages to stdout
    A % ZmumpsIDL % ICNTL(2)  = -1 ! No diagnostic and warning messages
    A % ZmumpsIDL % ICNTL(3)  = -1 ! No statistic messages
    A % ZmumpsIDL % ICNTL(4)  = 1  ! Print only error messages

    ! Set matrix format
    A % ZmumpsIDL % ICNTL(5)  = 0  ! Assembled matrix format
    A % ZmumpsIDL % ICNTL(18) = 0 ! Centralized matrix
    A % ZmumpsIDL % ICNTL(21) = 0 ! Centralized dense solution phase

    ! Check if solution of singular systems is ok
    A % ZmumpsIDL % ICNTL(24) = 0
    nullpiv = ListGetLogical(Solver % Values, 'Mumps Solve Singular', stat)
    IF (nullpiv) THEN
      A % ZmumpsIDL % ICNTL(24) = 1
      A % ZmumpsIDL % CNTL(1) = 1D-2     ! Pivoting threshold
      A % ZmumpsIDL % CNTL(3) = 1D-9     ! Null pivot detection threshold
      A % ZmumpsIDL % CNTL(5) = 1D6      ! Fixation value for null pivots
      A % ZmumpsIDL % CNTL(13) = 1       ! Do not use ScaLAPACK on the root node
      ! TODO: if needed, here set CNTL(3) and CNTL(5) as parameters for
      ! more accurate null pivot detection
    END IF

    ! Set permutation strategy for Mumps
    ptype = ListGetInteger(Solver % Values, 'Mumps Permutation Type', stat)
    IF (stat) A % ZmumpsIDL % ICNTL(6) = ptype

    ! TODO: Change this if system is larger than local.
    ! For larger than local systems define global->local numbering
    n = A % NumberofRows / 2
    nz = (A % Rows(A % NumberOfRows+1)-1) / 4
    A % ZmumpsIDL % N  = n
    A % ZmumpsIDL % nz = nz
    ! A % mumpsIDL % nz_loc = nz

    ! Allocate rows and columns for MUMPS
    ALLOCATE( A % ZmumpsIDL % IRN(nz), &
              A % ZmumpsIDL % JCN(nz), &
              A % ZmumpsIDL % A(nz), STAT=allocstat)
    IF (allocstat /= 0) THEN
      CALL Fatal('MumpsLocal_Factorize', &
            'Memory allocation for MUMPS row and column indices failed.')
    END IF

    ! Set matrix for Mumps (unsymmetric case)
    IF (A % ZmumpsIDL % sym == 0) THEN
      nzloc = 0
      DO i=1,A % NumberOfRows,2
         DO j=A % Rows(i),A % Rows(i+1)-1,2
            nzloc = nzloc + 1
            A % ZmumpsIDL % IRN(nzloc) = i/2+1
         END DO
       END DO

       ! Set columns and values
       nzloc = 0
       DO i=1,A % NumberOfRows,2
         DO j=A % Rows(i),A % Rows(i+1)-1,2
           nzloc = nzloc + 1
           A % ZmumpsIDL % JCN(nzloc) = A % Cols(j)/2+1
           A % ZmumpsIDL % A(nzloc) = CMPLX(A % Values(j), -A % Values(j+1), KIND=dp )
         END DO
       END DO
    ELSE
      ! Set matrix for Mumps (symmetric case)
      nzloc = 0
      DO i=1,A % NumberOfRows,2
        DO j=A % Rows(i),A % Rows(i+1)-1,2
          ! Only output lower triangular part to Mumps
          IF (i<=A % Cols(j)) THEN
            nzloc = nzloc + 1
            A % ZmumpsIDL % IRN(nzloc) = i/2+1
            A % ZmumpsIDL % JCN(nzloc) = A % Cols(j)/2+1
            A % ZmumpsIDL % A(nzloc) = CMPLX(A % Values(j), -A % Values(j+1),KIND=dp)
          END IF
        END DO
      END DO
    END IF

    icntlft = ListGetInteger(Solver % Values, 'mumps percentage increase working space', stat)
    IF (stat) THEN
       A % ZmumpsIDL % ICNTL(14) = icntlft
    END IF

    A % ZmumpsIDL % JOB = 1 ! Perform analysis
    CALL ZMumps(A % ZmumpsIDL)
    CALL Flush(6)

    ! Check return status
    IF (A % ZmumpsIDL % INFO(1)<0) THEN
      CALL Fatal('MumpsLocal_Factorize','Mumps analysis phase failed')
    END IF

    A % ZmumpsIDL % JOB = 2 ! Perform factorization
    CALL ZMumps(A % ZmumpsIDL)
    CALL Flush(6)

    ! Check return status
    IF (A % ZmumpsIDL % INFO(1)<0) THEN
      CALL Fatal('ZMumpsLocal_Factorize','Mumps factorize phase failed')
    END IF

    ! Allocate RHS
    ALLOCATE(A % ZmumpsIDL % RHS(A % ZmumpsIDL % N), STAT=allocstat)
    IF (allocstat /= 0) THEN
         CALL Fatal('ZMumpsLocal_Factorize', &
                 'Memory allocation for MUMPS solution vector and RHS failed.' )
    END IF
#else
   CALL Fatal( 'ZMumpsLocal_Factorize', 'MUMPS Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
 END SUBROUTINE ZMumpsLocal_Factorize
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Solve local nullspace using MUMPS direct solver. On exit, z will be
!> allocated and will hold the jth local nullspace vectors as z(j,:).
!------------------------------------------------------------------------------
  SUBROUTINE MumpsLocal_SolveNullSpace(Solver, A, z, nz)
!------------------------------------------------------------------------------
#ifdef HAVE_MUMPS
#  if defined(ELMER_HAVE_MPI_MODULE)
      USE mpi
#  endif
#endif
      IMPLICIT NONE

      TYPE(Solver_t) :: Solver
      TYPE(Matrix_t) :: A
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:), TARGET :: z
      INTEGER :: nz, nrhs

#ifdef HAVE_MUMPS
#  if defined(ELMER_HAVE_MPIF_HEADER)
      INCLUDE 'mpif.h'
#  endif

      INTEGER :: j,k,n, allocstat
      LOGICAL :: Factorize, FreeFactorize, stat

      REAL(KIND=dp), DIMENSION(:), POINTER :: rhs_tmp, nspace

      Factorize = ListGetLogical( Solver % Values,&
                                  'Linear System Refactorize', stat )
      IF (.NOT. stat) Factorize = .TRUE.

      ! Refactorize local matrix if needed
      IF ( Factorize .OR. (.NOT. ASSOCIATED(A % mumpsIDL)) ) THEN
            CALL MumpsLocal_Factorize(Solver, A)
      END IF

      ! Check symmetry condition
      IF (.NOT. (A % mumpsIDL % sym == 1 .OR. A % mumpsIDL % sym == 2)) THEN
         CALL Fatal( 'MumpsLocal_SolveNullSpace', &
           'Mumps null space computation not supported for unsymmetric matrices.')
      END IF

      ! Get the size of the local nullspace and allocate memory
      ! TODO: this is incorrect, should be number of columns
      n = A % NumberOfRows
      nz = A % mumpsIDL % INFOG(28)

      IF (nz == 0) THEN
            ALLOCATE(z(nz,n))
          RETURN
      END IF


      print*,parenv % mype, n,nz; flush(6)
      ALLOCATE(nspace(n*nz), STAT=allocstat)
      IF (allocstat /= 0) THEN
         CALL Fatal( 'MumpsLocal_SolveNullSpace', &
               'Storage allocation for local nullspace failed.')
      END IF

      rhs_tmp => A % mumpsIDL % RHS
      nrhs = A % mumpsIDL % NRHS
      A % mumpsIDL % RHS => nspace
      A % mumpsIDL % NRHS = nz

      ! Solve for the complete local nullspace
      A % mumpsIDL % JOB = 3
      A % mumpsIDL % ICNTL(25) = -1
      CALL DMumps(A % mumpsIDL)
      ! Check return status
      IF (A % mumpsIDL % INFO(1)<0) THEN
          CALL Fatal('MumpsLocal_SolveNullSpace','Mumps nullspace solution failed')
      END IF

      ! Restore pointer to local solution vector
      A % mumpsIDL % ICNTL(25) = 0
      A % mumpsIDL % RHS => rhs_tmp
      A % mumpsIDL % NRHS = nrhs

      FreeFactorize = ListGetLogical( Solver % Values, &
            'Linear System Free Factorization', stat )
      IF (.NOT. stat) FreeFactorize = .TRUE.

      IF (Factorize .AND. FreeFactorize) THEN
          CALL MumpsLocal_Free(A)
       END IF

      ALLOCATE(z(nz,n), STAT=allocstat)
      IF (allocstat /= 0) THEN
          CALL Fatal( 'MumpsLocal_SolveNullSpace', &
                           'Storage allocation for local nullspace failed.')
      END IF

      ! Copy computed nullspace to z and deallocate nspace
      DO j=1,nz
            ! z is row major, cannot use DCOPY here
        DO k=1,n
              z(j,k)=nspace(n*(j-1)+k)
        END DO
      END DO
      DEALLOCATE(nspace)
#else
   CALL Fatal( 'MumpsLocal_SolveNullSpace', 'MUMPS Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE MumpsLocal_SolveNullSpace
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Free local Mumps variables and solver internal allocations
!------------------------------------------------------------------------------
  SUBROUTINE MumpsLocal_Free(A)
!------------------------------------------------------------------------------
        IMPLICIT NONE

        TYPE(Matrix_t) :: A

#ifdef HAVE_MUMPS
        IF (ASSOCIATED(A % mumpsIDL)) THEN
          ! Deallocate Mumps structures
          DEALLOCATE( A % mumpsIDL % irn, A % mumpsIDL % jcn, &
              A % mumpsIDL % a, A % mumpsIDL % rhs)

          ! Free Mumps internal allocations
          A % mumpsIDL % job = -2
          CALL DMumps(A % mumpsIDL)
          DEALLOCATE(A % mumpsIDL)

          A % mumpsIDL => NULL()
        END IF
        IF (ASSOCIATED(A % ZmumpsIDL)) THEN
          ! Deallocate Mumps structures
          DEALLOCATE( A % ZmumpsIDL % irn, A % ZmumpsIDL % jcn, &
              A % ZmumpsIDL % a, A % ZmumpsIDL % rhs)

          ! Free Mumps internal allocations
          A % ZmumpsIDL % job = -2
          CALL ZMumps(A % ZmumpsIDL)
          DEALLOCATE(A % ZmumpsIDL)

          A % ZmumpsIDL => NULL()
        END IF
#else
     CALL Fatal( 'MumpsLocal_Free', 'MUMPS Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE MumpsLocal_Free
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Solves a linear system using SuperLU direct solver.
!------------------------------------------------------------------------------
  SUBROUTINE SuperLU_SolveSystem( Solver,A,x,b,Free_Fact )
!------------------------------------------------------------------------------
  LOGICAL, OPTIONAL :: Free_fact
  TYPE(Matrix_t) :: A
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp), TARGET :: x(*), b(*)

#ifdef HAVE_SUPERLU
      LOGICAL :: stat, Factorize, FreeFactorize
      integer :: n, nnz, nrhs, iinfo, iopt, nprocs

      interface
        subroutine solve_superlu( iopt, nprocs, n, nnz, nrhs, values, cols, &
                     rows, b, ldb, factors, iinfo )
            use types
            integer :: iopt, nprocs, n, nnz, nrhs, cols(*), rows(*), ldb, iinfo
            real(kind=dp) :: values(*), b(*)
            integer(kind=addrint) :: factors
        end subroutine solve_superlu
      end interface

      IF ( PRESENT(Free_Fact) ) THEN
        IF ( Free_Fact ) THEN
          IF ( A % SuperLU_Factors/= 0 ) THEN
            iopt = 3
            CALL Solve_SuperLU( iopt, nprocs, n, nnz, nrhs, A % Values, &
               A % Cols, A % Rows, x, n, A % SuperLU_Factors, iinfo )
            A % SuperLU_Factors = 0
          END IF
          RETURN
        END IF
      END IF

      n = A % NumberOfRows
      nrhs = 1
      x(1:n) = b(1:n)
      nnz = A % Rows(n+1)-1

      nprocs = ListGetInteger( Solver % Values, & 
              'Linear System Number of Threads', stat )
      IF ( .NOT. stat ) nprocs = 1
!
      Factorize = ListGetLogical( Solver % Values, &
         'Linear System Refactorize', stat )
      IF ( .NOT. stat ) Factorize = .TRUE.

      IF ( Factorize .OR. A % SuperLU_Factors==0 ) THEN

        IF ( A % SuperLU_Factors/= 0 ) THEN
          iopt = 3
          call Solve_SuperLU( iopt, nprocs, n, nnz, nrhs, A % Values, A % Cols, &
                     A % Rows, x, n, A % SuperLU_Factors, iinfo )
          A % SuperLU_Factors=0
        END IF

        ! First, factorize the matrix. The factors are stored in *factors* handle.
        iopt = 1
        call Solve_SuperLU( iopt, nprocs, n, nnz, nrhs, A % Values, A % Cols, &
              A % Rows, x, n, A % SuperLU_Factors, iinfo )
 
        if (iinfo .eq. 0) then
           write (*,*) 'Factorization succeeded'
        else
           write(*,*) 'INFO from factorization = ', iinfo
        endif
      END IF

      !
      ! Second, solve the system using the existing factors.
      iopt = 2
      call Solve_SuperLU( iopt, nprocs, n, nnz, nrhs, A % Values, A % Cols, &
           A % Rows,  x, n, A % SuperLU_Factors, iinfo )
!
      if (iinfo .eq. 0) then
         write (*,*) 'Solve succeeded'
      else
         write(*,*) 'INFO from triangular solve = ', iinfo
      endif

! Last, free the storage allocated inside SuperLU
      FreeFactorize = ListGetLogical( Solver % Values, &
          'Linear System Free Factorization', stat )
      IF ( .NOT. stat ) FreeFactorize = .TRUE.

      IF ( Factorize .AND. FreeFactorize ) THEN
        iopt = 3
        call Solve_SuperLU( iopt, nprocs, n, nnz, nrhs, A % Values, A % Cols, &
                    A % Rows, x, n, A % SuperLU_Factors, iinfo )
        A % SuperLU_Factors = 0
      END IF
#endif      
!------------------------------------------------------------------------------
  END SUBROUTINE SuperLU_SolveSystem
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Permon solver
!------------------------------------------------------------------------------
  SUBROUTINE Permon_SolveSystem( Solver,A,x,b,Free_Fact )
!------------------------------------------------------------------------------
#ifdef HAVE_FETI4I
  USE feti4i
#  if defined(ELMER_HAVE_MPI_MODULE)
  USE mpi
#  endif
#endif

  LOGICAL, OPTIONAL :: Free_Fact
  TYPE(Matrix_t) :: A
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp), TARGET :: x(*), b(*)

#ifdef HAVE_FETI4I
#  if defined(ELMER_HAVE_MPIF_HEADER)
  INCLUDE 'mpif.h'
#  endif

  INTEGER, ALLOCATABLE :: Owner(:)
  INTEGER :: i,j,n,nd,ip,ierr,icntlft,nzloc
  LOGICAL :: Factorize, FreeFactorize, stat, matsym, matspd, scaled

  INTEGER :: n_dof_partition

  INTEGER, ALLOCATABLE :: memb(:), DirichletInds(:), Neighbours(:), IntOptions(:)
  REAL(KIND=dp), ALLOCATABLE :: DirichletVals(:), RealOptions(:)
  INTEGER :: Comm_active, Group_active, Group_world

  REAL(KIND=dp), ALLOCATABLE :: dbuf(:)


  n_dof_partition = A % NumberOfRows

!  INTERFACE
!     FUNCTION Permon_InitSolve(n, gnum, nd, dinds, dvals, n_n, n_ranks) RESULT(handle) BIND(c,name='permon_initsolve') 
!        USE, INTRINSIC :: ISO_C_BINDING
!        TYPE(C_PTR) :: handle
!        INTEGER(C_INT), VALUE :: n, nd, n_n
!        REAL(C_DOUBLE) :: dvals(*)
!        INTEGER(C_INT) :: gnum(*), dinds(*), n_ranks(*)
!     END FUNCTION Permon_Initsolve
!
!     SUBROUTINE Permon_Solve( handle, x, b ) BIND(c,name='permon_solve')
!        USE, INTRINSIC :: ISO_C_BINDING
!        REAL(C_DOUBLE) :: x(*), b(*)
!        TYPE(C_PTR), VALUE :: handle
!     END SUBROUTINE Permon_solve
!  END INTERFACE

  IF ( PRESENT(Free_Fact) ) THEN
    IF ( Free_Fact ) THEN
      RETURN
    END IF
  END IF

  Factorize = ListGetLogical( Solver % Values, 'Linear System Refactorize', stat )
  IF ( .NOT. stat ) Factorize = .TRUE.

  IF ( Factorize .OR. .NOT.C_ASSOCIATED(A % PermonSolverInstance) ) THEN
    IF ( C_ASSOCIATED(A % PermonSolverInstance) ) THEN
       CALL Fatal( 'Permon', 're-entry not implemented' )
    END IF

    nd = COUNT(A % ConstrainedDOF)
    ALLOCATE(DirichletInds(nd), DirichletVals(nd))
    j = 0
    DO i=1,A % NumberOfRows
      IF(A % ConstrainedDOF(i)) THEN
        j = j + 1
        DirichletInds(j) = i; DirichletVals(j) = A % Dvalues(i)
      END IF
    END DO

    !TODO sequential case not working
    n = 0
    ALLOCATE(neighbours(Parenv % PEs))
    DO i=1,ParEnv % PEs
      IF( ParEnv % IsNeighbour(i) .AND. i-1/=ParEnv % myPE) THEN
        n = n + 1
        neighbours(n) = i-1
      END IF
    END DO

    !A % PermonSolverInstance = Permon_InitSolve( SIZE(A % ParallelInfo % GlobalDOFs), &
    !     A % ParallelInfo % GlobalDOFs, nd,  DirichletInds, DirichletVals, n, neighbours )

    IF( n_dof_partition /=  SIZE(A % ParallelInfo % GlobalDOFs) ) THEN
      CALL Fatal( 'Permon', &
        'inconsistency: A % NumberOfRows /=  SIZE(A % ParallelInfo % GlobalDOFs' )
    END IF

    ALLOCATE(IntOptions(10))
    ALLOCATE(RealOptions(1))

    CALL FETI4ISetDefaultIntegerOptions(IntOptions)
    CALL FETI4ISetDefaultRealOptions(RealOptions)

    CALL FETI4ICreateInstance(A % PermonSolverInstance, A % PermonMatrix, &
      A % NumberOfRows, b, A % ParallelInfo % GlobalDOFs, &
      n, neighbours, &
      nd, DirichletInds, DirichletVals, &
      IntOptions, RealOptions)
  END IF

  !CALL Permon_Solve( A % PermonSolverInstance, x, b )
  CALL FETI4ISolve(A % PermonSolverInstance, n_dof_partition, x)
#else
   CALL Fatal( 'Permon_SolveSystem', 'Permon Solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE Permon_SolveSystem
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Solves a linear system using Pardiso direct solver (from either MKL or
!> official Pardiso distribution. If possible, MKL-version is used).
!------------------------------------------------------------------------------
  SUBROUTINE Pardiso_SolveSystem( Solver,A,x,b,Free_fact )
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp), TARGET :: x(*), b(*)
    LOGICAL, OPTIONAL :: Free_fact

! MKL version of Pardiso (interface is different)
#if defined(HAVE_MKL)
    INTERFACE
      SUBROUTINE pardiso(pt, maxfct, mnum, mtype, phase, n, &
                           values, rows, cols, perm, nrhs, iparm, msglvl, b, x, ierror)
        USE Types
        IMPLICIT NONE
        REAL(KIND=dp) :: values(*), b(*), x(*)
        INTEGER(KIND=AddrInt) :: pt(*)
        INTEGER :: perm(*), nrhs, iparm(*), msglvl, ierror
        INTEGER :: maxfct, mnum, mtype, phase, n, rows(*), cols(*)
      END SUBROUTINE pardiso

      SUBROUTINE pardisoinit(pt, mtype, iparm)
        USE Types
        IMPLICIT NONE
        INTEGER(KIND=AddrInt) :: pt(*)
        INTEGER :: mtype
        INTEGER :: iparm(*)
      END SUBROUTINE pardisoinit
    END INTERFACE

    INTEGER maxfct, mnum, mtype, phase, n, nrhs, ierror, msglvl
    INTEGER, POINTER :: Iparm(:)
    INTEGER i, j, k, nz, idum(1), nzutd
    LOGICAL :: Found, matsym, matpd
    REAL*8  :: ddum(1)

    LOGICAL :: Factorize, FreeFactorize
    INTEGER :: tlen, allocstat
    CHARACTER(:), ALLOCATABLE :: mat_type

    REAL(KIND=dp), POINTER :: values(:)
    INTEGER, POINTER  :: rows(:), cols(:)

    ! Check if system needs to be refactorized
    Factorize = ListGetLogical( Solver % Values, &
                  'Linear System Refactorize', Found )
    IF ( .NOT. Found ) Factorize = .TRUE.

    ! Set matrix type for Pardiso
    mat_type = ListGetString(Solver % Values,'Linear System Matrix Type',Found)
    
    IF (Found) THEN
      SELECT CASE(mat_type)
      CASE('positive definite')
        mtype = 2
      CASE('symmetric indefinite')
        mtype = -2
      CASE('structurally symmetric')
        mtype = 1
      CASE('nonsymmetric', 'general')
        mtype = 11
      CASE DEFAULT
        mtype = 11
      END SELECT
    ELSE
      ! Check if matrix is symmetric or spd
      matsym = ListGetLogical(Solver % Values, &
                            'Linear System Symmetric', Found)

      matpd = ListGetLogical(Solver % Values, &
                    'Linear System Positive Definite', Found)

      IF (matsym) THEN
        IF (matpd) THEN
          ! Matrix is symmetric positive definite
           mtype = 2
         ELSE
          ! Matrix is structurally symmetric (can't handle indefinite systems!!!!)
          mtype = 1
        END IF
      ELSE
        ! Matrix is unsymmetric
        mtype = 11
      END IF
    END IF

    ! Free factorization if requested
    IF ( PRESENT(Free_Fact) ) THEN
      IF ( Free_Fact ) THEN
        IF(ASSOCIATED(A % PardisoId)) THEN
          phase     = -1           ! release internal memory
          n = A % NumberOfRows
          maxfct = 1
          mnum   = 1
          nrhs   = 1
          msglvl = 0
          CALL pardiso(A % PardisoId, maxfct, mnum, mtype, phase, n, &
              ddum, idum, idum, idum, nrhs, A % PardisoParam, msglvl, ddum, ddum, ierror)
          DEALLOCATE(A % PardisoId, A % PardisoParam)
          A % PardisoId => NULL()
          A % PardisoParam => NULL()
        END IF
        RETURN
      END IF
    END IF

    ! Get number of rows and number of nonzero elements
    n = A % Numberofrows
    nz = A % Rows(n+1)-1

    ! Copy upper triangular part of symmetric positive definite system
    IF ( ABS(mtype) == 2 ) THEN
      nzutd = 0
      DO i=1,n
        nzutd = nzutd + A % Rows(i+1)-A % Diag(i)
      END DO
      
      ALLOCATE( values(nzutd), cols(nzutd), rows(n+1), STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal('Pardiso_SolveSystem', &
           'Memory allocation for row and column indices failed')
      END IF

      ! Copy upper triangular part of A to Pardiso structure
      Rows(1)=1
      DO i=1,n
        ! Set up row pointers and copy values
        nzutd = A % Rows(i+1)-A % Diag(i)
        Rows(i+1) = Rows(i)+nzutd
        DO j=0,nzutd-1
          Cols(Rows(i)+j)=A % Cols(A%Diag(i)+j)
          Values(Rows(i)+j)=A % Values(A%Diag(i)+j)
        END DO
      END DO
      
    ELSE
      Cols => A % Cols
      Rows => A % Rows
      Values => A % Values
    END IF

    ! Set up parameters
    msglvl    = 0 ! Do not write out any info
    maxfct    = 1 ! Set up space for 1 matrix at most
    mnum      = 1 ! Matrix to use in the solution phase (1st and only one)
    nrhs      = 1 ! Use only one RHS
    iparm => A % PardisoParam

    ! Compute factorization if necessary
    IF ( Factorize .OR. .NOT.ASSOCIATED(A % PardisoID) ) THEN
      ! Free factorization
      IF (ASSOCIATED(A % PardisoId)) THEN
        phase = -1 ! release internal memory
        CALL pardiso (A % PardisoId, maxfct, mnum, mtype, phase, n, &
                    ddum, idum, idum, idum, nrhs, iparm, msglvl, ddum, ddum, ierror)
        DEALLOCATE(A % PardisoId, A % PardisoParam)
        A % PardisoId => NULL()
        A % PardisoParam => NULL()
      END IF

      ! Allocate pardiso internal structures
      ALLOCATE(A % PardisoId(64), A % PardisoParam(64), STAT=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal('Pardiso_SolveSystem', &
                     'Memory allocation for Pardiso failed')
      END IF
      iparm => A % PardisoParam
      iparm = 0
      A % PardisoId = 0
 
      ! Set up scaling values for solver based on matrix type
      CALL pardisoinit(A % PardisoId, mtype, iparm)

      ! Set up rest of parameters explicitly
      iparm(1)=1      ! Do not use solver default parameters
      iparm(2)=2      ! Minimize fill-in with nested dissection from Metis
      iparm(3)=0      ! Reserved
      iparm(4)=0      ! Always compute factorization
      iparm(5)=0      ! No user input permutation
      iparm(6)=0      ! Write solution vector to x
      iparm(8)=5      ! Number of iterative refinement steps
      iparm(18)=-1    ! Report nnz(L) and nnz(U)
      iparm(19)=0     ! Do not report Mflops
      iparm(27)=0     ! Do not check sparse matrix representation
      iparm(35)=0     ! Use Fortran style indexing
      iparm(60)=0     ! Use in-core version of Pardiso

      ! Perform analysis
      phase = 11            ! Analysis
      CALL pardiso(A % PardisoId, maxfct, mnum, mtype, phase, n, &
            values, rows, cols, idum, nrhs, iparm, msglvl, ddum, ddum, ierror)

      IF (ierror /= 0) THEN
        WRITE(*,'(A,I0)') 'MKL Pardiso: ERROR=', ierror
        CALL Fatal('Pardiso_SolveSystem','Error during analysis phase')
      END IF

      ! Perform factorization
      phase = 22            ! Factorization
      CALL pardiso (A % pardisoId, maxfct, mnum, mtype, phase, n, &
           values, rows, cols, idum, nrhs, iparm, msglvl, ddum, ddum, ierror)

      IF (ierror /= 0) THEN
        WRITE(*,'(A,I0)') 'MKL Pardiso: ERROR=', ierror
        CALL Fatal('Pardiso_SolveSystem','Error during factorization phase')
      END IF
    END IF ! Compute factorization

    ! Perform solve
    phase = 33            ! Solve, iterative refinement
    CALL pardiso(A % PardisoId, maxfct, mnum, mtype, phase, n, &
                values, rows, cols, idum, nrhs, iparm, msglvl, b, x, ierror)

    IF (ierror /= 0) THEN
      WRITE(*,'(A,I0)') 'MKL Pardiso: ERROR=', ierror
      CALL Fatal('Pardiso_SolveSystem','Error during solve phase')
    END IF

    ! Release memory if needed
    FreeFactorize = ListGetLogical( Solver % Values, &
                       'Linear System Free Factorization', Found )
    IF ( .NOT. Found ) FreeFactorize = .TRUE.

    IF ( Factorize .AND. FreeFactorize ) THEN
      phase     = -1           ! release internal memory
      CALL pardiso (A % PardisoId, maxfct, mnum, mtype, phase, n, &
            ddum, idum, idum, idum, nrhs, iparm, msglvl, ddum, ddum, ierror)

      DEALLOCATE(A % PardisoId, A % PardisoParam)
      A % PardisoId => NULL()
      A % PardisoParam => NULL()
    END IF
    IF (ABS(mtype) == 2) DEALLOCATE(Values, Rows, Cols)

! Distribution version of Pardiso
#elif defined(HAVE_PARDISO)
!..   All other variables

      INTERFACE
        SUBROUTINE pardiso(pt, maxfct, mnum, mtype, phase, n, &
          values, rows, cols, idum, nrhs, iparm, msglvl, b, x, ierror, dparm)
          USE Types
          REAL(KIND=dp) :: values(*), b(*), x(*), dparm(*)
          INTEGER(KIND=AddrInt) :: pt(*)
          INTEGER :: idum(*), nrhs, iparm(*), msglvl, ierror
          INTEGER :: maxfct, mnum, mtype, phase, n, rows(*), cols(*)
        END SUBROUTINE pardiso

        SUBROUTINE pardisoinit(pt,mtype,solver,iparm,dparm,ierror)
          USE Types
          INTEGER :: mtype, iparm(*),ierror,solver
          REAL(KIND=dp) :: dparm(*)
          INTEGER(KIND=AddrInt) :: pt(*)
        END SUBROUTINE pardisoinit
      END INTERFACE

      INTEGER maxfct, mnum, mtype, phase, n, nrhs, ierror, msglvl
      INTEGER, POINTER :: Iparm(:)
      INTEGER i, j, k, nz, idum(1)
      LOGICAL :: Found, Symm, Posdef
      REAL*8  waltime1, waltime2, ddum(1), dparm(64)

      LOGICAL :: Factorize, FreeFactorize
      INTEGER :: tlen
      CHARACTER(LEN=16) :: threads

      REAL(KIND=dp), POINTER :: values(:)
      INTEGER, POINTER  :: rows(:), cols(:)

      IF ( PRESENT(Free_Fact) ) THEN
        IF ( Free_Fact ) THEN
          IF(ASSOCIATED(A % PardisoId)) THEN
            phase     = -1           ! release internal memory
            n = A % NumberOfRows
            mtype  = 11
            maxfct = 1
            mnum   = 1
            nrhs   = 1
            msglvl = 0
            CALL pardiso(A % PardisoId, maxfct, mnum, mtype, phase, n, &
             ddum, idum, idum, idum, nrhs, A % PardisoParam, msglvl, ddum, ddum, ierror,dparm)
            DEALLOCATE(A % PardisoId, A % PardisoParam)
            A % PardisoId => NULL()
            A % PardisoParam => NULL()
          END IF
          RETURN
        END IF
      END IF

!  .. Setup Pardiso control parameters und initialize the solvers
!     internal address pointers. This is only necessary for the FIRST
!     call of the PARDISO solver.
!
      Factorize = ListGetLogical( Solver % Values, &
         'Linear System Refactorize', Found )
      IF ( .NOT. Found ) Factorize = .TRUE.

      symm = ListGetLogical( Solver % Values, &
         'Linear System Symmetric', Found )

      posdef = ListGetLogical( Solver % Values, &
         'Linear System Positive Definite', Found )

      mtype = 11

      cols => a % cols
      rows => a % rows
      values => a % values
      n = A % Numberofrows

      IF ( posdef ) THEN
        nz = A % rows(n+1)-1
        allocate( values(nz), cols(nz), rows(n+1) )
        k = 1
        do i=1,n
         rows(i)=k
         do j=a % rows(i), a % rows(i+1)-1
           if ( a % cols(j)>=i .and. a % values(j) /= 0._dp ) then
             cols(k) = a % cols(j)
             values(k) = a % values(j)
             k = k + 1
           end if
         end do
        end do
        rows(n+1)=k
        mtype     = 2
      ELSE IF ( symm ) THEN
        mtype = 1
      END IF

      msglvl    = 0       ! with statistical information
      maxfct    = 1
      mnum      = 1
      nrhs      = 1
      iparm => A % PardisoParam


      IF ( Factorize .OR. .NOT.ASSOCIATED(A % PardisoID) ) THEN
        IF(ASSOCIATED(A % PardisoId)) THEN
          phase     = -1           ! release internal memory
          CALL pardiso (A % PardisoId, maxfct, mnum, mtype, phase, n, &
           ddum, idum, idum, idum, nrhs, iparm, msglvl, ddum, ddum, ierror,dparm)
          DEALLOCATE(A % PardisoId, A % PardisoParam)
          A % PardisoId => NULL()
          A % PardisoParam => NULL()
        END IF

        ALLOCATE(A % PardisoId(64), A % PardisoParam(64))
        iparm => A % PardisoParam
        CALL pardisoinit(A % PardisoId, mtype, 0, iparm, dparm, ierror )

!..     Reordering and Symbolic Factorization, This step also allocates
!       all memory that is necessary for the factorization
!
        phase     = 11      ! only reordering and symbolic factorization
        msglvl    = 0       ! with statistical information
        maxfct    = 1
        mnum      = 1
        nrhs      = 1

!  ..   Numbers of Processors ( value of OMP_NUM_THREADS )
        CALL envir( 'OMP_NUM_THREADS', threads, tlen )

        iparm(3) = 1
        IF ( tlen>0 ) &
          READ(threads(1:tlen),*) iparm(3)
        IF (iparm(3)<=0) Iparm(3) = 1

        CALL pardiso(A % PardisoId, maxfct, mnum, mtype, phase, n, &
          values, rows, cols, idum, nrhs, iparm, msglvl, ddum, ddum, ierror, dparm)

        IF (ierror /= 0) THEN
          WRITE(*,*) 'The following ERROR was detected: ', ierror
          STOP EXIT_ERROR
        END IF

!..     Factorization.
        phase     = 22  ! only factorization
        CALL pardiso (A % pardisoId, maxfct, mnum, mtype, phase, n, &
         values, rows, cols, idum, nrhs, iparm, msglvl, ddum, ddum, ierror, dparm)

        IF (ierror /= 0) THEN
           WRITE(*,*) 'The following ERROR was detected: ', ierror
          STOP EXIT_ERROR
        ENDIF
      END IF

!..   Back substitution and iterative refinement
      phase     = 33  ! only factorization
      iparm(8)  = 0   ! max number of iterative refinement steps
                      ! (if perturbations occur in the factorization
                      ! phase, two refinement steps are taken)

      CALL pardiso(A % PardisoId, maxfct, mnum, mtype, phase, n, &
       values, rows, cols, idum, nrhs, iparm, msglvl, b, x, ierror, dparm)

!..   Termination and release of memory
      FreeFactorize = ListGetLogical( Solver % Values, &
          'Linear System Free Factorization', Found )
      IF ( .NOT. Found ) FreeFactorize = .TRUE.

      IF ( Factorize .AND. FreeFactorize ) THEN
        phase     = -1           ! release internal memory
        CALL pardiso (A % PardisoId, maxfct, mnum, mtype, phase, n, &
          ddum, idum, idum, idum, nrhs, iparm, msglvl, ddum, ddum, ierror, dparm)

        DEALLOCATE(A % PardisoId, A % PardisoParam)
        A % PardisoId => NULL()
        A % PardisoParam => NULL()
      END IF

      IF(posdef) DEALLOCATE( values, rows, cols )
#else
      CALL Fatal( 'Parsido_SolveSystem', 'Pardiso solver has not been installed.' )
#endif

!------------------------------------------------------------------------------
  END SUBROUTINE Pardiso_SolveSystem
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Solves a linear system using Cluster Pardiso direct solver from MKL
!------------------------------------------------------------------------------
  SUBROUTINE CPardiso_SolveSystem( Solver,A,x,b,Free_fact )
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp), TARGET :: x(*), b(*)
    LOGICAL, OPTIONAL :: Free_fact

! Cluster Pardiso
#if defined(HAVE_MKL) && defined(HAVE_CPARDISO)
    INTERFACE
        SUBROUTINE cluster_sparse_solver(pt, maxfct, mnum, mtype, phase, n, &
              values, rows, cols, perm, nrhs, iparm, msglvl, b, x, comm, ierror)
            USE Types
            REAL(KIND=dp) :: values(*), b(*), x(*)
            INTEGER(KIND=AddrInt) :: pt(*)
            INTEGER :: perm(*), nrhs, iparm(*), msglvl, ierror
            INTEGER :: maxfct, mnum, mtype, phase, n, rows(*), cols(*), comm
        END SUBROUTINE cluster_sparse_solver
    END INTERFACE

    INTEGER :: phase, n, ierror
    INTEGER, POINTER :: Iparm(:)
    INTEGER i, j, k, nz, nzutd, idum(1), nl, nt
    LOGICAL :: Found, matsym, matpd
    REAL(kind=dp) :: ddum(1)
    REAL(kind=dp), POINTER, DIMENSION(:) :: dbuf

    LOGICAL :: Factorize, FreeFactorize
    INTEGER :: tlen, allocstat

    REAL(KIND=dp), POINTER CONTIG :: values(:)
    INTEGER, POINTER CONTIG :: rows(:), cols(:)

    ! Free factorization if requested
    IF ( PRESENT(Free_Fact) ) THEN
        IF ( Free_Fact ) THEN
            CALL CPardiso_Free(A)
        END IF

        RETURN
    END IF

    ! Check if system needs to be refactorized
    Factorize = ListGetLogical( Solver % Values, &
                                'Linear System Refactorize', Found )
    IF ( .NOT. Found ) Factorize = .TRUE.

    ! Compute factorization if necessary
    IF ( Factorize .OR. .NOT.ASSOCIATED(A % CPardisoID) ) THEN
        CALL CPardiso_Factorize(Solver, A)
    END IF

    ! Get global start and end of domain
    nl = A % CPardisoId % iparm(41)
    nt = A % CPardisoId % iparm(42)

    ! Gather RHS
    A % CPardisoId % rhs = 0D0
    DO i=1,A % NumberOfRows
        A % CPardisoId % rhs(A % Gorder(i)-nl+1) = b(i)
    END DO

    ! Perform solve
    phase = 33      ! Solve, iterative refinement
    CALL cluster_sparse_solver(A % CPardisoId % ID, &
          A % CPardisoId % maxfct, A % CPardisoId % mnum, &
          A % CPardisoId % mtype,  phase, A % CPardisoId % n, &
          A % CPardisoID % aa, A % CPardisoID % ia, A % CPardisoID % ja, idum, &
          A % CPardisoId % nrhs, A % CPardisoID % iparm, &
          A % CPardisoId % msglvl, A % CPardisoId % rhs, &
          A % CPardisoId % x,  A % Comm, ierror)

    IF (ierror /= 0) THEN
        WRITE(*,'(A,I0)') 'MKL CPardiso: ERROR=', ierror
        CALL Fatal('CPardiso_SolveSystem','Error during solve phase')
    END IF

    ! Distribute solution
    DO i=1,A % NumberOfRows
        x(i)=A % CPardisoId % x(A % Gorder(i)-nl+1)
    END DO

    ! Release memory if needed
    FreeFactorize = ListGetLogical( Solver % Values, &
                                    'Linear System Free Factorization', Found )
    IF ( .NOT. Found ) FreeFactorize = .TRUE.

    IF ( Factorize .AND. FreeFactorize ) THEN
        CALL CPardiso_Free(A)
    END IF

#else
    CALL Fatal( 'CParsido_SolveSystem', 'Cluster Pardiso solver has not been installed.' )
#endif
!------------------------------------------------------------------------------
  END SUBROUTINE CPardiso_SolveSystem
!------------------------------------------------------------------------------

#if defined(HAVE_MKL) && defined(HAVE_CPARDISO)
  SUBROUTINE CPardiso_Factorize(Solver, A)
    IMPLICIT NONE
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A

    INTERFACE
        SUBROUTINE cluster_sparse_solver(pt, maxfct, mnum, mtype, phase, n, &
              values, rows, cols, perm, nrhs, iparm, msglvl, b, x, comm, ierror)
            USE Types
            REAL(KIND=dp) :: values(*), b(*), x(*)
            INTEGER(KIND=AddrInt) :: pt(*)
            INTEGER :: perm(*), nrhs, iparm(*), msglvl, ierror
            INTEGER :: maxfct, mnum, mtype, phase, n, rows(*), cols(*), comm
        END SUBROUTINE cluster_sparse_solver
    END INTERFACE

    LOGICAL :: matsym, matpd, Found
    INTEGER :: i, j, k, rind, lrow, rptr, rsize, lind, tind
    INTEGER :: allocstat, n, nz, nzutd, nl, nt, nOwned, nhalo, ierror
    INTEGER :: phase, idum(1)
    REAL(kind=dp) :: ddum(1)
    REAL(kind=dp), DIMENSION(:), POINTER :: aa
    INTEGER, DIMENSION(:), POINTER CONTIG :: iparm, ia, ja, owner, dsize, iperm, Order

    INTEGER :: fid
    CHARACTER(:), ALLOCATABLE :: mat_type

    ! Free old factorization if necessary
    IF (ASSOCIATED(A % CPardisoId)) THEN
        CALL CPardiso_Free(A)
    END IF

    ! Allocate Pardiso structure
    ALLOCATE(A % CPardisoId, STAT=allocstat)
    IF (allocstat /= 0) THEN
    CALL Fatal('CPardiso_Factorize', &
                   'Memory allocation for CPardiso failed')
    END IF
    ! Allocate control structures
    ALLOCATE(A % CPardisoId% ID(64), A % CPardisoId % IParm(64), STAT=allocstat)
    IF (allocstat /= 0) THEN
        CALL Fatal('CPardiso_Factorize', &
                   'Memory allocation for CPardiso failed')
    END IF

    iparm => A % CPardisoId % IParm
    ! Initialize control parameters and solver Id's
    DO i=1,64
        iparm(i)=0
    END DO
    DO i=1,64
        A % CPardisoId % ID(i)=0
    END DO

    ! Set matrix type for CPardiso
    mat_type = ListGetString(Solver % Values,'Linear System Matrix Type',Found)
    IF (Found) THEN
      SELECT CASE(mat_type)
      CASE('positive definite')
        A % CPardisoID % mtype = 2
      CASE('symmetric indefinite')
        A % CPardisoID % mtype = -2
      CASE('structurally symmetric')
        A % CPardisoID % mtype = 1
      CASE('nonsymmetric', 'general')
        A % CPardisoID % mtype = 11
      CASE DEFAULT
        A % CPardisoID % mtype = 11
      END SELECT
    ELSE
      ! Check if matrix is symmetric or spd
      matsym = ListGetLogical(Solver % Values, &
                            'Linear System Symmetric', Found)

      matpd = ListGetLogical(Solver % Values, &
                    'Linear System Positive Definite', Found)

      IF (matsym) THEN
        IF (matpd) THEN
          ! Matrix is symmetric positive definite
          A % CPardisoID % mtype = 2
        ELSE
          ! Matrix is structurally symmetric
          A % CPardisoID % mtype = 1
        END IF
      ELSE
        ! Matrix is nonsymmetric
        A % CPardisoID % mtype = 11
      END IF
    END IF

    ! Set up continuous numbering for the whole computation domain
    n = SIZE(A % ParallelInfo % GlobalDOFs)
    ALLOCATE(A % Gorder(n), Owner(n), STAT=allocstat)
    IF (allocstat /= 0) THEN
         CALL Fatal('CPardiso_Factorize', &
                    'Memory allocation for CPardiso global numbering failed')
    END IF
    CALL ContinuousNumbering(A % ParallelInfo, A % Perm, A % Gorder, Owner, nOwn=nOwned)

    ! Compute the number of global dofs
    CALL MPI_ALLREDUCE(nOwned, A % CPardisoId % n, &
                       1, MPI_INTEGER, MPI_SUM, A % Comm, ierror)
    DEALLOCATE(Owner)

    ! Find bounds of domain
    nl = A % Gorder(1)
    nt = A % Gorder(1)
    DO i=2,n
        ! NOTE: Matrix is structurally symmetric
        rind = A % Gorder(i)
        nl = MIN(rind, nl)
        nt = MAX(rind, nt)
    END DO

    ! Allocate temp storage for global numbering
    ALLOCATE(Order(n), iperm(n), STAT=allocstat)
    IF (allocstat /= 0) THEN
        CALL Fatal('CPardiso_Factorize', &
                    'Memory allocation for CPardiso global numbering failed')
    END IF

    ! Sort global numbering to build matrix
    Order(1:n) = A % Gorder(1:n)
    DO i=1,n
        iperm(i)=i
    END DO
    CALL SortI(n, Order, iperm)

    ! Allocate storage for CPardiso matrix
    nhalo = (nt-nl+1)-n
    nz = A % Rows(A % NumberOfRows+1)-1
    IF (ABS(A % CPardisoID % mtype) == 2) THEN
        nzutd = ((nz-n)/2)+1 + n
        ALLOCATE(A % CPardisoID % ia(nt-nl+2), &
                 A % CPardisoId % ja(nzutd+nhalo), &
                 A % CPardisoId % aa(nzutd+nhalo), &
                 STAT=allocstat)
    ELSE
        ALLOCATE(A % CPardisoID % ia(nt-nl+2), &
                 A % CPardisoId % ja(nz+nhalo), &
                 A % CPardisoId % aa(nz+nhalo), &
                STAT=allocstat)
    END IF
    IF (allocstat /= 0) THEN
        CALL Fatal('CPardiso_Factorize', &
                   'Memory allocation for CPardiso matrix failed')
    END IF
    ia => A % CPardisoId % ia
    ja => A % CPardisoId % ja
    aa => A % CPardisoId % aa

    ! Build distributed CRS matrix
    ia(1) = 1
    lrow = 1      ! Next row to add
    rptr = 1      ! Pointer to next row to add, equals ia(lrow)
    lind = Order(1)-1 ! Row pointer for the first round
    
    ! Add rows of matrix 
    DO i=1,n
      ! Add empty rows until the beginning of the row to add
      ! (first round adds nothing due to choice of lind)
      tind = Order(i)
      rsize = (tind-lind)-1
      
      ! Put zeroes to the diagonal
      DO j=1,rsize
        ia(lrow+j)=rptr+j
        ja(rptr+(j-1))=lind+j
        aa(rptr+(j-1))=0D0
      END DO
      ! Set up row pointers
      rptr = rptr + rsize
      lrow = lrow + rsize
      
      ! Add next row
      rind = iperm(i)
      lind = A % rows(rind)
      tind = A % rows(rind+1)
      IF (ABS(A % CPardisoId % mtype) == 2) THEN
        ! Add only upper triangular elements for symmetric matrices
        rsize = 0
        DO j=lind, tind-1
          IF (A % Gorder(A % Cols(j)) >= Order(i)) THEN
            ja(rptr+rsize)=A % Gorder(A % Cols(j))
            aa(rptr+rsize)=A % values(j)
            rsize = rsize + 1
          END IF
        END DO

      ELSE
        rsize = tind-lind
        DO j=lind, tind-1
          ja(rptr+(j-lind))=A % Gorder(A % Cols(j))
          aa(rptr+(j-lind))=A % values(j)
        END DO
      END IF
        
      ! Sort column indices
      CALL SortF(rsize, ja(rptr:rptr+rsize), aa(rptr:rptr+rsize))
        
      ! Set up row pointers
      rptr = rptr + rsize
      lrow = lrow + 1
      ia(lrow) = rptr

      lind = Order(i) ! Store row index for next round
    END DO

    ! Deallocate temp storage
    DEALLOCATE(Order, iperm)

    ! Set up parameters
    A % CPardisoId % msglvl    = 0 ! Do not write out = 0 / write out = 1 info
    A % CPardisoId % maxfct    = 1 ! Set up space for 1 matrix at most
    A % CPardisoId % mnum      = 1 ! Matrix to use in the solution phase (1st and only one)
    A % CPardisoId % nrhs      = 1 ! Use only one RHS
    ALLOCATE(A % CPardisoId % rhs(nt-nl+1), &
             A % CPardisoId % x(nt-nl+1), STAT=allocstat)
    IF (allocstat /= 0) THEN
        CALL Fatal('CPardiso_Factorize', &
                   'Memory allocation for CPardiso rhs and solution vector x failed')
    END IF

    ! Set up parameters explicitly
    iparm(1)=1          ! Do not use = 1 / use = 0 solver default parameters
    iparm(2)=2          ! Minimize fill-in with OpenMP nested dissection
    iparm(5)=0          ! No user input permutation
    iparm(6)=0          ! Write solution vector to x
    iparm(8)=0          ! Number of iterative refinement steps
    IF (A % CPardisoID % mtype ==11 .OR. &
        A % CPardisoID % mtype == 13) THEN
      iparm(10)=13      ! Perturbation value 10^-iparm(10) in case of small pivots
      iparm(11)=1       ! Use scalings from symmetric weighted matching
      iparm(13)=1       ! Use permutations from nonsymmetric weighted matching
    ELSE
      iparm(10)=8       ! Perturbation value 10^-iparm(10) in case of small pivots
      iparm(11)=0       ! Do not use scalings from symmetric weighted matching
      iparm(13)=0       ! Do not use permutations from symmetric weighted matching
    END IF
    
    iparm(21)=1         ! Do not use Bunch Kaufman pivoting
    iparm(27)=0         ! Do not check sparse matrix representation
    iparm(28)=0         ! Use double precision
    iparm(35)=0         ! Use Fortran indexing

    ! CPardiso matrix input format
    iparm(40) = 2       ! Distributed solution phase, distributed solution vector
    iparm(41) = nl      ! Beginning of solution domain
    iparm(42) = nt      ! End of solution domain

    ! Perform analysis
    phase = 11      ! Analysis
    CALL cluster_sparse_solver(A % CPardisoId % ID, &
          A % CPardisoId % maxfct, A % CPardisoId % mnum, &
          A % CPardisoId % mtype, phase, A % CPardisoId % n, &
          aa, ia, ja, idum, A % CPardisoId % nrhs, iparm, &
          A % CPardisoId % msglvl, &
          ddum, ddum, A % Comm, ierror)
    IF (ierror /= 0) THEN
        WRITE(*,'(A,I0)') 'MKL CPardiso: ERROR=', ierror
        CALL Fatal('CPardiso_SolveSystem','Error during analysis phase')
    END IF

    ! Perform factorization
    phase = 22      ! Factorization
    CALL cluster_sparse_solver(A % CPardisoId % ID, &
          A % CPardisoId % maxfct, A % CPardisoId % mnum, &
          A % CPardisoId % mtype, phase, A % CPardisoId % n, &
          aa, ia, ja, idum, A % CPardisoId % nrhs, iparm, &
          A % CPardisoId % msglvl, &
          ddum, ddum,  A % Comm, ierror)
    IF (ierror /= 0) THEN
        WRITE(*,'(A,I0)') 'MKL CPardiso: ERROR=', ierror
        CALL Fatal('CPardiso_SolveSystem','Error during factorization phase')
    END IF
  END SUBROUTINE CPardiso_Factorize


  SUBROUTINE CPardiso_Free(A)
    IMPLICIT NONE

    TYPE(Matrix_t) :: A
    INTERFACE
        SUBROUTINE cluster_sparse_solver(pt, maxfct, mnum, mtype, phase, n, &
                           values, rows, cols, perm, nrhs, iparm, &
                           msglvl, b, x, comm, ierror)
            USE Types
            REAL(KIND=dp) :: values(*), b(*), x(*)
            INTEGER(KIND=AddrInt) :: pt(*)
            INTEGER :: perm(*), nrhs, iparm(*), msglvl, ierror
            INTEGER :: maxfct, mnum, mtype, phase, n, rows(*), cols(*), comm
        END SUBROUTINE cluster_sparse_solver
    END INTERFACE

    INTEGER :: ierror, phase, idum(1)
    REAL(kind=dp) :: ddum(1)

    IF(ASSOCIATED(A % CPardisoId)) THEN
        phase     = -1           ! release internal memory
        CALL cluster_sparse_solver(A % CPardisoId % ID, &
              A % CPardisoID % maxfct, A % CPardisoID % mnum, &
              A % CPardisoID % mtype, phase, A % CPardisoID % n, &
              A % CPardisoId % aa, A % CPardisoId % ia, A % CPardisoId % ja, &
              idum, A % CPardisoID % nrhs, A % CPardisoID % IParm, &
              A % CPardisoID % msglvl, ddum, ddum, A % Comm, ierror)

        ! Deallocate global ordering
        DEALLOCATE(A % Gorder)

        ! Deallocate control structures
        DEALLOCATE(A % CPardisoId % ID)
        DEALLOCATE(A % CPardisoID % IParm)
        ! Deallocate matrix and permutation
        DEALLOCATE(A % CPardisoId % ia, A % CPardisoId % ja, &
                   A % CPardisoId % aa, A % CPardisoID % rhs, &
                   A % CPardisoID % x)
        ! Deallocate CPardiso structure
        DEALLOCATE(A % CPardisoID)
    END IF
  END SUBROUTINE CPardiso_Free
#endif


!------------------------------------------------------------------------------
  SUBROUTINE DirectSolver( A,x,b,Solver,Free_Fact )
!------------------------------------------------------------------------------

    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: x(*),b(*)
    TYPE(Matrix_t) :: A
    LOGICAL, OPTIONAL :: Free_Fact
!------------------------------------------------------------------------------

    LOGICAL :: GotIt
    CHARACTER(:), ALLOCATABLE :: Method
!------------------------------------------------------------------------------

    IF ( PRESENT(Free_Fact) ) THEN
      IF ( Free_Fact ) THEN
        CALL BandSolver( A, x, b, Free_Fact )
        CALL ComplexBandSolver( A, x, b, Free_Fact )
#ifdef HAVE_MUMPS
        CALL FreeMumpsFactorizations(A)
        CALL MumpsLocal_SolveSystem( Solver, A, x, b, Free_Fact )
#endif
#if defined(HAVE_MKL) || defined(HAVE_PARDISO)
        CALL Pardiso_SolveSystem( Solver, A, x, b, Free_Fact )
#endif
#if defined(HAVE_MKL) && defined(HAVE_CPARDISO)
        CALL CPardiso_SolveSystem( Solver, A, x, b, Free_Fact )
#endif
#ifdef HAVE_SUPERLU
        CALL SuperLU_SolveSystem( Solver, A, x, b, Free_Fact )
#endif
#ifdef HAVE_UMFPACK
        CALL Umfpack_SolveSystem( Solver, A, x, b, Free_Fact )
#endif
#ifdef HAVE_CHOLMOD
        CALL SPQR_SolveSystem( Solver, A, x, b, Free_Fact )
        CALL Cholmod_SolveSystem( Solver, A, x, b, Free_Fact )
#endif
#ifdef HAVE_FETI4I
        CALL Permon_SolveSystem( Solver, A, x, b, Free_Fact )
#endif
        RETURN
      END IF
    END IF

    Method=ListGetString(Solver % Values,'Linear System Direct Method',GotIt)
    IF ( .NOT. GotIt ) Method = 'banded'
    
    
    CALL Info('DirectSolver','Using direct method: '//Method,Level=9)

#if !defined (HAVE_UMFPACK) && defined (HAVE_MUMPS)
    IF ( Method == 'umfpack' .OR. Method == 'big umfpack' ) THEN
      CALL Warn( 'CheckLinearSolverOptions', 'UMFPACK solver not installed, using MUMPS instead!' )
      Method = 'mumps'
    END IF
#endif

    SELECT CASE(Method)
      CASE( 'banded', 'symmetric banded' )
        IF ( .NOT. A % Complex ) THEN
           CALL BandSolver( A, x, b )
        ELSE
           CALL ComplexBandSolver( A, x, b )
        END IF

      CASE( 'umfpack', 'big umfpack' )
        CALL Umfpack_SolveSystem( Solver, A, x, b )

      CASE( 'cholmod' )
        CALL Cholmod_SolveSystem( Solver, A, x, b )

      CASE( 'spqr' )
        CALL SPQR_SolveSystem( Solver, A, x, b )

      CASE( 'smumps', 'cmumps' )
        IF( A % Complex ) THEN
          CALL CMumps_SolveSystem( Solver, A, x, b )
        ELSE
          CALL SMumps_SolveSystem( Solver, A, x, b )
        END IF

      CASE( 'mumps', 'dmumps', 'zmumps' )
        IF( A % Complex ) THEN
          CALL ZMumps_SolveSystem( Solver, A, x, b )
        ELSE
          CALL Mumps_SolveSystem( Solver, A, x, b )
        END IF

      CASE( 'mumpslocal' )
        IF( A % Complex ) THEN
          CALL ZMumpsLocal_SolveSystem( Solver, A, x, b )
        ELSE
          CALL MumpsLocal_SolveSystem( Solver, A, x, b )
        END IF
          
      CASE( 'superlu' )
        CALL SuperLU_SolveSystem( Solver, A, x, b )

      CASE( 'permon' )
        CALL Permon_SolveSystem( Solver, A, x, b )

      CASE( 'pardiso' )
        CALL Pardiso_SolveSystem( Solver, A, x, b )

      CASE( 'cpardiso' )
        CALL CPardiso_SolveSystem( Solver, A, x, b )

      CASE DEFAULT
        CALL Fatal( 'DirectSolver', 'Unknown direct solver method.' )
    END SELECT

    ! We should be able to trust that a direct strategy will return a converged
    ! linear system.
    IF( ASSOCIATED( Solver % Variable ) ) THEN
      Solver % Variable % LinConverged = 1
    END IF
    
!------------------------------------------------------------------------------
  END SUBROUTINE DirectSolver
!------------------------------------------------------------------------------

END MODULE DirectSolve

!> \} ElmerLib
