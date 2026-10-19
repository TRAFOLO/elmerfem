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
! *  Authors: TRAFOLO Sp. z o.o.
! *  Email:   info@trafolo.eu
! *  Web:     https://trafolo.eu
! *  Address: Pl. Solny 15,
! *           50-062 Wrocław, Poland
! *
! *  Original Date: July 22, 2022
! *  Last mod Date: July  9, 2024  - added proxmity loss
! *                 April 30, 2026 - added HAm magnitude for harmonic cases
! *                                  added iGSE assuming only major loop
! *                                  replaced BTmax with native spatial DeltaB tracking
! *                 September 30, 2026 - Steinmetz thermal factor per node (DEV-1537)
! *                 October 18, 2026 - no core loss from the switch-on step; start time without core loss (DEV-1607)
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Module for reading and processing losses computed by electromagnetic solvers.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE ProcessFields_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------

  TYPE(ValueList_t), POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: LossMethodStr
  LOGICAL :: Found

!------------------------------------------------------------------------------

  SolverParams => GetSolverParams()

  LossMethodStr = GetString( SolverParams, 'Transient Core Loss Method', Found )
  IF (.NOT. Found) LossMethodStr = 'gse'

  CALL ListAddString( SolverParams, &
         NextFreeKeyword('Exported Variable ',SolverParams), '-dg LoadWkg' )
  CALL ListAddString( SolverParams, &
         NextFreeKeyword('Exported Variable ',SolverParams), '-dg LoadWm3' )
  CALL ListAddString( SolverParams, &
         NextFreeKeyword('Exported Variable ',SolverParams), '-dg BT' )

  IF (TransientSimulation) THEN
    ! Export Native Vector Components for Bounding Box tracking (kept in RAM)
    CALL ListAddString( SolverParams, &
         NextFreeKeyword('Exported Variable ',SolverParams), '-nooutput -dg Bmax[Bmax:3]' )
    CALL ListAddString( SolverParams, &
         NextFreeKeyword('Exported Variable ',SolverParams), '-nooutput -dg Bmin[Bmin:3]' )

    ! Export calculated spatial DeltaB
    CALL ListAddString( SolverParams, &
         NextFreeKeyword('Exported Variable ',SolverParams), '-dg DeltaB' )
         
    ! Universal averaging variables (always needed for copper/windings)
    CALL ListAddString( SolverParams, &
         NextFreeKeyword('Exported Variable ',SolverParams), '-dg LoadWkgAvg' )
    CALL ListAddString( SolverParams, &
         NextFreeKeyword('Exported Variable ',SolverParams), '-dg LoadWm3Avg' )

    IF (LossMethodStr == 'igse') THEN
      CALL ListAddString( SolverParams, &
           NextFreeKeyword('Exported Variable ',SolverParams), '-nooutput -dg iGSE_Base1_Wkg' )
      CALL ListAddString( SolverParams, &
           NextFreeKeyword('Exported Variable ',SolverParams), '-nooutput -dg iGSE_Base2_Wkg' )
      CALL ListAddString( SolverParams, &
           NextFreeKeyword('Exported Variable ',SolverParams), '-dg iGSE_Base1_WkgAvg' )
      CALL ListAddString( SolverParams, &
           NextFreeKeyword('Exported Variable ',SolverParams), '-dg iGSE_Base2_WkgAvg' )
    END IF
  END IF

END SUBROUTINE ProcessFields_init
!------------------------------------------------------------------------------
SUBROUTINE ProcessFields( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER :: Visit = 0, nF, i, n, dim, NSIZE, jPL
  LOGICAL :: Found, HAmFound = .FALSE.
  REAL(KIND=dp) :: FreqPower, FieldPower, LossCoeff, LossCoeff2, k1, k2
  REAL(KIND=dp) :: Told = 0.0d0, Tstart = 0.0d0
  REAL(KIND=dp), POINTER :: VarLCk1(:), VarLC2k2(:), VarFdP(:), VarFqP(:),  &
                            VarDens(:), dBdt(:), VarBoldX(:), VarBoldY(:),  &
                            VarBoldZ(:), VarTF(:)

  TYPE(Variable_t),POINTER :: VarJH, VarBT, TmpVar, VarLoadWkg, VarLoadWm3, &
                              VarDeltaB, VarLoadWkgAvg, VarLoadWm3Avg, VarLossLin, &
                              VarLossQuad, VarHAm, VarPL, VarPLelem, &
                              Var_iGSE_Base1, Var_iGSE_Base2, &
                              Var_iGSE_Base1Avg, Var_iGSE_Base2Avg, &
                              VarBmaxComp, VarBminComp
  TYPE(ValueList_t), POINTER :: Material, SolverParams
  TYPE(Element_t), POINTER :: Element
  CHARACTER(LEN=15) :: coeffStr
  CHARACTER(LEN=MAX_NAME_LEN) :: LossMethodStr
  ! Transient core loss of a 'layered' laminate stack: its loss law takes the
  ! steel's flux density, B_s = B_par/F + B_perp (B/F without a normal), on the
  ! metal fraction F of the volume, as the harmonic loss of CalcFields does.
  ! Per DG node: the stacking factor (1 outside layered stacks), the normal,
  ! and the flux density of the loss law; the exported B stays the stack's.
  REAL(KIND=dp), POINTER :: LamFill(:), LamNrm(:,:), LossBX(:), LossBY(:), LossBZ(:), &
                            LossBAbs(:)
  LOGICAL, POINTER :: LamHasN(:)
  LOGICAL :: LamAny = .FALSE.

  SAVE Visit, Mesh, VarLCk1, VarLC2k2, VarFdP, VarFqP, VarDens, VarTF, LossCoeff, &
       LossCoeff2, Told, Tstart, VarBoldX, VarBoldY, VarBoldZ, dBdt, LossMethodStr, &
       LamFill, LamNrm, LamHasN, LamAny, LossBX, LossBY, LossBZ, LossBAbs

!------------------------------------------------------------------------------
!    Check if version number output is requested
!------------------------------------------------------------------------------

  Mesh => Solver % Mesh
  dim = CoordinateSystemDimension()
  SolverParams => GetSolverParams()

  VarBT => VariableGet( Mesh % Variables, 'BT', ThisOnly = .TRUE. )
  IF(.NOT. ASSOCIATED( VarBT ) ) THEN
    CALL Fatal('ProcessFields','Variable BT does not exist in Elmer mesh!')
  END IF

  VarHAm => VariableGet( Mesh % Variables, 'HAm', ThisOnly = .TRUE. )
  IF(ASSOCIATED( VarHAm ) ) HAmFound = .TRUE.

  VarLoadWkg => VariableGet( Mesh % Variables, 'LoadWkg', ThisOnly = .TRUE. )
  IF(.NOT. ASSOCIATED( VarLoadWkg ) ) THEN
    CALL Fatal('ProcessFields','Variable LoadWkg does not exist in Elmer mesh!')
  END IF

  VarLoadWm3 => VariableGet( Mesh % Variables, 'LoadWm3', ThisOnly = .TRUE. )
  IF(.NOT. ASSOCIATED( VarLoadWm3 ) ) THEN
    CALL Fatal('ProcessFields','Variable LoadWm3 does not exist in Elmer mesh!')
  END IF

  VarJH => VariableGet( Mesh % Variables, 'joule heating e', ThisOnly = .TRUE. )
  IF(.NOT. ASSOCIATED( VarJH ) ) THEN
    CALL Fatal('ProcessFields','Variable joule heating e does not exist in Elmer mesh!')
  END IF

  NSIZE = SIZE( VarLoadWkg % Values )

  IF (VISIT .EQ. 0) THEN
    ALLOCATE( VarDens(NSIZE) )

    DO i = 1, GetNOFActive()
      Element => Mesh % Elements(i)
      Material => GetMaterial( Element )

      VarDens(VarLoadWkg % Perm(Element % DGIndexes ) ) = GetReal( Material, 'Density', Found, UElement = Element )
      IF (.NOT. Found) THEN
        CALL Warn('ProcessFields', 'Density not found for Body: '//TRIM(i2s(Element % BodyId)))
        CALL Warn('ProcessFields', 'Setting default value 1.2 kg/m3')
        VarDens(VarLoadWkg % Perm(Element % DGIndexes ) ) = 1.2d0
      END IF
    END DO

    ! Steinmetz thermal correction poly(T), per node at the current Temperature field (the initial
    ! condition, or the previous thermal solve restarted for temperature iterations > 1). It can't
    ! live in the Harmonic Loss coefficients: CalcFields evaluates those with ListGetFun(..., Freq),
    ! which puts the frequency where the temperature should be (DEV-1537).
    ALLOCATE( VarTF(NSIZE) )
    VarTF = 1.0d0
    DO i = 1, GetNOFActive()
      Element => Mesh % Elements(i)
      Material => GetMaterial( Element )
      IF ( ListCheckPresent( Material, 'Harmonic Loss Temperature Factor' ) ) THEN
        VarTF(VarLoadWkg % Perm(Element % DGIndexes ) ) = &
            GetReal( Material, 'Harmonic Loss Temperature Factor', Found, UElement = Element )
      END IF
    END DO
    WRITE(Message, '(A,ES12.5,A,ES12.5)') 'Core loss temperature factor: min ', MINVAL(VarTF), &
        ', max ', MAXVAL(VarTF)
    CALL Info('ProcessFields', Message, Level=3 )
  END IF

  ! Reset LoadWkg for the current timestep and unconditionally add Joule & Proximity Loss if present
  VarLoadWkg % Values = VarJH % Values / VarDens

  VarPL => VariableGet( Mesh % Variables, 'proximity loss e', ThisOnly = .TRUE. )
  IF( ASSOCIATED( VarPL ) ) THEN
    VarLoadWkg % Values = VarLoadWkg % Values + VarPL % Values / VarDens
  END IF

  ! Transient winding homogenization writes an elemental (-elem) field named
  ! 'Proximity Loss': one W/m3 value per element. Broadcast it onto the
  ! element's DG nodes so the visualized, integrated and averaged loads (and
  ! through LoadWkgAvg the thermal heat source) include the proximity part.
  ! The TYPE guard is essential: in harmonic runs the DG export of CalcFields
  ! also registers a NODAL-averaged variable with this exact name, whose Perm
  ! must not be indexed by element - the harmonic case is already covered by
  ! the 'proximity loss e' fold above.
  VarPLelem => VariableGet( Mesh % Variables, 'proximity loss', ThisOnly = .TRUE. )
  IF( ASSOCIATED( VarPLelem ) ) THEN
   IF( VarPLelem % TYPE == Variable_on_elements ) THEN
    DO i = 1, GetNOFActive()
      Element => Mesh % Elements(i)
      jPL = VarPLelem % Perm(Element % ElementIndex)
      IF (jPL <= 0) CYCLE
      VarLoadWkg % Values(VarLoadWkg % Perm(Element % DGIndexes)) = &
          VarLoadWkg % Values(VarLoadWkg % Perm(Element % DGIndexes)) + &
          VarPLelem % Values(jPL) / VarDens(VarLoadWkg % Perm(Element % DGIndexes))
    END DO
   END IF
  END IF

  IF (Transient) THEN
    IF (VISIT .EQ. 0)  THEN
      ALLOCATE( VarBoldX(NSIZE), VarBoldY(NSIZE), VarBoldZ(NSIZE), dbdt(NSIZE) )
      ALLOCATE( LamFill(NSIZE), LamNrm(3,NSIZE), LamHasN(NSIZE), LossBX(NSIZE), LossBY(NSIZE), &
          LossBZ(NSIZE), LossBAbs(NSIZE) )
      LamFill = 1.0d0
      LamNrm = 0.0d0
      LamHasN = .FALSE.
      DO i = 1, GetNOFActive()
        Element => Mesh % Elements(i)
        Material => GetMaterial( Element )
        IF ( GetString( Material, 'Laminate Stack Model', Found ) /= 'layered' ) CYCLE
        CALL LayeredStackOfElement()
      END DO
      LamAny = ANY( LamFill < 1.0d0 )
      IF ( LamAny ) CALL Info( 'ProcessFields', 'Core loss of the layered laminate stacks '// &
          'from the steel flux density on their metal fraction', Level=3 )

      LossMethodStr = GetString( SolverParams, 'Transient Core Loss Method', Found )
      IF (.NOT. Found) LossMethodStr = 'gse'
      CALL Info ('ProcessFields', 'Transient Core Loss Method: '//TRIM(LossMethodStr), Level=3 )

      Tstart = GetCReal( SolverParams,'Averaging Start Time',Found )
      WRITE(coeffStr , '(ES12.5)') Tstart
      CALL Info ('ProcessFields', 'Averaging Start Time: '//TRIM(coeffStr), Level=3 )
    END IF

    VarBT % Values = 0.0d0

    DO i = 1,dim
      TmpVar => VariableGet( Mesh % Variables, 'magnetic flux density e '//TRIM(i2s(i)), ThisOnly = .TRUE. )
      VarBT % Values = VarBT % Values + TmpVar % Values * TmpVar % Values
    END DO
    VarBT % Values = SQRT(VarBT % Values)

    ! The flux density of the loss law, and at step 1 its history: no dB/dt for
    ! step 1, against the zero initial field it is the switch-on, not the waveform
    CALL LossFluxDensity()
    IF (VISIT .EQ. 0) THEN
      VarBoldX = LossBX
      VarBoldY = LossBY
      VarBoldZ = LossBZ
    END IF

    IF ( HAmFound ) THEN
      VarHAm % Values = 0.0d0
      DO i = 1,dim
        TmpVar => VariableGet( Mesh % Variables, 'magnetic field strength e '//TRIM(i2s(i)), ThisOnly = .TRUE. )
        IF( .NOT. ASSOCIATED( TmpVar ) ) CALL Fatal('ProcessFields','Enable Magnetic Field Strength!')
        VarHAm % Values = VarHAm % Values + TmpVar % Values * TmpVar % Values
      END DO
      VarHAm % Values = SQRT(VarHAm % Values)
    END IF


    IF ( GetLogical( SolverParams, 'Calculate Harmonic Loss', Found ) ) THEN
      IF (VISIT .EQ. 0) THEN
        ALLOCATE( VarLCk1(NSIZE), VarLC2k2(NSIZE), VarFdP(NSIZE), VarFqP(NSIZE) )

        DO i = 1, GetNOFActive()
          Element => Mesh % Elements(i)
          Material => GetMaterial( Element )

          FreqPower = ListGetCReal( Material,'Harmonic Loss Linear Frequency Exponent',Found )
          IF (Found) THEN
            FieldPower = ListGetCReal( Material,'Harmonic Loss Linear Exponent',Found )
            IF (.NOT. Found) CALL Fatal ('ProcessFields', 'Harmonic Loss Linear Exponent not found' )
          ELSE
            FreqPower = 1.0d0
            FieldPower = 1.0d0
          END IF

          IF ( LossMethodStr == 'igse' ) THEN
             ! iGSE implementation for k1 with integer fast-paths
             IF (abs(FreqPower-1.0d0) < 1.0d-6) THEN
                k1 = 1.0d0 / (2.0d0**(FieldPower + 1.0d0))
             ELSE IF (abs(FreqPower-2.0d0) < 1.0d-6) THEN
                k1 = 2.0d0 / ( (PI**2) * (2.0d0**FieldPower) )
             ELSE
                k1 = PI**(0.5d0 - FreqPower) * gamma(1.0d0 + FreqPower/2.0d0) / &
                     (2.0d0**FieldPower * gamma(0.5d0 + 0.5d0*FreqPower))
             END IF
          ELSE
             ! Original GSE implementation
             IF (abs(FreqPower-1.0d0) < 1.0d-6) THEN
               k1 = FieldPower / 4.0d0
             ELSE
               k1 = (2.0d0 * PI)**(1.0d0 - FreqPower) * 2.0d0 * gamma(0.5d0 - 0.5d0 * FreqPower) * &
                  gamma(1.0d0 + 0.5d0 * FieldPower) / gamma(0.5d0 * (1.0d0 - FreqPower + FieldPower)) / &
                  (gamma(0.5d0 - 0.5d0 * FreqPower) * gamma(0.5d0 * (1.0d0 + FreqPower)) + 3.0d0 * PI / cos(0.5d0 * FreqPower * PI))
             END IF
          END IF

          k2 = 0.5d0 / PI**2

          LossCoeff = ListGetCReal( Material,'Harmonic Loss Linear Coefficient',Found )
          IF (Found) THEN
            LossCoeff2 = ListGetCReal( Material,'Harmonic Loss Quadratic Coefficient',Found )
            IF (.NOT. Found) CALL Fatal ('ProcessFields', 'Harmonic Loss Quadratic Coefficient not found' )
          ELSE
            LossCoeff = 0.0d0
            LossCoeff2 = 0.0d0
          END IF

          VarLCk1(VarLoadWkg % Perm(Element % DGIndexes ) ) = LossCoeff * k1 * &
              VarTF(VarLoadWkg % Perm(Element % DGIndexes ) )
          VarLC2k2(VarLoadWkg % Perm(Element % DGIndexes ) ) = LossCoeff2 * k2 * &
              VarTF(VarLoadWkg % Perm(Element % DGIndexes ) )
          VarFdP(VarLoadWkg % Perm(Element % DGIndexes ) ) = FieldPower
          VarFqP(VarLoadWkg % Perm(Element % DGIndexes ) ) = FreqPower
        END DO

        CALL Info ('ProcessFields', 'Core loss coeffs: '//TRIM(coeffStr), Level=3 )
      END IF


      dBdt = SQRT( (LossBX - VarBoldX)**2 + (LossBY - VarBoldY)**2 + (LossBZ - VarBoldZ)**2 )/dt

      ! Apply Losses properly separated by method, on the metal fraction of a layered stack
      IF ( LossMethodStr == 'igse' ) THEN
        ! Base 1: iGSE Hysteresis Temporal Base
        Var_iGSE_Base1 => VariableGet( Mesh % Variables, 'iGSE_Base1_Wkg', ThisOnly = .TRUE. )
        IF (ASSOCIATED(Var_iGSE_Base1)) THEN
          Var_iGSE_Base1 % Values = LamFill * VarLCk1 * dBdt**VarFqP
        END IF

        ! Base 2: Eddy Current Base
        Var_iGSE_Base2 => VariableGet( Mesh % Variables, 'iGSE_Base2_Wkg', ThisOnly = .TRUE. )
        IF (ASSOCIATED(Var_iGSE_Base2)) THEN
          Var_iGSE_Base2 % Values = LamFill * VarLC2k2 * dBdt**2
        END IF
      ELSE
        ! Original GSE accumulates directly into LoadWkg (alongside Joule losses)
        IF ( LamAny ) THEN
          VarLoadWkg % Values = VarLoadWkg % Values + LamFill * ( &
                                VarLCk1 * dBdt**VarFqP * &
                                LossBAbs**(VarFdP-VarFqP) + &
                                VarLC2k2 * dBdt**2 )
        ELSE
          ! without layered stacks the terms as before, to the last bit
          VarLoadWkg % Values = VarLoadWkg % Values + &
                                VarLCk1 * dBdt**VarFqP * &
                                VarBT % Values**(VarFdP-VarFqP) + &
                                VarLC2k2 * dBdt**2
        END IF
      END IF

    ELSE
      CALL Info('ProcessFields','Core loss will not be computed!', Level=3 )
    END IF

    ! Compute LoadWm3 explicitly for the main variable
    VarLoadWm3 % Values = VarDens * VarLoadWkg % Values

    ! ==========================================
    ! UNIVERSAL AVERAGING BLOCK (Runs for all)
    ! ==========================================
    
    ! 1. Always average standard LoadWkg and LoadWm3
    VarLoadWkgAvg => VariableGet( Mesh % Variables, 'LoadWkgAvg', ThisOnly = .TRUE. )
    IF(.NOT. ASSOCIATED( VarLoadWkgAvg ) ) THEN
      CALL Fatal('ProcessFields','Variable LoadWkgAvg does not exist in Elmer mesh!')
    ELSE
      IF (Told + dt/2 > Tstart) THEN
        VarLoadWkgAvg % Values = ((Told-Tstart) * VarLoadWkgAvg % Values + dt * VarLoadWkg % Values)/(Told-Tstart+dt)
      ELSE
        VarLoadWkgAvg % Values = 0.0d0
      END IF
    END IF
    
    VarLoadWm3Avg => VariableGet( Mesh % Variables, 'LoadWm3Avg', ThisOnly = .TRUE. )
    IF(.NOT. ASSOCIATED( VarLoadWm3Avg ) ) THEN
      CALL Fatal('ProcessFields','Variable LoadWm3Avg does not exist in Elmer mesh!')
    ELSE
      VarLoadWm3Avg % Values = VarDens * VarLoadWkgAvg % Values
    END IF

    ! 2. If iGSE is active, ALSO average the isolated base variables
    IF ( LossMethodStr == 'igse' ) THEN
      ! Base 1 Averaging
      Var_iGSE_Base1Avg => VariableGet( Mesh % Variables, 'iGSE_Base1_WkgAvg', ThisOnly = .TRUE. )
      IF(.NOT. ASSOCIATED( Var_iGSE_Base1Avg ) ) THEN
        CALL Fatal('ProcessFields','Variable iGSE_Base1_WkgAvg does not exist in Elmer mesh!')
      ELSE
        IF (Told + dt/2 > Tstart) THEN
          Var_iGSE_Base1Avg % Values = ((Told-Tstart) * Var_iGSE_Base1Avg % Values + dt * Var_iGSE_Base1 % Values)/(Told-Tstart+dt)
        ELSE
          Var_iGSE_Base1Avg % Values = 0.0d0
        END IF
      END IF

      ! Base 2 Averaging
      Var_iGSE_Base2Avg => VariableGet( Mesh % Variables, 'iGSE_Base2_WkgAvg', ThisOnly = .TRUE. )
      IF(.NOT. ASSOCIATED( Var_iGSE_Base2Avg ) ) THEN
        CALL Fatal('ProcessFields','Variable iGSE_Base2_WkgAvg does not exist in Elmer mesh!')
      ELSE
        IF (Told + dt/2 > Tstart) THEN
          Var_iGSE_Base2Avg % Values = ((Told-Tstart) * Var_iGSE_Base2Avg % Values + dt * Var_iGSE_Base2 % Values)/(Told-Tstart+dt)
        ELSE
          Var_iGSE_Base2Avg % Values = 0.0d0
        END IF
      END IF
    END IF

    ! Save current B vectors (of the loss law) for next dB/dt derivative calculation
    VarBoldX = LossBX
    VarBoldY = LossBY
    VarBoldZ = LossBZ

    ! ==========================================
    ! UNIVERSAL BOUNDING BOX & DELTA B TRACKING
    ! ==========================================
    VarDeltaB => VariableGet( Mesh % Variables, 'DeltaB', ThisOnly = .TRUE. )
    
    IF (ASSOCIATED(VarDeltaB)) THEN
      VarDeltaB % Values = 0.0d0
      
      DO i = 1, dim
        ! Get current vector component, of the loss law's flux density
        TmpVar => VariableGet( Mesh % Variables, 'magnetic flux density e '//TRIM(i2s(i)), ThisOnly = .TRUE. )

        ! Get bounding box components
        VarBmaxComp => VariableGet( Mesh % Variables, 'Bmax '//TRIM(i2s(i)), ThisOnly = .TRUE. )
        VarBminComp => VariableGet( Mesh % Variables, 'Bmin '//TRIM(i2s(i)), ThisOnly = .TRUE. )

        IF (ASSOCIATED(VarBmaxComp) .AND. ASSOCIATED(VarBminComp) .AND. ASSOCIATED(TmpVar)) THEN
          IF (Visit == 0 .OR. Told + dt/2 <= Tstart) THEN
            ! Pin min/max to the current field until averaging starts, and at step 1 (not to B = 0)
            VarBmaxComp % Values = LossComponent(i)
            VarBminComp % Values = LossComponent(i)
          ELSE
            ! Track spatial maximum and minimum bounds
            VarBmaxComp % Values = MAX(VarBmaxComp % Values, LossComponent(i))
            VarBminComp % Values = MIN(VarBminComp % Values, LossComponent(i))
          END IF

          ! Calculate the squared spatial difference for this axis and add to running DeltaB sum
          VarDeltaB % Values = VarDeltaB % Values + (VarBmaxComp % Values - VarBminComp % Values)**2
        END IF
      END DO
      
      ! Finalize true spatial DeltaB magnitude
      VarDeltaB % Values = SQRT(VarDeltaB % Values)
    END IF

  ELSE ! HARMONIC
    IF (Visit>0) THEN
      CALL Fatal ('ProcessFields', 'For Harmonic B-field solver should be called once per simulation')
    END IF

    ! ----------------------------------------------------
    ! Harmonic Case: Calculate BT (Magnetic Flux Density)
    ! ----------------------------------------------------
    VarBT % Values = 0.0d0
    DO i = 1,dim
      TmpVar => VariableGet( Mesh % Variables, 'magnetic flux density re e '//TRIM(i2s(i)), ThisOnly = .TRUE. )
      VarBT % Values = VarBT % Values + TmpVar % Values * TmpVar % Values
    END DO

    DO i = 1,dim
      TmpVar => VariableGet( Mesh % Variables, 'magnetic flux density im e '//TRIM(i2s(i)), ThisOnly = .TRUE. )
      VarBT % Values = VarBT % Values + TmpVar % Values * TmpVar % Values
    END DO
    VarBT % Values = SQRT(VarBT % Values)

    ! ----------------------------------------------------
    ! Harmonic Case: Calculate HAm (Magnetic Field Strength)
    ! ----------------------------------------------------
    IF ( HAmFound ) THEN
      VarHAm % Values = 0.0d0
      DO i = 1,dim
        TmpVar => VariableGet( Mesh % Variables, 'magnetic field strength re e '//TRIM(i2s(i)), ThisOnly = .TRUE. )
        IF( .NOT. ASSOCIATED( TmpVar ) ) CALL Fatal('ProcessFields','Missing Magnetic Field Strength re!')
        VarHAm % Values = VarHAm % Values + TmpVar % Values * TmpVar % Values
      END DO

      DO i = 1,dim
        TmpVar => VariableGet( Mesh % Variables, 'magnetic field strength im e '//TRIM(i2s(i)), ThisOnly = .TRUE. )
        IF( .NOT. ASSOCIATED( TmpVar ) ) CALL Fatal('ProcessFields','Missing Magnetic Field Strength im!')
        VarHAm % Values = VarHAm % Values + TmpVar % Values * TmpVar % Values
      END DO
      VarHAm % Values = SQRT(VarHAm % Values)
    END IF
  
    ! ----------------------------------------------------
    ! Add Harmonic Losses
    ! ----------------------------------------------------
    VarLossLin => VariableGet( Mesh % Variables, 'harmonic loss linear e', ThisOnly = .TRUE. )
    VarLossQuad => VariableGet( Mesh % Variables, 'harmonic loss quadratic e', ThisOnly = .TRUE. )

    IF (ASSOCIATED( VarLossLin ) .AND. ASSOCIATED( VarLossQuad ) ) THEN
      VarLoadWkg % Values = VarLoadWkg % Values + VarTF * (VarLossLin % Values + VarLossQuad % Values)
    ELSE
      CALL Info('ProcessFields','Core loss will not be computed!', Level=3 )
    END IF

    VarLoadWm3 % Values = VarDens * VarLoadWkg % Values
  END IF

!------------------------------------------------------------------------------

  Visit = Visit + 1
  Told = Told + dt

CONTAINS

  ! Stacking factor and normal of the layered stack of Element at its DG nodes
  SUBROUTINE LayeredStackOfElement()
    REAL(KIND=dp) :: Fill, Nrm(3)
    REAL(KIND=dp), POINTER :: NormalArray(:,:)
    LOGICAL :: HasNormal, GotFill
    INTEGER :: k, j
    Fill = GetCReal( Material, 'Laminate Stacking Factor', GotFill )
    IF ( .NOT. GotFill ) Fill = 1.0d0
    IF ( Fill <= 0.0d0 .OR. Fill > 1.0d0 ) CALL Fatal( 'ProcessFields', &
        '> Laminate Stacking Factor < must be in (0,1]!' )
    Nrm = 0.0d0
    NormalArray => ListGetConstRealArray( Material, 'Laminate Normal', HasNormal )
    IF ( HasNormal ) THEN
      IF ( SIZE(NormalArray,1) < 3 ) CALL Fatal( 'ProcessFields', '> Laminate Normal < needs three components!' )
      Nrm = NormalArray(1:3,1)
      IF ( SQRT(SUM(Nrm**2)) < EPSILON(1.0d0) ) CALL Fatal( 'ProcessFields', '> Laminate Normal < is a zero vector!' )
      Nrm = Nrm / SQRT(SUM(Nrm**2))
    END IF
    DO k = 1, SIZE( Element % DGIndexes )
      j = VarLoadWkg % Perm( Element % DGIndexes(k) )
      IF ( j <= 0 ) CYCLE
      LamFill(j) = Fill
      LamNrm(:,j) = Nrm
      LamHasN(j) = HasNormal
    END DO
  END SUBROUTINE LayeredStackOfElement

  ! The flux density of the loss law at the DG nodes: B, in a layered stack the
  ! steel's B_par/F + B_perp; its magnitude as BT outside layered stacks
  SUBROUTINE LossFluxDensity()
    TYPE(Variable_t), POINTER :: BVar
    REAL(KIND=dp) :: bn
    INTEGER :: j
    BVar => VariableGet( Mesh % Variables, 'magnetic flux density e 1', ThisOnly = .TRUE. )
    LossBX = BVar % Values
    BVar => VariableGet( Mesh % Variables, 'magnetic flux density e 2', ThisOnly = .TRUE. )
    LossBY = BVar % Values
    BVar => VariableGet( Mesh % Variables, 'magnetic flux density e 3', ThisOnly = .TRUE. )
    LossBZ = BVar % Values
    LossBAbs = VarBT % Values
    IF ( .NOT. LamAny ) RETURN
    DO j = 1, NSIZE
      IF ( LamFill(j) >= 1.0d0 ) CYCLE
      IF ( LamHasN(j) ) THEN
        bn = LossBX(j) * LamNrm(1,j) + LossBY(j) * LamNrm(2,j) + LossBZ(j) * LamNrm(3,j)
        LossBX(j) = (LossBX(j) - bn * LamNrm(1,j)) / LamFill(j) + bn * LamNrm(1,j)
        LossBY(j) = (LossBY(j) - bn * LamNrm(2,j)) / LamFill(j) + bn * LamNrm(2,j)
        LossBZ(j) = (LossBZ(j) - bn * LamNrm(3,j)) / LamFill(j) + bn * LamNrm(3,j)
      ELSE
        LossBX(j) = LossBX(j) / LamFill(j)
        LossBY(j) = LossBY(j) / LamFill(j)
        LossBZ(j) = LossBZ(j) / LamFill(j)
      END IF
      LossBAbs(j) = SQRT( LossBX(j)**2 + LossBY(j)**2 + LossBZ(j)**2 )
    END DO
  END SUBROUTINE LossFluxDensity

  ! Component i of the flux density of the loss law
  FUNCTION LossComponent( i ) RESULT( Bi )
    INTEGER :: i
    REAL(KIND=dp) :: Bi(NSIZE)
    SELECT CASE( i )
    CASE( 1 )
      Bi = LossBX
    CASE( 2 )
      Bi = LossBY
    CASE DEFAULT
      Bi = LossBZ
    END SELECT
  END FUNCTION LossComponent
!------------------------------------------------------------------------------
END SUBROUTINE ProcessFields
!------------------------------------------------------------------------------