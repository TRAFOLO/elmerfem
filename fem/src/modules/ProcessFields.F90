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
                            VarBoldZ(:)

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

  SAVE Visit, Mesh, VarLCk1, VarLC2k2, VarFdP, VarFqP, VarDens, LossCoeff, &
       LossCoeff2, Told, Tstart, VarBoldX, VarBoldY, VarBoldZ, dBdt, LossMethodStr

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
      VarBoldX = 0.0d0
      VarBoldY = 0.0d0
      VarBoldZ = 0.0d0

      LossMethodStr = GetString( SolverParams, 'Transient Core Loss Method', Found )
      IF (.NOT. Found) LossMethodStr = 'gse'
      CALL Info ('ProcessFields', 'Transient Core Loss Method: '//TRIM(LossMethodStr), Level=3 )
    END IF

    VarBT % Values = 0.0d0

    DO i = 1,dim
      TmpVar => VariableGet( Mesh % Variables, 'magnetic flux density e '//TRIM(i2s(i)), ThisOnly = .TRUE. )
      VarBT % Values = VarBT % Values + TmpVar % Values * TmpVar % Values
    END DO
    VarBT % Values = SQRT(VarBT % Values)

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

          VarLCk1(VarLoadWkg % Perm(Element % DGIndexes ) ) = LossCoeff * k1
          VarLC2k2(VarLoadWkg % Perm(Element % DGIndexes ) ) = LossCoeff2 * k2
          VarFdP(VarLoadWkg % Perm(Element % DGIndexes ) ) = FieldPower
          VarFqP(VarLoadWkg % Perm(Element % DGIndexes ) ) = FreqPower
        END DO

        CALL Info ('ProcessFields', 'Core loss coeffs: '//TRIM(coeffStr), Level=3 )

        Tstart = GetCReal( SolverParams,'Averaging Start Time',Found )
        IF( .NOT. Found ) THEN
          CALL Info ('ProcessFields', 'Averaging Start Time: using default value 0.0s', Level=3 )
          Tstart = 0.0d0
        ELSE
          WRITE(coeffStr , '(F9.5)') Tstart
          CALL Info ('ProcessFields', 'Averaging Start Time: '//TRIM(coeffStr), Level=3 )
        END IF
      END IF


      TmpVar => VariableGet( Mesh % Variables, 'magnetic flux density e 1', ThisOnly = .TRUE. )
      dBdt = (TmpVar % Values - VarBoldX)**2

      TmpVar => VariableGet( Mesh % Variables, 'magnetic flux density e 2', ThisOnly = .TRUE. )
      dBdt = dBdt + (TmpVar % Values - VarBoldY)**2

      TmpVar => VariableGet( Mesh % Variables, 'magnetic flux density e 3', ThisOnly = .TRUE. )
      dBdt = dBdt + (TmpVar % Values - VarBoldZ)**2

      dBdt = SQRT(dBdt)/dt

      ! Apply Losses properly separated by method
      IF ( LossMethodStr == 'igse' ) THEN
        ! Base 1: iGSE Hysteresis Temporal Base
        Var_iGSE_Base1 => VariableGet( Mesh % Variables, 'iGSE_Base1_Wkg', ThisOnly = .TRUE. )
        IF (ASSOCIATED(Var_iGSE_Base1)) THEN
          Var_iGSE_Base1 % Values = VarLCk1 * dBdt**VarFqP
        END IF
        
        ! Base 2: Eddy Current Base
        Var_iGSE_Base2 => VariableGet( Mesh % Variables, 'iGSE_Base2_Wkg', ThisOnly = .TRUE. )
        IF (ASSOCIATED(Var_iGSE_Base2)) THEN
          Var_iGSE_Base2 % Values = VarLC2k2 * dBdt**2
        END IF
      ELSE
        ! Original GSE accumulates directly into LoadWkg (alongside Joule losses)
        VarLoadWkg % Values = VarLoadWkg % Values + &
                              VarLCk1 * dBdt**VarFqP * &
                              VarBT % Values**(VarFdP-VarFqP) + &
                              VarLC2k2 * dBdt**2
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

    ! Save current B vectors for next dB/dt derivative calculation
    TmpVar => VariableGet( Mesh % Variables, 'magnetic flux density e 1', ThisOnly = .TRUE. )
    VarBoldX = TmpVar % Values
    TmpVar => VariableGet( Mesh % Variables, 'magnetic flux density e 2', ThisOnly = .TRUE. )
    VarBoldY = TmpVar % Values
    TmpVar => VariableGet( Mesh % Variables, 'magnetic flux density e 3', ThisOnly = .TRUE. )
    VarBoldZ = TmpVar % Values

    ! ==========================================
    ! UNIVERSAL BOUNDING BOX & DELTA B TRACKING
    ! ==========================================
    VarDeltaB => VariableGet( Mesh % Variables, 'DeltaB', ThisOnly = .TRUE. )
    
    IF (ASSOCIATED(VarDeltaB)) THEN
      VarDeltaB % Values = 0.0d0
      
      DO i = 1, dim
        ! Get current vector component
        TmpVar => VariableGet( Mesh % Variables, 'magnetic flux density e '//TRIM(i2s(i)), ThisOnly = .TRUE. )
        
        ! Get bounding box components
        VarBmaxComp => VariableGet( Mesh % Variables, 'Bmax '//TRIM(i2s(i)), ThisOnly = .TRUE. )
        VarBminComp => VariableGet( Mesh % Variables, 'Bmin '//TRIM(i2s(i)), ThisOnly = .TRUE. )

        IF (ASSOCIATED(VarBmaxComp) .AND. ASSOCIATED(VarBminComp) .AND. ASSOCIATED(TmpVar)) THEN
          IF (Told + dt/2 <= Tstart) THEN
            ! Pin min/max to current field vector before averaging kicks in
            VarBmaxComp % Values = TmpVar % Values
            VarBminComp % Values = TmpVar % Values
          ELSE
            ! Track spatial maximum and minimum bounds
            VarBmaxComp % Values = MAX(VarBmaxComp % Values, TmpVar % Values)
            VarBminComp % Values = MIN(VarBminComp % Values, TmpVar % Values)
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
      VarLoadWkg % Values = VarLoadWkg % Values + VarLossLin % Values + VarLossQuad % Values
    ELSE
      CALL Info('ProcessFields','Core loss will not be computed!', Level=3 )
    END IF

    VarLoadWm3 % Values = VarDens * VarLoadWkg % Values
  END IF

!------------------------------------------------------------------------------

  Visit = Visit + 1
  Told = Told + dt
!------------------------------------------------------------------------------
END SUBROUTINE ProcessFields
!------------------------------------------------------------------------------