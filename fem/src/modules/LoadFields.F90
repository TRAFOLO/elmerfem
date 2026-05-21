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
! *  Last mod Date: April 30, 2026 - added iGSE post-processing combining logic
! *                                  restored array vectorization for DG mapping
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Module for reading and processing losses computed by electromagnetic solvers.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE LoadFields_init( Model,Solver,dt,TransientSimulation )
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
  CHARACTER(LEN=MAX_NAME_LEN) :: EMSim, LossMethodStr
  TYPE(ValueList_t), POINTER :: SolverParams
  LOGICAL :: Found

!------------------------------------------------------------------------------

  SolverParams => GetSolverParams()

  EMSim = GetString(Solver % Values, 'EM Simulation Type', Found)
  
  IF (Found) THEN
    SELECT CASE (EMSim)
    CASE ('transient')
      LossMethodStr = GetString( SolverParams, 'Transient Core Loss Method', Found )
      IF (.NOT. Found) LossMethodStr = 'gse'

      CALL ListAddString( SolverParams, &
           NextFreeKeyword('Exported Variable ',SolverParams), '-dg LoadWkgAvg' )
      CALL ListAddString( SolverParams, &
           NextFreeKeyword('Exported Variable ',SolverParams), '-dg LoadWm3Avg' )
           
      ! Declare variables as -nooutput so they load into RAM from the restart file 
      ! but are NOT written to the final thermal VTU output.
      IF (LossMethodStr == 'igse') THEN
        CALL ListAddString( SolverParams, &
             NextFreeKeyword('Exported Variable ',SolverParams), '-nooutput -dg iGSE_Base1_WkgAvg' )
        CALL ListAddString( SolverParams, &
             NextFreeKeyword('Exported Variable ',SolverParams), '-nooutput -dg iGSE_Base2_WkgAvg' )
        CALL ListAddString( SolverParams, &
             NextFreeKeyword('Exported Variable ',SolverParams), '-nooutput -dg DeltaB' )
      END IF
      
    CASE ('harmonic')
      CALL ListAddString( SolverParams, &
             NextFreeKeyword('Exported Variable ',SolverParams), '-dg LoadWkg' )
      CALL ListAddString( SolverParams, &
             NextFreeKeyword('Exported Variable ',SolverParams), '-dg LoadWm3' )
      CALL ListAddString( SolverParams, &
             NextFreeKeyword('Exported Variable ',SolverParams), '-dg BT' )
    CASE ('none')
      CALL ListAddString( SolverParams, &
             NextFreeKeyword('Exported Variable ',SolverParams), '-dg LoadWkg' )
      CALL ListAddString( SolverParams, &
             NextFreeKeyword('Exported Variable ',SolverParams), '-dg LoadWm3' )
    CASE DEFAULT
        CALL Fatal ('LoadFields', 'Non existent electromagnetic simulation type!')
    END SELECT
  ELSE
    CALL Fatal ('LoadFields', 'EM Simulation Type not found! It should be harmonic or transient.')
  END IF

END SUBROUTINE LoadFields_init
!------------------------------------------------------------------------------
SUBROUTINE LoadFields( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
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
  TYPE(Mesh_t), POINTER :: Mesh
  CHARACTER(LEN=MAX_NAME_LEN) :: RestartFile, OutputName, EMSim, LossMethodStr
  INTEGER :: Visit = 0, nF, i, FieldSize, n, NSIZE
  INTEGER, POINTER :: Frequencies(:)
  LOGICAL :: Found
  
  TYPE(ValueList_t), POINTER :: SolverParams, Material
  TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp) :: FreqPower, FieldPower
  REAL(KIND=dp), POINTER :: VarFdP(:), VarFqP(:), VarDens(:)
  
  TYPE(Variable_t),POINTER :: TmpVar1, TmpVar2, TmpVar3, VarLoadWkg, VarLoadWm3, VarBT
  TYPE(Variable_t),POINTER :: VarBase1, VarBase2, VarDeltaB
  
  SAVE Visit, EMSim

!------------------------------------------------------------------------------
!    Check if version number output is requested
!------------------------------------------------------------------------------

  Mesh => Solver % Mesh
  SolverParams => GetSolverParams()
  
  IF ( Visit == 0 ) THEN

    EMSim = GetString(Solver % Values, 'EM Simulation Type', Found)
    
    IF (Found) THEN
      SELECT CASE (EMSim)
      CASE ('harmonic')
        WRITE( Message, *) 'Electromagnetic results will be loaded as Harmonic'
        CALL Info( 'LoadFields', Message, Level=4 )
        
        Frequencies =>  ListGetIntegerArray( Solver % Values, 'Frequencies', Found )
        IF(.NOT.ASSOCIATED(Frequencies)) THEN
          CALL Fatal ('LoadFields', 'Frequencies array not provided for harmonic mode')
        ELSE
          WRITE( Message, *) 'Frequencies: ', Frequencies
          CALL Info( 'LoadFields', Message, Level=4 )
        END IF
        
        nF = SIZE(Frequencies)
        
        WRITE( Message, *) 'Field size: ', FieldSize
        CALL Info( 'LoadFields', Message, Level=4 )
        
         ALLOCATE( TmpVar1 )
         ALLOCATE( TmpVar2 )
         ALLOCATE( TmpVar3 )
            
        DO i=1,nF
          RestartFile = 'h'//TRIM(I2S(Frequencies(i)))//'.result'
          OutputName = TRIM(Mesh % Name) // '/' // TRIM(RestartFile)
          IF(ParEnv % PEs > 1) OutputName = TRIM(OutputName) // '.' // TRIM(i2s(ParEnv % MyPe))
          
          WRITE( Message, *) 'Opening file: ', TRIM(OutputName)
          CALL LoadRestartFile( OutputName, i, Mesh)
          
          VarLoadWkg => VariableGet( Mesh % Variables, 'LoadWkg', ThisOnly = .TRUE. )
          IF(.NOT. ASSOCIATED( VarLoadWkg ) ) THEN
            CALL Fatal('LoadFields','Variable LoadWkg does not exist in Elmer mesh!')
          END IF
          
          VarLoadWm3 => VariableGet( Mesh % Variables, 'LoadWm3', ThisOnly = .TRUE. )
          IF(.NOT. ASSOCIATED( VarLoadWm3 ) ) THEN
            CALL Fatal('LoadFields','Variable LoadWm3 does not exist in Elmer mesh!')
          END IF
          
          VarBT => VariableGet( Mesh % Variables, 'BT', ThisOnly = .TRUE. )
          IF(.NOT. ASSOCIATED( VarBT ) ) THEN
            CALL Fatal('LoadFields','Variable BT does not exist in Elmer mesh!')
          END IF
          
          IF(SIZE(VarLoadWkg % Values) .NE. SIZE(TmpVar1 % Values)) THEN
            ALLOCATE( TmpVar1 % Values(SIZE(VarLoadWkg % Values)) )
            ALLOCATE( TmpVar2 % Values(SIZE(VarLoadWkg % Values)) )
            ALLOCATE( TmpVar3 % Values(SIZE(VarLoadWkg % Values)) )
            TmpVar1 % Values = 0.0
            TmpVar2 % Values = 0.0
            TmpVar3 % Values = 0.0
          END IF
          
          TmpVar1 % Values = TmpVar1 % Values + VarLoadWkg % Values
          TmpVar2 % Values = TmpVar2 % Values + VarLoadWm3 % Values
          TmpVar3 % Values = TmpVar3 % Values + VarBT % Values
        END DO
        
        VarLoadWkg % Values = TmpVar1 % Values
        VarLoadWm3 % Values = TmpVar2 % Values
        VarBT % Values = TmpVar3 % Values
        
      CASE ('transient')
        LossMethodStr = GetString( SolverParams, 'Transient Core Loss Method', Found )
        IF (.NOT. Found) LossMethodStr = 'gse'

        WRITE( Message, *) 'Electromagnetic results will be loaded as Transient'
        CALL Info( 'LoadFields', Message, Level=4 )
        
          RestartFile = 'h1.result'
          OutputName = TRIM(Mesh % Name) // '/' // TRIM(RestartFile)
          IF(ParEnv % PEs > 1) OutputName = TRIM(OutputName) // '.' // TRIM(i2s(ParEnv % MyPe))
          
          WRITE( Message, *) 'Opening file: ', TRIM(OutputName)
          CALL LoadRestartFile( OutputName, i, Mesh)
          
          VarLoadWkg => VariableGet( Mesh % Variables, 'LoadWkgAvg', ThisOnly = .TRUE. )
          IF(.NOT. ASSOCIATED( VarLoadWkg ) ) THEN
            CALL Fatal('LoadFields','Variable LoadWkgAvg does not exist in Elmer mesh!')
          END IF
          
          VarLoadWm3 => VariableGet( Mesh % Variables, 'LoadWm3Avg', ThisOnly = .TRUE. )
          IF(.NOT. ASSOCIATED( VarLoadWm3 ) ) THEN
            CALL Fatal('LoadFields','Variable LoadWm3Avg does not exist in Elmer mesh!')
          END IF
          
          ! ----------------------------------------------------
          ! Apply iGSE Core Loss Math (if enabled)
          ! ----------------------------------------------------
          IF (LossMethodStr == 'igse') THEN
            VarBase1 => VariableGet( Mesh % Variables, 'iGSE_Base1_WkgAvg', ThisOnly = .TRUE. )
            IF(.NOT. ASSOCIATED( VarBase1 ) ) CALL Fatal('LoadFields','iGSE_Base1_WkgAvg missing!')

            VarBase2 => VariableGet( Mesh % Variables, 'iGSE_Base2_WkgAvg', ThisOnly = .TRUE. )
            IF(.NOT. ASSOCIATED( VarBase2 ) ) CALL Fatal('LoadFields','iGSE_Base2_WkgAvg missing!')

            VarDeltaB => VariableGet( Mesh % Variables, 'DeltaB', ThisOnly = .TRUE. )
            IF(.NOT. ASSOCIATED( VarDeltaB ) ) CALL Fatal('LoadFields','DeltaB missing!')

            NSIZE = SIZE( VarLoadWkg % Values )
            ALLOCATE( VarFdP(NSIZE), VarFqP(NSIZE), VarDens(NSIZE) )

            ! Iterate to map the material exponents & density per active element to DG fields
            DO i = 1, GetNOFActive()
              Element => Mesh % Elements(i)
              Material => GetMaterial( Element )

              ! Strict Density Check - assign directly to DG array to avoid rank mismatch
              VarDens(VarLoadWkg % Perm(Element % DGIndexes ) ) = GetReal( Material, 'Density', Found, UElement = Element )
              IF (.NOT. Found) THEN
                CALL Fatal('LoadFields', 'Density not found for Body: '//TRIM(i2s(Element % BodyId)))
              END IF

              ! Get Exponents - fallback to 0.0 for air/windings
              FreqPower = ListGetCReal( Material,'Harmonic Loss Linear Frequency Exponent',Found )
              IF (.NOT. Found) FreqPower = 0.0d0
              
              FieldPower = ListGetCReal( Material,'Harmonic Loss Linear Exponent',Found )
              IF (.NOT. Found) FieldPower = 0.0d0

              ! Map scalars dynamically onto the element's Discontinuous Galerkin array indices
              VarFdP(VarLoadWkg % Perm(Element % DGIndexes ) ) = FieldPower
              VarFqP(VarLoadWkg % Perm(Element % DGIndexes ) ) = FreqPower
            END DO

            ! Vectorized math: Add the iGSE core loss to the total loss array (which currently holds just Joule)
            VarLoadWkg % Values = VarLoadWkg % Values + &
                                  VarBase1 % Values * (VarDeltaB % Values**(VarFdP - VarFqP)) + &
                                  VarBase2 % Values

            ! Explicitly calculate the W/m3 equivalent
            VarLoadWm3 % Values = VarDens * VarLoadWkg % Values

            DEALLOCATE( VarFdP, VarFqP, VarDens )
            
            WRITE( Message, *) 'iGSE core losses successfully combined with Joule losses.'
            CALL Info( 'LoadFields', Message, Level=4 )
          END IF

      CASE ('none')
        WRITE( Message, *) 'Electromagnetic results not used.'
        CALL Info( 'LoadFields', Message, Level=4 )
      CASE DEFAULT
        CALL Fatal ('LoadFields', 'Non existent electromagnetic simulation type!')
      END SELECT
    ELSE
        CALL Fatal ('LoadFields', 'EM Simulation Type not found! It should be harmonic or transient.')
    END IF
  ELSE
    WRITE( Message, *) 'Transient mode'
    CALL Info( 'LoadFields', Message, Level=4 )
  END IF

!------------------------------------------------------------------------------

  Visit = Visit + 1
!------------------------------------------------------------------------------
END SUBROUTINE LoadFields
!------------------------------------------------------------------------------