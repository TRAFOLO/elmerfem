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
! *  Authors:   Eelis Takala(Trafotek Oy) and Juha Ruokolainen(CSC)
! *  Emails:    eelis.takala@trafotek.fi and Juha.Ruokolainen@csc.fi
! *  Web:       http://www.trafotek.fi and http://www.csc.fi/elmer
! *  Addresses: Trafotek Oy
! *             Kaarinantie 700
! *             Turku
! *
! *             and
! *
! *             CSC - IT Center for Science Ltd.
! *             Keilaranta 14
! *             02101 Espoo, Finland 
! *
! *  Original Date: October 2015
! *
! *****************************************************************************/
 
MODULE CircuitUtils

    USE DefUtils
    IMPLICIT NONE

    ! Largest number of strand pieces one element may be cut into. An element
    ! straddles as many strands as the layout is fine; FoilSheetPieces is fatal
    ! above this instead of running off the buffer.
    INTEGER, PARAMETER :: MaxFoilSheetPieces = 1024

    !> Element quadrature of the foil sheet strand integrals, shared by the
    !> harmonic and the transient kernel and by the strand checks.
    TYPE FoilSheetQuad_t
      LOGICAL :: Exact = .TRUE.
      INTEGER :: nItem = 0
      INTEGER :: pCell(MaxFoilSheetPieces) = 0, pSeg(MaxFoilSheetPieces) = 0
      REAL(KIND=dp) :: pVol(MaxFoilSheetPieces) = 0._dp
      REAL(KIND=dp) :: pBary(4,MaxFoilSheetPieces) = 0._dp
      TYPE(GaussIntegrationPoints_t) :: IP
    END TYPE FoilSheetQuad_t

CONTAINS

!------------------------------------------------------------------------------
  FUNCTION GetCircuitModelDepth() RESULT (Depth)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    TYPE(Valuelist_t), POINTER :: simulation
    REAL(KIND=dp) :: depth
    LOGICAL :: Found, CSymmetry, Parallel
    INTEGER :: NoSlices
    
    CSymmetry = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric )

    simulation => GetSimulation()
    depth = GetConstReal(simulation, 'Circuit Model Depth', Found)

    IF( Found ) THEN
      NoSlices = ListGetInteger(simulation,'Number of Slices',Found)
      IF(NoSlices > 1) THEN
        IF( CurrentModel % Solver % Parallel ) depth = depth / NoSlices
      END IF
    ELSE
      depth = 1._dp
      IF (CSymmetry) depth = 2._dp * pi
    END IF
        
!------------------------------------------------------------------------------
  END FUNCTION GetCircuitModelDepth
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION GetComponentVoltageFactor(CompInd) RESULT (VoltageFactor)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: CompInd
    REAL(KIND=dp) :: VoltageFactor
    TYPE(Valuelist_t), POINTER :: CompParams
    LOGICAL :: Found
     
    CompParams => CurrentModel % Components(CompInd) % Values
    IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('GetComponentVoltageFactor',&
                                                        'Component parameters not found')
    VoltageFactor = GetConstReal(CompParams, 'Circuit Equation Voltage Factor', Found)
    IF (.NOT. Found) VoltageFactor = 1._dp
!------------------------------------------------------------------------------
  END FUNCTION GetComponentVoltageFactor
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION GetComponentParams(Element) RESULT (ComponentParams)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: i
    TYPE(Element_t), POINTER :: Element
    TYPE(Valuelist_t), POINTER :: ComponentParams, EntityParams
    LOGICAL :: Found
    
    EntityParams => GetBC(Element)
    IF (.NOT. ASSOCIATED(EntityParams)) THEN

      EntityParams => GetBodyParams( Element )
      IF (.NOT. ASSOCIATED(EntityParams)) CALL Fatal ('GetCompParams', 'Body Parameters not found')

    END IF
   
    i = GetInteger(EntityParams, 'Component', Found)
    
    IF (.NOT. Found) THEN
      ComponentParams => Null()
    ELSE
      ComponentParams => CurrentModel % Components(i) % Values
    END IF
   
!------------------------------------------------------------------------------
  END FUNCTION GetComponentParams
!------------------------------------------------------------------------------


! Get the current associated to the component from the solution of the
! constraint problem having lagrange multiplier vector. 
!------------------------------------------------------------------------------
  FUNCTION GetComponentCurrent(CompId,Found) RESULT ( Curr ) 
    INTEGER :: CompId
    LOGICAL :: Found
    COMPLEX(KIND=dp) :: Curr
    
    INTEGER :: i,j
    TYPE(CircuitVariable_t), POINTER :: iVar
    TYPE(Variable_t), POINTER :: LagrangeVar
    REAL(KIND=dp) :: CurrIm, CurrRe
    TYPE(Circuit_t), POINTER :: Circuit
    CHARACTER(LEN=MAX_NAME_LEN) :: str 
       
    Found = .FALSE.
    Curr = 0.0_dp

    IF(CurrentModel % n_Circuits == 0) RETURN
    
    CurrRe = 0.0_dp
    CurrIm = 0.0_dp
    
    DO i = 1, CurrentModel % n_Circuits     
      Circuit => CurrentModel % Circuits(i)
      
      str = LagrangeMultiplierName( CurrentModel % ASolver ) 
      LagrangeVar => VariableGet( CurrentModel % Mesh % Variables, str, ThisOnly = .TRUE.)
      IF(.NOT. ASSOCIATED(LagrangeVar) ) RETURN           
      
      DO j = 1, SIZE(Circuit % Components)
        ivar => Circuit % Components(j) % ivar
        IF(.NOT. ASSOCIATED(ivar)) CYCLE            
        IF(.NOT. iVar % isIvar ) CYCLE
        IF(iVar % BodyId /= CompId ) CYCLE
        IF(iVar % ValueId > 0 ) THEN
          Found = .TRUE.
          CurrRe = LagrangeVar % Values(iVar % ValueId)
        END IF
        IF(iVar % ImValueId > 0 ) THEN
          CurrIm = LagrangeVar % Values(iVar % ImValueId)
        END IF        
        IF(Found) EXIT
      END DO
      IF(Found) EXIT
    END DO

    IF(.NOT. Found) THEN
      CALL Fatal('GetComponentCurrent','Got circuits but no current for component: '//I2S(CompId))
    END IF
      
    !PRINT *,'Curr:',CompId,CurrRe,CurrIm
    Curr = CMPLX(CurrRe,CurrIm,KIND=dp)

  END FUNCTION GetComponentCurrent


  ! Get the current associated to the component from the solution of the
! constraint problem having lagrange multiplier vector. 
!------------------------------------------------------------------------------
  FUNCTION GetComponentArea(CompId,Found) RESULT ( Area ) 
    INTEGER :: CompId
    LOGICAL :: Found
    REAL(KIND=dp) :: Area
    
    INTEGER :: i,j
    TYPE(CircuitVariable_t), POINTER :: iVar
    TYPE(Variable_t), POINTER :: LagrangeVar
    REAL(KIND=dp) :: CurrIm, CurrRe
    TYPE(Circuit_t), POINTER :: Circuit
    CHARACTER(LEN=MAX_NAME_LEN) :: str 
       
    Found = .FALSE.
    Area = 0.0_dp

    IF(CurrentModel % n_Circuits == 0) RETURN
    
    DO i = 1, CurrentModel % n_Circuits     
      Circuit => CurrentModel % Circuits(i)      
      DO j = 1, SIZE(Circuit % Components)
        ivar => Circuit % Components(j) % ivar
        IF(iVar % BodyId /= CompId ) CYCLE

        Area = Circuit % Components(j) % ElArea
        Found = .TRUE.
        EXIT
      END DO
      IF(Found) EXIT
    END DO

    IF(.NOT. Found) THEN
      CALL Fatal('GetComponentArea','Got circuits but no area for component: '//I2S(CompId))
    END IF
      
  END FUNCTION GetComponentArea

  
!------------------------------------------------------------------------------
  FUNCTION GetComponentId(Element) RESULT (ComponentId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: ComponentId
    TYPE(Element_t), POINTER :: Element
    TYPE(Valuelist_t), POINTER :: BodyParams
    LOGICAL :: Found
    
    BodyParams => GetBodyParams( Element )
    IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('GetCompParams', 'Body Parameters not found')
   
    ComponentId = GetInteger(BodyParams, 'Component', Found)
    
!------------------------------------------------------------------------------
  END FUNCTION GetComponentId
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE GetWPotential(Wbase)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: Wbase(:)

    CALL GetLocalSolution(Wbase,'W Potential')
    IF(.NOT. ANY(Wbase/=0._dp)) CALL GetLocalSolution(Wbase,'W')
!------------------------------------------------------------------------------
  END SUBROUTINE GetWPotential
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE GetWPotentialVar(pVar)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Variable_t), POINTER :: pVar

    pVar => VariableGet( CurrentModel % Mesh % Variables,'W Potential')
    IF(.NOT. ASSOCIATED(pVar) ) THEN
      pVar => VariableGet( CurrentModel % Mesh % Variables,'W')
    END IF
    IF(ASSOCIATED(pVar)) THEN
      CALL Info('GetWPotentialVar','Using gradient of field to define direction: '&
          //TRIM(pVar % Name),Level=7)
    ELSE
      CALL Warn('GetWPotentialVar','Could not obtain variable for potential "W"')
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GetWPotentialVar
!------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
! Nodal direction potential of a coil component, normalized so that grad of it
! has unit circulation along every current path of the coil.
! Open coil: the electrode potential W of WPotentialSolver. Closed coil (a full
! ring, no electrodes, no single valued W): the CoilSolver cut potential, which
! is CoilPot or CoilPotB depending on PotSelect exactly as in the CoilSolver's
! own LocalFluxMatrix, divided by the multiplier ScalePotential applied to it.
!------------------------------------------------------------------------------
  SUBROUTINE GetCoilWBase(Element, nn, CompParams, Wbase, Wpot)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    INTEGER :: nn
    TYPE(ValueList_t), POINTER :: CompParams
    REAL(KIND=dp) :: Wbase(:)
    TYPE(Variable_t), POINTER, OPTIONAL :: Wpot

    TYPE(Variable_t), POINTER, SAVE :: PotVarA => NULL(), PotVarB => NULL(), &
        PotSelect => NULL()
    LOGICAL, SAVE :: CoilVarsVisited = .FALSE.
    TYPE(Variable_t), POINTER :: PotVar
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(KIND=dp) :: Mult, Circ
    INTEGER :: i, j
    LOGICAL :: Found, UseB
    CHARACTER(*), PARAMETER :: Caller = 'GetCoilWBase'
!------------------------------------------------------------------------------

    IF( .NOT. ListGetLogical( CompParams,'Coil Closed', Found ) ) THEN
      IF( PRESENT( Wpot ) ) THEN
        IF( ASSOCIATED( Wpot ) ) THEN
          CALL GetLocalSolution( Wbase, UElement=Element, UVariable=Wpot, Found=Found )
          RETURN
        END IF
      END IF
      CALL GetWPotential( Wbase )
      RETURN
    END IF

    IF( .NOT. CoilVarsVisited ) THEN
      CoilVarsVisited = .TRUE.
      PotVarA => VariableGet( CurrentModel % Mesh % Variables,'CoilPot' )
      PotVarB => VariableGet( CurrentModel % Mesh % Variables,'CoilPotB' )
      PotSelect => VariableGet( CurrentModel % Mesh % Variables,'PotSelect' )
    END IF

    IF( .NOT. ( ASSOCIATED(PotVarA) .AND. ASSOCIATED(PotVarB) .AND. &
        ASSOCIATED(PotSelect) ) ) CALL Fatal(Caller,&
        'closed coil component needs the CoilSolver with Coil Closed')

    NodeIndexes => Element % NodeIndexes

    ! Same branch as the CoilSolver's LocalFluxMatrix: the cut of CoilPotB is
    ! elsewhere, so it is the smooth branch where PotSelect marks it.
    UseB = .TRUE.
    DO i=1,nn
      j = PotSelect % Perm( NodeIndexes(i) )
      IF( j == 0 ) CALL Fatal(Caller,&
          'closed coil component needs the CoilSolver with Coil Closed')
      IF( PotSelect % Values(j) <= 0.0_dp ) THEN
        UseB = .FALSE.
        EXIT
      END IF
    END DO

    IF( UseB ) THEN
      PotVar => PotVarB
      Mult = ListGetConstReal( CompParams,'Coil Potential Multiplier B', Found )
    ELSE
      PotVar => PotVarA
      Mult = ListGetConstReal( CompParams,'Coil Potential Multiplier', Found )
    END IF

    IF( .NOT. Found .OR. Mult == 0.0_dp ) CALL Fatal(Caller,&
        'closed coil component needs the CoilSolver with Coil Closed')

    ! The cut potential jumps by 1 but does not circulate by 1. The keyword is
    ! absent only while ComputeCoilCirculation is measuring that difference.
    Circ = ListGetConstReal( CompParams,'Coil Circulation', Found )
    IF( Found ) Mult = Mult * Circ

    DO i=1,nn
      j = PotVar % Perm( NodeIndexes(i) )
      IF( j == 0 ) CALL Fatal(Caller,&
          'closed coil component needs the CoilSolver with Coil Closed')
      Wbase(i) = PotVar % Values(j) / Mult
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE GetCoilWBase
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AddComponentsToBodyLists()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    LOGICAL :: Found
    INTEGER :: i, j, k
    LOGICAL, SAVE :: Visited=.FALSE.
    ! Components and Bodies:
    ! ----------------------  
    INTEGER :: BodyId, BoundaryId
    INTEGER, POINTER :: BodyAssociations(:) => Null()
    INTEGER, POINTER :: BCAssociations(:) => Null()
    TYPE(Valuelist_t), POINTER :: BodyParams, BCParams, ComponentParams
     
    IF (Visited) RETURN

    Visited = .TRUE.
    DO i = 1, SIZE(CurrentModel % Components)
      ComponentParams => CurrentModel % Components(i) % Values

      IF( ListGetLogical( ComponentParams,'Passive Component', Found ) ) CYCLE 
      
      IF (.NOT. ASSOCIATED(ComponentParams)) CALL Fatal ('AddComponentsToBodyList', &
                                                         'Component parameters not found!')
      BodyAssociations => ListGetIntegerArray(ComponentParams, 'Body', Found)

      IF (.NOT. Found) BodyAssociations => ListGetIntegerArray(ComponentParams, 'Master Bodies', Found)

      IF (.NOT. Found) BCAssociations => ListGetIntegerArray(ComponentParams, 'Master BCs', Found)
      
      IF (.NOT. Found) CYCLE

      IF (ASSOCIATED(BodyAssociations)) THEN
        DO j = 1, SIZE(BodyAssociations)
          BodyId = BodyAssociations(j)
          BodyParams => CurrentModel % Bodies(BodyId) % Values
          IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('AddComponentsToBodyList', &
                                                        'Body parameters not found!')
          k = GetInteger(BodyParams, 'Component', Found)
          IF (Found) CALL Fatal ('AddComponentsToBodyList', &
                                'Body '//TRIM(i2s(BodyId))//' associated to two components!')
          CALL listAddInteger(BodyParams, 'Component', i)
          BodyParams => Null()
        END DO
      END IF
      IF (ASSOCIATED(BCAssociations)) THEN
        DO j = 1, SIZE(BCAssociations)
          BoundaryId = BCAssociations(j)
          BCParams => CurrentModel % BCs(BoundaryId) % Values
          IF (.NOT. ASSOCIATED(BCParams)) CALL Fatal ('AddComponentsToBodyList', &
                         'Boundary Condition parameters not found!')

          k = GetInteger(BCParams, 'Component', Found)

          IF (Found) CALL Fatal ('AddComponentsToBodyList', &
            'Boundary Condition '//TRIM(i2s(BoundaryId))//' associated to two components!')

          CALL ListAddInteger(BCParams, 'Component', i)
          BCParams => Null()
        END DO
      END IF
      BodyAssociations => Null()
      BCAssociations => Null()
    END DO

    DO i = 1, CurrentModel % NumberOfBodies
      BodyParams => CurrentModel % Bodies(i) % Values
      IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('AddComponentsToBodyList', &
                          'Body parameters not found!')
      j = GetInteger(BodyParams, 'Component', Found)
      IF (.NOT. Found) CYCLE

      WRITE(Message,'(A)') '"Body '//TRIM(I2S(i))//'" associated to "Component '//TRIM(I2S(j))//'"' 
      CALL Info('AddComponentsToBodyList',Message,Level=5)
      BodyParams => Null()
    END DO

    DO i = 1, CurrentModel % NumberOfBCs
      BCParams => CurrentModel % BCs(i) % Values
      IF (.NOT. ASSOCIATED(BCParams)) CALL Fatal ('AddComponentsToBodyList', &
                  'Boundary Condition parameters not found!')
      j = GetInteger(BCParams, 'Component', Found)
      IF (.NOT. Found) CYCLE

      WRITE(Message,'(A)') '"Boundary Condition '//TRIM(I2S(i))// &
          '" associated to "Component '//TRIM(I2S(j))//'"' 
      CALL Info('AddComponentsToBodyList',Message,Level=5)
      BCParams => Null()
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE AddComponentsToBodyLists
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE CheckComponentVariables()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    LOGICAL :: Found
    INTEGER :: i, j, k
    LOGICAL, SAVE :: Visited=.FALSE.
    TYPE(Valuelist_t), POINTER :: ComponentParams
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType, VarName
   
    IF(Visited) RETURN
    IF(CurrentModel % NumberOfComponents == 0) RETURN
    
    Visited = .TRUE.

    j = 0
    DO i = 1, CurrentModel % NumberOfComponents      
      ComponentParams => CurrentModel % Components(i) % Values                  
      IF( ListGetLogical( ComponentParams,'Passive Component', Found ) ) CYCLE 
      CoilType = GetString(ComponentParams, 'Coil Type', Found)
      IF(.NOT. Found) CYCLE

      SELECT CASE (CoilType)
      CASE ('stranded')
        VarName = 'Circuit Current Variable Id'
      CASE ('massive','foil winding','flat wire','foil sheet')
        VarName = 'Circuit Voltage Variable Id'
      CASE DEFAULT
        CYCLE
      END SELECT

      IF(.NOT. ListCheckPresent(ComponentParams, VarName) ) THEN
        CALL Warn('CheckComponentOwners','Coil type given for component '//I2S(i)//' but no: '//TRIM(VarName))
        j = j + 1
      END IF
    END DO

    IF(j > 0) THEN
      CALL Warn('CheckComponentOwners','Could not find variables for '//I2S(j)//' coils!')
      CALL Fatal('CheckComponentOwners','Check your circuit settings!')
    END IF

  END SUBROUTINE CheckComponentVariables
      

  
!------------------------------------------------------------------------------
  FUNCTION GetComponentBodyIds(Id) RESULT (BodyIds)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    LOGICAL :: Found
    INTEGER :: Id
    INTEGER, POINTER :: BodyIds(:)
    TYPE(Valuelist_t), POINTER :: ComponentParams
    
    ComponentParams => CurrentModel % Components(Id) % Values
    
    IF (.NOT. ASSOCIATED(ComponentParams)) CALL Fatal ('GetComponentBodyIds', &
                      'Component parameters not found!')
    BodyIds => ListGetIntegerArray(ComponentParams, 'Body', Found)
    IF (.NOT. Found) BodyIds => ListGetIntegerArray(ComponentParams, 'Master Bodies', Found)
    IF (.NOT. Found) BodyIds => Null()
    
!------------------------------------------------------------------------------
  END FUNCTION GetComponentBodyIds
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION GetComponentHomogenizationBodyIds(Id) RESULT (BodyIds)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    LOGICAL :: Found
    INTEGER :: Id
    INTEGER, POINTER :: BodyIds(:)
    TYPE(Valuelist_t), POINTER :: ComponentParams
    
    ComponentParams => CurrentModel % Components(Id) % Values
    
    IF (.NOT. ASSOCIATED(ComponentParams)) CALL Fatal ('GetComponentHomogenizationBodyIds', &
                          'Component parameters not found!')
    BodyIds => ListGetIntegerArray(ComponentParams, 'Homogenization Parameters Body', Found)
    IF (.NOT. Found) BodyIds => GetComponentBodyIds(Id)

!------------------------------------------------------------------------------
  END FUNCTION GetComponentHomogenizationBodyIds
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Read one transient-homogenization ladder triplet (y0, alpha, Sigma) attached
! to a Component. The SIF keywords are
!     <name> y0          = Real ...
!     <name> alpha       = Real ...
!     <name> Sigma(n,n)  = Real ...
! where <name> is one of "Nu 11", "Nu 22", "Sigma 33". Hard-fails on any
! missing keyword; that's deliberate — silent DC fallbacks mask Python emitter
! bugs.
!------------------------------------------------------------------------------
  SUBROUTINE GetTransientHomogenizationLadder(CompParams, name, n_ladder, y0, alpha, SigmaMat)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: CompParams
    CHARACTER(*), INTENT(IN)   :: name
    INTEGER, INTENT(IN)        :: n_ladder
    REAL(KIND=dp), INTENT(OUT) :: y0, alpha
    REAL(KIND=dp), INTENT(OUT) :: SigmaMat(n_ladder, n_ladder)

    LOGICAL :: Found
    INTEGER :: i, j, DummyNodeIdx(1)
    REAL(KIND=dp), POINTER :: Hwrk(:,:,:)
    CHARACTER(LEN=:), ALLOCATABLE :: key
!------------------------------------------------------------------------------
    IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal('GetTransientHomogenizationLadder', &
        'Component parameters not associated when reading "' // TRIM(name) // '"')

    key = TRIM(name) // ' y0'
    y0 = GetConstReal(CompParams, key, Found)
    IF (.NOT. Found) CALL Fatal('GetTransientHomogenizationLadder', &
        'Missing keyword "' // key // '" on Component (required when Transient Homogenization = True)')

    key = TRIM(name) // ' alpha'
    alpha = GetConstReal(CompParams, key, Found)
    IF (.NOT. Found) CALL Fatal('GetTransientHomogenizationLadder', &
        'Missing keyword "' // key // '" on Component (required when Transient Homogenization = True)')

    key = TRIM(name) // ' Sigma'
    Hwrk => NULL()
    DummyNodeIdx = 1
    CALL ListGetRealArray(CompParams, key, Hwrk, 1, DummyNodeIdx, Found)
    IF (.NOT. Found) CALL Fatal('GetTransientHomogenizationLadder', &
        'Missing keyword "' // key // '(n,n)" on Component (required when Transient Homogenization = True)')

    IF (SIZE(Hwrk, 1) /= n_ladder .OR. SIZE(Hwrk, 2) /= n_ladder) THEN
      CALL Fatal('GetTransientHomogenizationLadder', &
          '"' // key // '" matrix shape does not match Homogenization Ladder Order')
    END IF

    DO i = 1, n_ladder
      DO j = 1, n_ladder
        SigmaMat(i, j) = Hwrk(i, j, 1)
      END DO
    END DO

    IF (ASSOCIATED(Hwrk)) DEALLOCATE(Hwrk)
!------------------------------------------------------------------------------
  END SUBROUTINE GetTransientHomogenizationLadder
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Foster form of a homogenization ladder: y(s) = y0 + sum_k r_k/(1 + s T_k).
!> Reads '<name> Residues(N)' and '<name> Taus(N)'. For N = 1 it also accepts the
!> older triplet '<name> alpha' / '<name> Sigma(1,1)', which is the same thing
!> with r_1 = alpha and T_1 = Sigma(1,1), so existing SIFs keep working bit for
!> bit. '<name> y0' is required in both forms.
!------------------------------------------------------------------------------
  SUBROUTINE GetFosterLadder(CompParams, name, n_ladder, y0, res, taus)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: CompParams
    CHARACTER(*), INTENT(IN)   :: name
    INTEGER, INTENT(IN)        :: n_ladder
    REAL(KIND=dp), INTENT(OUT) :: y0
    REAL(KIND=dp), INTENT(OUT) :: res(n_ladder), taus(n_ladder)

    LOGICAL :: Found, FoundR, FoundT
    INTEGER :: k, DummyNodeIdx(1)
    REAL(KIND=dp), POINTER :: Hwrk(:,:,:)
    REAL(KIND=dp) :: SigmaMat(1,1), alpha
    CHARACTER(LEN=:), ALLOCATABLE :: key
!------------------------------------------------------------------------------
    IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal('GetFosterLadder', &
        'Component parameters not associated when reading "' // TRIM(name) // '"')
    IF (n_ladder < 1 .OR. n_ladder > 6) CALL Fatal('GetFosterLadder', &
        'Homogenization Ladder Order must be between 1 and 6')

    key = TRIM(name) // ' y0'
    y0 = GetConstReal(CompParams, key, Found)
    IF (.NOT. Found) CALL Fatal('GetFosterLadder', 'Missing keyword "' // key // '" on Component')

    res  = 0._dp
    taus = 1._dp

    key = TRIM(name) // ' Residues'
    Hwrk => NULL(); DummyNodeIdx = 1
    CALL ListGetRealArray(CompParams, key, Hwrk, 1, DummyNodeIdx, FoundR)
    IF (FoundR) THEN
      IF (SIZE(Hwrk,1)*SIZE(Hwrk,2) /= n_ladder) CALL Fatal('GetFosterLadder', &
          '"' // key // '" does not have Homogenization Ladder Order entries')
      DO k = 1, n_ladder
        IF (SIZE(Hwrk,1) >= k) THEN
          res(k) = Hwrk(k,1,1)
        ELSE
          res(k) = Hwrk(1,k,1)
        END IF
      END DO
      DEALLOCATE(Hwrk)
    END IF

    key = TRIM(name) // ' Taus'
    Hwrk => NULL()
    CALL ListGetRealArray(CompParams, key, Hwrk, 1, DummyNodeIdx, FoundT)
    IF (FoundT) THEN
      IF (SIZE(Hwrk,1)*SIZE(Hwrk,2) /= n_ladder) CALL Fatal('GetFosterLadder', &
          '"' // key // '" does not have Homogenization Ladder Order entries')
      DO k = 1, n_ladder
        IF (SIZE(Hwrk,1) >= k) THEN
          taus(k) = Hwrk(k,1,1)
        ELSE
          taus(k) = Hwrk(1,k,1)
        END IF
      END DO
      DEALLOCATE(Hwrk)
    END IF

    IF (FoundR .AND. FoundT) RETURN

    IF (n_ladder /= 1) CALL Fatal('GetFosterLadder', &
        'Missing "' // TRIM(name) // ' Residues(N)" / "' // TRIM(name) // &
        ' Taus(N)"; the alpha/Sigma triplet is only accepted for order 1')

    ! Order 1 fallback: the original triplet.
    CALL GetTransientHomogenizationLadder(CompParams, name, 1, y0, alpha, SigmaMat)
    res(1)  = alpha
    taus(1) = SigmaMat(1,1)
!------------------------------------------------------------------------------
  END SUBROUTINE GetFosterLadder
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Foster coefficients of the foil sheet in-plane reluctivity, computed from the
!> physics instead of read from the SIF, so that harmonic and transient share one
!> source of truth.
!>
!>   mu_e(s)/mu0 = (1-ff) + ff [ sum_{n<=N} c_n/(1+s tau_n) + tail ],
!>   c_n = 2/((n-1/2)pi)^2,  tau_n = tau0/((n-1/2)pi)^2,  tail = 1 - sum c_n.
!>
!> The truncation remainder is folded into the constant so that mu_e is exact at
!> DC; the price is that the s -> infinity "fully excluded" limit becomes
!> approximate, which is outside the band of interest, exactly as for the tail
!> inductance of the strand skin ladder. nu_e = 1/mu_e is then rational of
!> degree N: its poles are the roots of the numerator of mu_e and its residues
!> follow analytically, r_k = T_k D(s_k)/(mu0 Nm'(s_k)).
!------------------------------------------------------------------------------
  SUBROUTINE FoilSheetNuFoster(tau0, ff, n_ladder, nu_inf, res, taus)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(IN)  :: tau0, ff
    INTEGER, INTENT(IN)        :: n_ladder
    REAL(KIND=dp), INTENT(OUT) :: nu_inf, res(n_ladder), taus(n_ladder)

    INTEGER, PARAMETER :: MAXN = 6
    REAL(KIND=dp) :: c(MAXN), tau(MAXN), den(0:MAXN), num(0:MAXN), prod(0:MAXN)
    REAL(KIND=dp) :: cmat(MAXN,MAXN), wr(MAXN), wi(MAXN), vdum(1,1), work(8*MAXN)
    REAL(KIND=dp) :: mu0, cst, sk, dval, dnum, leadD, leadN
    INTEGER :: n, k, i, j, info

    IF (n_ladder < 1 .OR. n_ladder > MAXN) &
        CALL Fatal('FoilSheetNuFoster', &
            'Foil sheet: "Homogenization Ladder Order" must be between 1 and 6!')
    mu0 = 4.0d-7 * PI

    DO n = 1, n_ladder
      c(n)   = 2._dp / (((REAL(n,dp) - 0.5_dp)*PI)**2)
      tau(n) = tau0  / (((REAL(n,dp) - 0.5_dp)*PI)**2)
    END DO
    cst = (1._dp - ff) + ff * (1._dp - SUM(c(1:n_ladder)))

    ! den(s) = prod (1 + s tau_n), ascending coefficients.
    den = 0._dp; den(0) = 1._dp
    DO n = 1, n_ladder
      DO i = n, 1, -1
        den(i) = den(i) + tau(n)*den(i-1)
      END DO
    END DO

    ! num(s) = cst*den + ff * sum_n c_n * prod_{m/=n} (1 + s tau_m)
    num = cst * den
    DO n = 1, n_ladder
      prod = 0._dp; prod(0) = 1._dp
      j = 0
      DO i = 1, n_ladder
        IF (i == n) CYCLE
        j = j + 1
        DO k = j, 1, -1
          prod(k) = prod(k) + tau(i)*prod(k-1)
        END DO
      END DO
      num(0:n_ladder) = num(0:n_ladder) + ff*c(n)*prod(0:n_ladder)
    END DO

    leadD = den(n_ladder)
    leadN = num(n_ladder)
    IF (ABS(leadN) <= 0._dp) CALL Fatal('FoilSheetNuFoster', &
        'Foil sheet: degenerate reluctivity polynomial, lower "Homogenization Ladder Order"!')
    nu_inf = (leadD/leadN)/mu0

    ! Roots of num via the companion matrix of its monic form.
    cmat = 0._dp
    DO j = 1, n_ladder
      cmat(1,j) = -num(n_ladder-j)/leadN
    END DO
    DO i = 2, n_ladder
      cmat(i,i-1) = 1._dp
    END DO
    CALL DGEEV('N','N', n_ladder, cmat, MAXN, wr, wi, vdum, 1, vdum, 1, work, 8*MAXN, info)
    IF (info /= 0) CALL Fatal('FoilSheetNuFoster', &
        'Foil sheet: cannot factor the reluctivity polynomial, lower "Homogenization Ladder Order"!')

    DO k = 1, n_ladder
      IF (ABS(wi(k)) > 1.0d-8*MAX(ABS(wr(k)),1._dp)) &
          CALL Fatal('FoilSheetNuFoster', &
              'Foil sheet: complex pole in the reluctivity ladder, lower "Homogenization Ladder Order"!')
      sk = wr(k)
      IF (sk >= 0._dp) CALL Fatal('FoilSheetNuFoster', &
          'Foil sheet: unstable pole in the reluctivity ladder, lower "Homogenization Ladder Order"!')
      taus(k) = -1._dp/sk
      dval = 0._dp
      DO i = n_ladder, 0, -1
        dval = dval*sk + den(i)
      END DO
      dnum = 0._dp
      DO i = n_ladder, 1, -1
        dnum = dnum*sk + REAL(i,dp)*num(i)
      END DO
      res(k) = taus(k) * dval / (mu0 * dnum)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE FoilSheetNuFoster
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION FindSolverWithKey(key) RESULT (Solver)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    
    CHARACTER(*) :: key

    LOGICAL :: Found
    INTEGER :: i
    TYPE(Solver_t), POINTER :: Solver
    
    ! Look for the solver we attach the circuit equations to:
    ! -------------------------------------------------------
    Found = .FALSE.
    DO i=1, CurrentModel % NumberOfSolvers
      Solver => CurrentModel % Solvers(i)
      IF(ListCheckPresent(Solver % Values, key)) THEN 
        Found = .TRUE. 
        EXIT
      END IF
    END DO
    
    IF (.NOT. Found) CALL Fatal('FindSolverWithKey', & 
       TRIM(Key)//' keyword not found in any of the solvers!')

!------------------------------------------------------------------------------
  END FUNCTION FindSolverWithKey
!------------------------------------------------------------------------------




  
!------------------------------------------------------------------------------
!> Nodal values of the stacking coordinate and the across coordinate of a
!> flat wire element, taken from the Alpha and Beta direction fields.
!------------------------------------------------------------------------------
  SUBROUTINE GetFlatWireLocalFields(StackAlongAlpha, Element, n, sStack, sAcross)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    LOGICAL :: StackAlongAlpha
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n
    REAL(KIND=dp) :: sStack(:), sAcross(:)
    REAL(KIND=dp) :: alpha(n), beta(n)

    CALL GetLocalSolution(alpha, 'Alpha', UElement=Element)
    CALL GetLocalSolution(beta, 'Beta', UElement=Element)
    IF (StackAlongAlpha) THEN
      sStack(1:n) = alpha
      sAcross(1:n) = beta
    ELSE
      sStack(1:n) = beta
      sAcross(1:n) = alpha
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE GetFlatWireLocalFields
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> 1-based turn cell index of a point with normalized coordinates sStack and
!> sAcross in [0,1]. Cells are numbered across first, then along the stack.
!------------------------------------------------------------------------------
  FUNCTION FlatWireCellIndex(nStack, nAcross, sStack, sAcross) RESULT(k)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nStack, nAcross, k
    REAL(KIND=dp) :: sStack, sAcross
    INTEGER :: ks, ka

    ks = MIN(nStack,  MAX(1, FLOOR(sStack  * nStack)  + 1))
    ka = MIN(nAcross, MAX(1, FLOOR(sAcross * nAcross) + 1))
    k = (ks-1) * nAcross + ka
!------------------------------------------------------------------------------
  END FUNCTION FlatWireCellIndex
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Foil sheet: 1-based cell index (along the stacking direction) and segment
!> index (across it) of a point with normalized coordinates in [0,1].
!------------------------------------------------------------------------------
  SUBROUTINE FoilSheetStrand(nCells, nSegments, sStack, sAcross, k, j)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nCells, nSegments, k, j
    REAL(KIND=dp) :: sStack, sAcross

    k = MIN(nCells,    MAX(1, FLOOR(sStack  * nCells)    + 1))
    j = MIN(nSegments, MAX(1, FLOOR(sAcross * nSegments) + 1))
!------------------------------------------------------------------------------
  END SUBROUTINE FoilSheetStrand
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Local reluctivity tensor directions of a foil sheet. Local direction 1 is
!> Alpha, 2 is Beta and 3 is Gamma, the current direction, which lies in the
!> turn plane whichever way the turns are stacked. The stacking normal keeps
!> 1/mu0 and the two in-plane directions carry the complex stack permeability,
!> so 'Stacking Direction' alone decides which tensor component is which:
!> alpha stacking (foils, flatwise flat wire) homogenizes Beta and Gamma,
!> beta stacking (edgewise flat wire) homogenizes Alpha and Gamma.
!------------------------------------------------------------------------------
  SUBROUTINE FoilSheetNuDirections(StackAlongAlpha, dStack, dPlane)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    LOGICAL :: StackAlongAlpha
    INTEGER :: dStack, dPlane(2)

    IF (StackAlongAlpha) THEN
      dStack = 1
      dPlane = [2, 3]
    ELSE
      dStack = 2
      dPlane = [1, 3]
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE FoilSheetNuDirections
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Reluctivity keyword of a local direction: 'Nu 11', 'Nu 22' or 'Nu 33'.
!------------------------------------------------------------------------------
  FUNCTION FoilSheetNuKey(d) RESULT(key)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: d
    CHARACTER(LEN=5) :: key

    WRITE(key,'(A,I1,I1)') 'Nu ', d, d
!------------------------------------------------------------------------------
  END FUNCTION FoilSheetNuKey
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Offset of the strand (l,j) current dof inside the foil sheet voltage
!> variable, l being the sub-layer index along the stack. Layout: 0 = V,
!> 1..nCells = V_k, then the strands sub-layer by sub-layer. With one sub-layer
!> per turn the sub-layer index is the cell index and this is the old layout.
!------------------------------------------------------------------------------
  FUNCTION FoilSheetStrandDof(nCells, nSegments, l, j) RESULT(ind)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nCells, nSegments, l, j, ind

    ind = nCells + (l-1) * nSegments + j
!------------------------------------------------------------------------------
  END FUNCTION FoilSheetStrandDof
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Stacking coordinate edges of strand band l of a foil sheet. With one band per
!> turn the band fills the pitch and the current is smeared over it, which is
!> the validated model. With nSublayers > 1 the bands are the copper itself: the
!> m of them sit in the central 'Fill Factor' of the pitch, t/m apart, and the
!> two margins of (1-ff)/2 carry no strand. That spacing is the whole point -
!> the loop area between two sub-layers of a turn is what drives the intra-turn
!> circulating current, and over the pitch it would be 1/ff too large.
!------------------------------------------------------------------------------
  SUBROUTINE FoilSheetBand(nCells, nSublayers, ff, l, s1, s2)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nCells, nSublayers, l
    REAL(KIND=dp) :: ff, s1, s2
    INTEGER :: k, ls
    REAL(KIND=dp) :: pitch

    pitch = 1._dp / nCells
    IF (nSublayers == 1) THEN
      s1 = (l-1) * pitch
      s2 = l * pitch
    ELSE
      k  = (l-1) / nSublayers + 1
      ls = l - (k-1) * nSublayers
      s1 = (k-1) * pitch + 0.5_dp * (1._dp - ff) * pitch + (ls-1) * ff * pitch / nSublayers
      s2 = s1 + ff * pitch / nSublayers
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE FoilSheetBand
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Strand band of a stacking coordinate, 0 when the point is in the insulation
!> margin of its turn and carries no strand.
!------------------------------------------------------------------------------
  FUNCTION FoilSheetBandIndex(nCells, nSublayers, ff, s) RESULT(l)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nCells, nSublayers, l
    REAL(KIND=dp) :: ff, s
    INTEGER :: k, ls
    REAL(KIND=dp) :: u, marg

    k = MIN(nCells, MAX(1, FLOOR(s * nCells) + 1))
    IF (nSublayers == 1) THEN
      l = k
      RETURN
    END IF
    u = s * nCells - (k-1)
    marg = 0.5_dp * (1._dp - ff)
    IF (u < marg .OR. u > 1._dp - marg) THEN
      l = 0
      RETURN
    END IF
    ls = MIN(nSublayers, MAX(1, FLOOR((u - marg) * nSublayers / ff) + 1))
    l = (k-1) * nSublayers + ls
!------------------------------------------------------------------------------
  END FUNCTION FoilSheetBandIndex
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Turn cell of a sub-layer: sub-layers (k-1)*nSublayers+1 .. k*nSublayers are
!> the same conductor in parallel and share the voltage dof of turn cell k.
!------------------------------------------------------------------------------
  FUNCTION FoilSheetLayerCell(nSublayers, l) RESULT(k)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nSublayers, l, k

    k = (l-1) / nSublayers + 1
!------------------------------------------------------------------------------
  END FUNCTION FoilSheetLayerCell
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Euler potential form of the foil sheet strand direction at an integration
!> point: t = s * grad(stack) x grad(across), with s = +-1 fixed once per
!> component so that t points along grad(W). Swapping the two fields only flips
!> the cross product, which the sign absorbs, so the strand direction is the
!> same physical field for either stacking direction.
!>
!> This field is exactly solenoidal in the discrete weak sense - it is a cross
!> product of two gradients, so div t vanishes identically and its normal
!> component is continuous across element faces, because it only involves the
!> tangential gradients of the two continuous nodal fields. It is also exactly
!> tangential to every iso-surface of Alpha and of Beta, which are precisely the
!> strand interfaces, so a strand-wise constant current density c_kj * t carries
!> no current across them: div(c t) = grad(c) . t = 0 because grad(c) is a
!> combination of grad(Alpha) and grad(Beta), both perpendicular to t. That is
!> what makes the piecewise constant strand source consistent, provided the
!> strand index is taken at the integration point (FoilSheetStrand on the
!> interpolated Alpha and Beta), not once per element.
!------------------------------------------------------------------------------
  FUNCTION FoilSheetDirection(sStack, sAcross, dBasisdx, n, sgn) RESULT(t)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: sStack(:), sAcross(:), dBasisdx(:,:), sgn, t(3)
    INTEGER :: n
    REAL(KIND=dp) :: ga(3), gb(3)

    ga = MATMUL(sStack(1:n), dBasisdx(1:n,:))
    gb = MATMUL(sAcross(1:n), dBasisdx(1:n,:))
    t(1) = ga(2)*gb(3) - ga(3)*gb(2)
    t(2) = ga(3)*gb(1) - ga(1)*gb(3)
    t(3) = ga(1)*gb(2) - ga(2)*gb(1)
    t = sgn * t
!------------------------------------------------------------------------------
  END FUNCTION FoilSheetDirection
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> A foil sheet strand with no integration point of its own would leave a zero
!> row in the circuit matrix. Fail early and say which (cell, segment) is empty.
!------------------------------------------------------------------------------
  SUBROUTINE CheckFoilSheetStrands(Comp)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Component_t), POINTER :: Comp
    REAL(KIND=dp) :: w, wmax, bs1, bs2, bandflux
    INTEGER :: k, j, ind, nempty, nLayers

    IF (.NOT. ALLOCATED(Comp % StrandWeight)) RETURN

    wmax = 0._dp
    DO ind = 1, SIZE(Comp % StrandWeight)
      Comp % StrandWeight(ind) = ParallelReduction(Comp % StrandWeight(ind))
      wmax = MAX(wmax, Comp % StrandWeight(ind))
    END DO
    IF (wmax <= 0._dp) RETURN

    nempty = 0
    nLayers = Comp % nCells * Comp % nSublayers
    DO k = 1, nLayers
      DO j = 1, Comp % nSegments
        ind = (k-1) * Comp % nSegments + j
        w = Comp % StrandWeight(ind)
        IF (w > 1.0d-8 * wmax) CYCLE
        nempty = nempty + 1
        IF (nempty <= 20) CALL Error('CheckFoilSheetStrands', &
            'Foil sheet strand (sub-layer '//I2S(k)//', segment '//I2S(j)//') has no element!')
      END DO
    END DO

    IF (nempty > 0) THEN
      CALL Error('CheckFoilSheetStrands','Component '//I2S(Comp % ComponentId)//': '// &
          I2S(nempty)//' of '//I2S(nLayers * Comp % nSegments)//' strands are empty.')
      CALL Fatal('CheckFoilSheetStrands', &
          'Lower "Sheet Cells" / "Sheet Segments" / "Sheet Sublayers" or refine the coil mesh!')
    END IF

    ! With the Euler potential direction the strand flux is exactly
    ! Delta(alpha_k)*Delta(beta_j) = 1/(nCells*nSegments), so this ratio is 1 up
    ! to the quadrature error of the elements cut by a strand interface. It is
    ! the sharpest check that the strand bookkeeping and the direction sign are
    ! right, and it is mesh independent.
    ! The flux of a band is its own stacking width times its segment width, so
    ! this ratio is 1 up to quadrature error whatever the layout. With copper
    ! only bands the width carries the fill factor, which is why it is taken
    ! from FoilSheetBand rather than assumed uniform.
    CALL FoilSheetBand(Comp % nCells, Comp % nSublayers, Comp % FillFactor, 1, bs1, bs2)
    bandflux = (bs2 - bs1) / Comp % nSegments
    WRITE(Message,'(A,F14.10,A,F14.10)') 'Foil sheet strand flux / (dStack dAcross): min ', &
        MINVAL(Comp % StrandWeight) / bandflux, &
        '  max ', MAXVAL(Comp % StrandWeight) / bandflux
    CALL Info('CheckFoilSheetStrands', Message, Level=5)
!------------------------------------------------------------------------------
  END SUBROUTINE CheckFoilSheetStrands
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> tanh(u)/u for complex u, the one transcendental the foil sheet homogenization
!> needs: the sheet conductivity is ff*sigma*tanh(u)/u and the in-plane
!> permeability is mu0[(1-ff) + ff tanh(u)/u], with u = (1+i) t/(2 delta).
!> The series keeps it accurate as u -> 0, where tanh(u)/u is 0/0.
!------------------------------------------------------------------------------
  FUNCTION FoilSheetTanhOverU(u) RESULT(v)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX(KIND=dp) :: u, v

    IF (ABS(u) < 1.0d-4) THEN
      v = 1._dp - u*u/3._dp
    ELSE
      v = TANH(u)/u
    END IF
!------------------------------------------------------------------------------
  END FUNCTION FoilSheetTanhOverU
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Volume of a sub tetrahedron relative to its parent. The sub tet is given by
!> the barycentric coordinates of its four vertices in the parent, t(bary,vertex),
!> so the parent itself has the identity and volume 1.
!------------------------------------------------------------------------------
  FUNCTION TetVolumeFraction(t) RESULT(vol)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: t(4,4), vol
    REAL(KIND=dp) :: d1(3), d2(3), d3(3)

    d1 = t(2:4,2) - t(2:4,1)
    d2 = t(2:4,3) - t(2:4,1)
    d3 = t(2:4,4) - t(2:4,1)
    vol = ABS( d1(1)*(d2(2)*d3(3) - d2(3)*d3(2)) &
             - d1(2)*(d2(1)*d3(3) - d2(3)*d3(1)) &
             + d1(3)*(d2(1)*d3(2) - d2(2)*d3(1)) )
!------------------------------------------------------------------------------
  END FUNCTION TetVolumeFraction
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Clip a list of sub tetrahedra of a common parent by the half space f >= 0,
!> where f is affine on the parent and given by its four nodal values fNod.
!> A tet cut by a plane leaves a tet (one vertex inside), a wedge (three inside)
!> or a wedge (two inside); a wedge splits into three tets, so the result is
!> again a list of tets and the recursion over several planes stays closed.
!> Vertices are carried as barycentric coordinates of the parent, which keeps
!> the cut exact and lets every integrand be evaluated from the parent basis.
!> A plane through a vertex or an edge produces zero volume pieces, which are
!> kept here and dropped by the caller on their volume.
!------------------------------------------------------------------------------
  SUBROUTINE ClipTetsByHalfSpace(tets, ntet, fNod, MaxTet)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: ntet, MaxTet
    REAL(KIND=dp) :: tets(4,4,MaxTet), fNod(4)

    REAL(KIND=dp) :: out(4,4,MaxTet), fv(4), tt, a(4,3), b(4,3)
    INTEGER :: i, v, nout, nin, ins(4), outs(4), nio, noo, p, q, r, s

    nout = 0
    DO i = 1, ntet
      DO v = 1, 4
        fv(v) = SUM( tets(1:4,v,i) * fNod(1:4) )
      END DO

      nio = 0; noo = 0
      DO v = 1, 4
        IF (fv(v) >= 0._dp) THEN
          nio = nio + 1; ins(nio) = v
        ELSE
          noo = noo + 1; outs(noo) = v
        END IF
      END DO
      nin = nio

      SELECT CASE(nin)

      CASE(0)
        CYCLE

      CASE(4)
        nout = nout + 1
        IF (nout > MaxTet) CALL Fatal('ClipTetsByHalfSpace','Too many pieces!')
        out(:,:,nout) = tets(:,:,i)

      CASE(1)
        ! One vertex inside: a single corner tet.
        p = ins(1)
        nout = nout + 1
        IF (nout > MaxTet) CALL Fatal('ClipTetsByHalfSpace','Too many pieces!')
        out(:,1,nout) = tets(:,p,i)
        DO v = 1, 3
          q = outs(v)
          tt = fv(p) / (fv(p) - fv(q))
          out(:,v+1,nout) = (1._dp - tt) * tets(:,p,i) + tt * tets(:,q,i)
        END DO

      CASE(3)
        ! Three vertices inside: the tet minus the corner at the outside vertex,
        ! i.e. a wedge between the inside triangle and its cut image.
        s = outs(1)
        DO v = 1, 3
          p = ins(v)
          a(:,v) = tets(:,p,i)
          tt = fv(p) / (fv(p) - fv(s))
          b(:,v) = (1._dp - tt) * tets(:,p,i) + tt * tets(:,s,i)
        END DO
        CALL EmitWedge(out, nout, MaxTet, a, b)

      CASE(2)
        ! Two vertices inside: a wedge whose triangular ends sit at the two
        ! inside vertices.
        p = ins(1); q = ins(2)
        r = outs(1); s = outs(2)
        a(:,1) = tets(:,p,i)
        tt = fv(p) / (fv(p) - fv(r))
        a(:,2) = (1._dp - tt) * tets(:,p,i) + tt * tets(:,r,i)
        tt = fv(p) / (fv(p) - fv(s))
        a(:,3) = (1._dp - tt) * tets(:,p,i) + tt * tets(:,s,i)
        b(:,1) = tets(:,q,i)
        tt = fv(q) / (fv(q) - fv(r))
        b(:,2) = (1._dp - tt) * tets(:,q,i) + tt * tets(:,r,i)
        tt = fv(q) / (fv(q) - fv(s))
        b(:,3) = (1._dp - tt) * tets(:,q,i) + tt * tets(:,s,i)
        CALL EmitWedge(out, nout, MaxTet, a, b)

      END SELECT
    END DO

    ntet = nout
    IF (nout > 0) tets(:,:,1:nout) = out(:,:,1:nout)
!------------------------------------------------------------------------------
  END SUBROUTINE ClipTetsByHalfSpace
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Split a triangular wedge with ends a(:,1:3) and b(:,1:3), vertex v of one end
!> matching vertex v of the other, into three tetrahedra.
!------------------------------------------------------------------------------
  SUBROUTINE EmitWedge(out, nout, MaxTet, a, b)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nout, MaxTet
    REAL(KIND=dp) :: out(4,4,MaxTet), a(4,3), b(4,3)

    IF (nout + 3 > MaxTet) CALL Fatal('EmitWedge','Too many pieces!')

    nout = nout + 1
    out(:,1,nout) = a(:,1); out(:,2,nout) = a(:,2)
    out(:,3,nout) = a(:,3); out(:,4,nout) = b(:,3)

    nout = nout + 1
    out(:,1,nout) = a(:,1); out(:,2,nout) = a(:,2)
    out(:,3,nout) = b(:,3); out(:,4,nout) = b(:,2)

    nout = nout + 1
    out(:,1,nout) = a(:,1); out(:,2,nout) = b(:,2)
    out(:,3,nout) = b(:,3); out(:,4,nout) = b(:,1)
!------------------------------------------------------------------------------
  END SUBROUTINE EmitWedge
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Exact decomposition of a linear tetrahedron into its foil sheet strand
!> pieces. The stacking and across fields are linear inside the tet, so the part
!> of the tet that belongs to strand (k,j) is the convex polytope cut out by
!>    a_k <= stack <= a_(k+1),  b_j <= across <= b_(j+1),
!> obtained here by clipping against the at most four planes that cross the
!> element. One piece is returned per strand: its volume fraction of the parent
!> and the barycentric coordinates of its centroid. Every strand integrand of
!> the foil sheet kernel is affine on the parent tet - t is constant, the nodal
!> gradients are constant, the Whitney edge basis is linear - so the value at
!> the centroid times the volume integrates it exactly, and no quadrature error
!> is left in the strand volume fractions.
!>
!> The outermost cells are open ended on purpose, matching the clamping in
!> FoilSheetStrand, so nodes that fall slightly outside [0,1] still belong to
!> the first or last strand instead of being dropped.
!------------------------------------------------------------------------------
  SUBROUTINE FoilSheetPieces(nCells, nSegments, aNod, bNod, nPiece, pCell, pSeg, pVol, pBary, &
      nSublayers, ff)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nCells, nSegments, nPiece, pCell(:), pSeg(:)
    REAL(KIND=dp) :: aNod(4), bNod(4), pVol(:), pBary(:,:)
    INTEGER, OPTIONAL :: nSublayers
    REAL(KIND=dp), OPTIONAL :: ff

    INTEGER, PARAMETER :: MaxTet = 96
    REAL(KIND=dp) :: tets(4,4,MaxTet), da, db, vol, vtot, ctot(4), plane(4)
    REAL(KIND=dp) :: fillf, s1, s2
    INTEGER :: k, j, l, ls, msub, kmin, kmax, jmin, jmax, i, v, ntet

    msub = 1
    IF (PRESENT(nSublayers)) msub = nSublayers
    fillf = 1._dp
    IF (PRESENT(ff)) fillf = ff

    da = 1._dp / nCells
    db = 1._dp / nSegments

    CALL FoilSheetStrand(nCells, nSegments, MINVAL(aNod), MINVAL(bNod), kmin, jmin)
    CALL FoilSheetStrand(nCells, nSegments, MAXVAL(aNod), MAXVAL(bNod), kmax, jmax)

    nPiece = 0
    IF (msub > 1) THEN
      ! Copper only bands. Both band faces are hard planes - the insulation
      ! margins of a turn belong to no strand and are simply not emitted, so
      ! the pieces cover the fill factor of the element, not all of it.
      DO k = kmin, kmax
        DO ls = 1, msub
          l = (k-1) * msub + ls
          CALL FoilSheetBand(nCells, msub, fillf, l, s1, s2)
          DO j = jmin, jmax

            ntet = 1
            tets(:,:,1) = 0._dp
            DO v = 1, 4
              tets(v,v,1) = 1._dp
            END DO

            plane = aNod - s1
            CALL ClipTetsByHalfSpace(tets, ntet, plane, MaxTet)
            IF (ntet > 0) THEN
              plane = s2 - aNod
              CALL ClipTetsByHalfSpace(tets, ntet, plane, MaxTet)
            END IF
            IF (ntet > 0 .AND. j > 1) THEN
              plane = bNod - REAL(j-1,dp) * db
              CALL ClipTetsByHalfSpace(tets, ntet, plane, MaxTet)
            END IF
            IF (ntet > 0 .AND. j < nSegments) THEN
              plane = REAL(j,dp) * db - bNod
              CALL ClipTetsByHalfSpace(tets, ntet, plane, MaxTet)
            END IF
            IF (ntet <= 0) CYCLE

            vtot = 0._dp
            ctot = 0._dp
            DO i = 1, ntet
              vol = TetVolumeFraction(tets(:,:,i))
              IF (vol <= 0._dp) CYCLE
              vtot = vtot + vol
              DO v = 1, 4
                ctot = ctot + 0.25_dp * vol * tets(:,v,i)
              END DO
            END DO
            IF (vtot <= 1.0d-14) CYCLE

            IF (nPiece >= SIZE(pCell)) THEN
              WRITE(Message,'(A,I0,A,I0,A,I0,A,I0,A)') 'An element straddles more than ', &
                  SIZE(pCell), ' strands of the ', nCells, ' x ', msub, ' x ', nSegments, ' layout'
              CALL Error('FoilSheetPieces', Message)
              CALL Fatal('FoilSheetPieces', &
                  'Lower "Sheet Sublayers" / "Sheet Segments" or refine the coil mesh!')
            END IF
            nPiece = nPiece + 1
            pCell(nPiece) = l
            pSeg(nPiece)  = j
            pVol(nPiece)  = vtot
            pBary(:,nPiece) = ctot / vtot
          END DO
        END DO
      END DO
      RETURN
    END IF

    DO k = kmin, kmax
      DO j = jmin, jmax

        ntet = 1
        tets(:,:,1) = 0._dp
        DO v = 1, 4
          tets(v,v,1) = 1._dp
        END DO

        IF (k > 1) THEN
          plane = aNod - REAL(k-1,dp) * da
          CALL ClipTetsByHalfSpace(tets, ntet, plane, MaxTet)
        END IF
        IF (ntet > 0 .AND. k < nCells) THEN
          plane = REAL(k,dp) * da - aNod
          CALL ClipTetsByHalfSpace(tets, ntet, plane, MaxTet)
        END IF
        IF (ntet > 0 .AND. j > 1) THEN
          plane = bNod - REAL(j-1,dp) * db
          CALL ClipTetsByHalfSpace(tets, ntet, plane, MaxTet)
        END IF
        IF (ntet > 0 .AND. j < nSegments) THEN
          plane = REAL(j,dp) * db - bNod
          CALL ClipTetsByHalfSpace(tets, ntet, plane, MaxTet)
        END IF
        IF (ntet <= 0) CYCLE

        vtot = 0._dp
        ctot = 0._dp
        DO i = 1, ntet
          vol = TetVolumeFraction(tets(:,:,i))
          IF (vol <= 0._dp) CYCLE
          vtot = vtot + vol
          DO v = 1, 4
            ctot = ctot + 0.25_dp * vol * tets(:,v,i)
          END DO
        END DO
        IF (vtot <= 1.0d-14) CYCLE

        ! One element can straddle as many strands as the layout is fine, and
        ! sub-layers multiply that. Say so instead of running off the caller's
        ! buffer.
        IF (nPiece >= SIZE(pCell)) THEN
          WRITE(Message,'(A,I0,A,I0,A,I0,A)') 'An element straddles more than ', SIZE(pCell), &
              ' strands of the ', nCells, ' x ', nSegments, ' layout'
          CALL Error('FoilSheetPieces', Message)
          CALL Fatal('FoilSheetPieces', &
              'Lower "Sheet Sublayers" / "Sheet Segments" or refine the coil mesh!')
        END IF
        nPiece = nPiece + 1
        pCell(nPiece) = k
        pSeg(nPiece)  = j
        pVol(nPiece)  = vtot
        pBary(:,nPiece) = ctot / vtot
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE FoilSheetPieces
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Set up the strand quadrature of one element. By default the element is cut
!> exactly along the strand interfaces and each piece is integrated at its
!> centroid, which is exact because every strand integrand is affine on a linear
!> tet. 'Sheet Integration Points' > 0 selects an ordinary Gauss rule with the
!> strand picked at the point instead: the fallback for coil meshes that are not
!> linear tetrahedra. It is not a substitute for the clipping, because on a mesh
!> whose elements are as large as the strands it mis-estimates the strand volume
!> fractions by percents and leaves the assembled source that far from
!> solenoidal.
!------------------------------------------------------------------------------
  SUBROUTINE FoilSheetQuadrature(Comp, CompParams, Element, nn, sStack, sAcross, Q)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Component_t) :: Comp
    TYPE(ValueList_t), POINTER :: CompParams
    TYPE(Element_t), POINTER :: Element
    INTEGER :: nn
    REAL(KIND=dp) :: sStack(:), sAcross(:)
    TYPE(FoilSheetQuad_t) :: Q

    INTEGER :: ngp
    LOGICAL :: Found

    ngp = GetInteger(CompParams, 'Sheet Integration Points', Found)
    IF (.NOT. Found) ngp = 0

    Q % Exact = (ngp <= 0)
    IF (Q % Exact) THEN
      IF (nn /= 4 .OR. Element % TYPE % ElementCode /= 504) CALL Fatal('FoilSheetQuadrature', &
          'Exact strand clipping needs linear tetrahedra; set "Sheet Integration Points" '// &
          'to integrate this coil mesh with a Gauss rule instead!')
      CALL FoilSheetPieces(Comp % nCells, Comp % nSegments, sStack(1:4), sAcross(1:4), &
          Q % nItem, Q % pCell, Q % pSeg, Q % pVol, Q % pBary, Comp % nSublayers, Comp % FillFactor)
    ELSE
      Q % IP = GaussPoints(Element, np=ngp)
      Q % nItem = Q % IP % n
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE FoilSheetQuadrature
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Local coordinates of quadrature item t: the centroid of the clipped piece or
!> the Gauss point of the fallback rule.
!------------------------------------------------------------------------------
  SUBROUTINE FoilSheetQuadPoint(Q, t, u, v, w)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(FoilSheetQuad_t) :: Q
    INTEGER :: t
    REAL(KIND=dp) :: u, v, w

    IF (Q % Exact) THEN
      u = Q % pBary(2,t); v = Q % pBary(3,t); w = Q % pBary(4,t)
    ELSE
      u = Q % IP % U(t); v = Q % IP % V(t); w = Q % IP % W(t)
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE FoilSheetQuadPoint
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Integration weight of quadrature item t, given the Jacobian at its point.
!------------------------------------------------------------------------------
  FUNCTION FoilSheetQuadWeight(Q, t, detJ) RESULT(wgt)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(FoilSheetQuad_t) :: Q
    INTEGER :: t
    REAL(KIND=dp) :: detJ, wgt

    IF (Q % Exact) THEN
      wgt = Q % pVol(t) * detJ / 6._dp
    ELSE
      wgt = Q % IP % s(t) * detJ
    END IF
!------------------------------------------------------------------------------
  END FUNCTION FoilSheetQuadWeight
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Strand of quadrature item t: the band index along the stack and the segment
!> across it. A band index of zero means the point lies in the insulation margin
!> of a turn and belongs to no strand.
!------------------------------------------------------------------------------
  SUBROUTINE FoilSheetQuadStrand(Q, Comp, t, sStack, sAcross, Basis, nn, kc, js)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(FoilSheetQuad_t) :: Q
    TYPE(Component_t) :: Comp
    INTEGER :: t, nn, kc, js
    REAL(KIND=dp) :: sStack(:), sAcross(:), Basis(:)

    INTEGER :: idummy

    js = 0
    IF (Q % Exact) THEN
      kc = Q % pCell(t); js = Q % pSeg(t)
    ELSE
      kc = FoilSheetBandIndex(Comp % nCells, Comp % nSublayers, Comp % FillFactor, &
          SUM(sStack(1:nn)*Basis(1:nn)))
      IF (kc <= 0) RETURN
      CALL FoilSheetStrand(1, Comp % nSegments, 0._dp, SUM(sAcross(1:nn)*Basis(1:nn)), idummy, js)
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE FoilSheetQuadStrand
!------------------------------------------------------------------------------


END MODULE CircuitUtils


MODULE CircuitsMod

  USE DefUtils
  IMPLICIT NONE

CONTAINS 

!------------------------------------------------------------------------------
  SUBROUTINE AllocateCircuitsList()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: slen,n_Circuits
    CHARACTER(:), ALLOCATABLE :: cmd
    CHARACTER(LEN=MAX_NAME_LEN) :: name

    ! Read Circuit definitions from MATC:
    ! ----------------------------------
    n_Circuits = NINT(GetMatcReal("Circuits"))
    CurrentModel % n_Circuits = n_Circuits

    IF( ASSOCIATED( CurrentModel % Circuits ) ) THEN
      IF( SIZE( CurrentModel % Circuits ) == n_Circuits ) THEN
        CALL Info('AllocateCircuitList','Circuit list already allocated!')
      ELSE
        CALL Warn('AllocateCircuitList','Circuit of wrong size already allocated, deallocating this!')
      END IF
    END IF

    IF(.NOT. ASSOCIATED(CurrentModel % Circuits ) ) THEN
      ALLOCATE( CurrentModel % Circuits(n_Circuits) )
    END IF
      
!------------------------------------------------------------------------------
  END SUBROUTINE AllocateCircuitsList
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION CountNofCircVarsOfType(CId, Var_type) RESULT (nofc)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nofc, char_len, slen, CId, i
    CHARACTER(LEN=*) :: Var_type
    CHARACTER(:), ALLOCATABLE :: cmd
    CHARACTER(LEN=MAX_NAME_LEN) :: name
    TYPE(Circuit_t), POINTER :: Circuit
    
    Circuit => CurrentModel % Circuits(CId)
    
    nofc = 0
    
    char_len = LEN_TRIM(Var_type)
    DO i=1,Circuit % n
      slen = Matc('C.'//i2s(CId)//'.name.'//i2s(i),name)
      IF(name(1:char_len) == Var_type(1:char_len)) nofc = nofc + 1
    END DO

!------------------------------------------------------------------------------
  END FUNCTION CountNofCircVarsOfType
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION CountNofCircComponents(CId, nofvar) RESULT (nofc)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: nofc, nofvar, slen, CId, i, j, CompId, ibracket
    INTEGER :: ComponentIDs(nofvar)
    TYPE(Circuit_t), POINTER :: Circuit
    CHARACTER(:), ALLOCATABLE :: cmd
    CHARACTER(LEN=MAX_NAME_LEN) :: name

    nofc = 0
    ComponentIDs = -1
    
    Circuit => CurrentModel % Circuits(CId)
    
   
    DO i=1,Circuit % n
      slen = Matc('C.'//i2s(CId)//'.name.'//i2s(i),name)

      IF(isComponentName(name,slen)) THEN
        DO ibracket=1,slen
          IF(name(ibracket:ibracket)=='(') EXIT 
        END DO

        DO j=ibracket+1,slen
          IF(name(j:j)==')') EXIT 
        END DO

        READ(name(ibracket+1:j-1),*) CompId

        IF (.NOT. ANY(ComponentIDs == CompID)) nofc = nofc + 1
        ComponentIDs(i) = CompId
      END IF
    
    END DO

!------------------------------------------------------------------------------
  END FUNCTION CountNofCircComponents
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
FUNCTION isComponentName(name, len) RESULT(L)
!------------------------------------------------------------------------------
   CHARACTER(LEN=*) :: name
   INTEGER :: len
   LOGICAL :: L
   
   L = .FALSE.
   IF(len<12) RETURN
   IF(name(1:12)=='i_component(' .OR. &
      name(1:12)=='v_component(' .OR. &
      name(1:14)=='phi_component(') L=.TRUE.
!------------------------------------------------------------------------------
END FUNCTION isComponentName
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE ReadCircuitVariables(CId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: slen, ComponentId,i,j,CId, CompInd, nofc, ibracket
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(CircuitVariable_t), POINTER :: CVar
    LOGICAL :: LondonEquations = .FALSE.
    CHARACTER(:), ALLOCATABLE :: cmd
    CHARACTER(LEN=MAX_NAME_LEN) :: name

    Circuit => CurrentModel % Circuits(CId)

    nofc = SIZE(Circuit % Components)
    DO i=1,nofc
      Circuit % Components(i) % ComponentId = 0
    END DO

    CompInd = 0
    DO i=1,Circuit % n
      slen = Matc('C.'//i2s(CId)//'.name.'//i2s(i),name)
      Circuit % names(i) = name(1:slen)

      CVar => Circuit % CircuitVariables(i)
      CVar % isIvar = .FALSE.
      CVar % isVvar = .FALSE.
      CVar % Component => Null()

      IF(isComponentName(name,slen)) THEN
        DO ibracket=1,slen
          IF(name(ibracket:ibracket)=='(') EXIT 
        END DO

        DO j=ibracket+1,slen
          IF(name(j:j)==')') EXIT 
        END DO
        READ(name(ibracket+1:j-1),*) ComponentId

        CVar % BodyId = ComponentId
        
        DO j=1,nofc
          Cvar % Component => Circuit % Components(j)
          IF(CVar % Component % ComponentId==ComponentId) EXIT
        END DO

        IF(CVar % Component % ComponentID /= ComponentId ) THEN
          CompInd = CompInd + 1
          CVar % Component => Circuit % Components(CompInd)
        END IF

        Cvar % Component % ComponentId = ComponentId

        SELECT CASE (name(1:ibracket))
        CASE('i_component(')
          CVar % isIvar = .TRUE.
          CVar % Component % ivar => CVar
        CASE('v_component(')
          LondonEquations = ListGetLogical(CurrentModel % Components (ComponentId) % Values, &
                                           'London Equations', LondonEquations)
          IF (.NOT. LondonEquations) THEN
            CVar % isVvar = .TRUE.
            CVar % Component % vvar => CVar
          ELSE
            Cvar % Component => Null()
            CVar % isIvar = .FALSE.
            CVar % isVvar = .FALSE.
            CVar % dofs = 1
            CVar % pdofs = 0
            CVar % BodyId = 0
          END IF
        CASE('phi_component(')
          ! London equations lead to driving the a-formulation 
          ! with the so called node flux. Thus we replace 'v_component'
          ! variable with phi_component:
          ! (beta a, phi') + phi_component(1) (beta grad phi_0, grad phi') = i_component(1)
          !--------------------------------------------------------------------------------
          CVar % isVvar = .TRUE.
          CVar % Component % vvar => CVar 
        CASE DEFAULT
          CALL Fatal('Circuits_Init()', 'Circuit variable should be either i_component or v_component!')
        END SELECT
      ELSE
          CVar % isIvar = .FALSE.
          CVar % isVvar = .FALSE.
          CVar % dofs = 1
          CVar % pdofs = 0
          CVar % BodyId = 0
      END IF
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE ReadCircuitVariables
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION GetNofCircVariables(CId) RESULT(n)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: CId, n, slen 
    TYPE(Circuit_t), POINTER :: Circuit

    Circuit => CurrentModel % Circuits(CId)
    Circuit % n = NINT(GetMatcReal('C.'//i2s(CId)//'.variables'))
    n = Circuit % n
!------------------------------------------------------------------------------
  END FUNCTION GetNofCircVariables
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AllocateCircuit(CId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: CId,n
    TYPE(Circuit_t), POINTER :: Circuit

    Circuit => CurrentModel % Circuits(CId)
    
    n = Circuit % n
    
    ALLOCATE(Circuit % ComponentIds(n))
    ALLOCATE(Circuit % CircuitVariables(n), Circuit % Perm(n))
    ALLOCATE(Circuit % names(n), Circuit % source(n))
    ALLOCATE(Circuit % A(n,n), Circuit % B(n,n), &
             Circuit % Mre(n,n), Circuit % Mim(n,n)  )
    Circuit % ComponentIds = 0
    Circuit % names = ' '
    Circuit % A = 0._dp
    Circuit % B = 0._dp
    Circuit % Mre = 0._dp
    Circuit % Mim = 0._dp

!------------------------------------------------------------------------------
  END SUBROUTINE AllocateCircuit
!------------------------------------------------------------------------------

!-------------------------------------------------------------------
 SUBROUTINE SetBoundaryAreasToValueLists()
!-------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Valuelist_t), POINTER :: BC
    REAL(KIND=dp), ALLOCATABLE :: BoundaryAreas(:)
    REAL(KIND=dp) :: area
    INTEGER :: Active, t0, t, i, BCid, n, nBC
    INTEGER, POINTER :: ChildBCs(:)
    LOGICAL :: Found
    LOGICAL :: Parallel 

    Mesh => CurrentModel % Mesh

    Parallel = CurrentModel % Solver % Parallel

    ! Not all boundary elements are associated to a BC.
    ! Some may also be created by extrusion. 
    t0 = Mesh % NumberOfBulkElements
    nBC = 0
    DO t=1,Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(t0+t)
      IF(ASSOCIATED(Element % BoundaryInfo) ) THEN
        nBC = MAX(nBC, Element % BoundaryInfo % Constraint )
      END IF
    END DO

    nBC = MAX(nBC,CurrentModel % NumberOfBCs)
    nBC = ParallelReduction(nBC,2)
    
    ALLOCATE( BoundaryAreas(nBC) )
    BoundaryAreas = 0.0_dp
        
    DO i=1, CurrentModel % NumberOfBcs
       BC => CurrentModel % BCs(i) % Values
       IF (.NOT. ASSOCIATED(BC) ) CALL Fatal('SetBoundaryAreasToValueLists', 'Boundary not found!')
       CALL ListAddInteger(BC, 'Boundary Id', i)
    END DO
    
    Active = GetNOFBoundaryElements()
    DO t=1,Active
       Element => GetBoundaryElement(t)
       
       BC=>GetBC()
       IF (ASSOCIATED(BC) ) THEN
         BCid = GetInteger(BC, 'Boundary Id', Found)
       ELSE
         BCid = Element % BoundaryInfo % Constraint
       END IF

       IF( BCid > 0 ) THEN
         n = GetElementNOFNodes() 
         BoundaryAreas(BCid) = BoundaryAreas(BCid) + ElementAreaNoAxisTreatment(Mesh, Element, n)
       END IF
     END DO
     
     IF( Parallel ) THEN
       DO i=1, nBC
         BoundaryAreas(i) = ParallelReduction(BoundaryAreas(i))
       END DO
     END IF

     IF( InfoActive(25) ) THEN
       DO i=1,nBC
         PRINT *,'A(i)',i,i<=CurrentModel % NumberOfBCs,BoundaryAreas(i)
       END DO
     END IF
     
     DO i=1, CurrentModel % NumberOfBcs
       BC => CurrentModel % BCs(i) % Values
       IF (.NOT. ASSOCIATED(BC) ) CALL Fatal('ComputeCoilBoundaryAreas', 'Boundary not found!')
       BCid = GetInteger(BC, 'Boundary Id', Found)
       CALL ListAddConstReal(BC, 'Area', BoundaryAreas(BCid))
     END DO
     
     DO i=1, CurrentModel % NumberOfBodies
       BC => CurrentModel % Bodies(i) % Values
       ChildBCs => ListGetIntegerArray( BC,'Extruded Child BCs',Found ) 
       IF(Found) THEN
         !PRINT *,'Child BCs area:',i,BoundaryAreas(ChildBCs)
         area = SUM(BoundaryAreas(ChildBCs)) / SIZE( ChildBCs )         
         CALL ListAddConstReal(BC,'Extruded Child Area', area )         
       END IF
     END DO
    
!-------------------------------------------------------------------
 END SUBROUTINE SetBoundaryAreasToValueLists
!-------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE ReadComponents(CId)
!------------------------------------------------------------------------------
    USE CircuitUtils
    IMPLICIT NONE
    INTEGER :: CId, CompInd
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(Component_t), POINTER :: Comp
    TYPE(Valuelist_t), POINTER :: CompParams
    LOGICAL :: Found
    INTEGER :: ExtMaster

    CALL Info('ReadComponents','Reading component: '//I2S(Cid),Level=20)

    Circuit => CurrentModel % Circuits(CId)
    
    Circuit % CvarDofs = 0
    DO CompInd=1,Circuit % n_comp
      Comp => Circuit % Components(CompInd)
      Comp % nofcnts = 0
!        Comp % ComponentId = Circuits(p) % body(CompInd)
      Comp % BodyIds => GetComponentBodyIds(Comp % ComponentId)

      IF (.NOT. ASSOCIATED(Comp % ivar) ) THEN
        CALL FATAL('Circuits_Init', 'Current Circuit Variable is not found for Component '//i2s(Comp % ComponentId))
      ELSE IF (.NOT. ASSOCIATED(Comp % vvar) ) THEN
        CALL FATAL('Circuits_Init', 'Voltage Circuit Variable is not found for Component '//i2s(Comp % ComponentId))
      END IF

      CompParams => CurrentModel % Components (Comp % ComponentId) % Values
      IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('Circuits_Init', 'Component parameters not found!')
      
      Comp % CoilType = GetString(CompParams, 'Coil Type', Found)
      IF (.NOT. Found) THEN
        CALL Info('Circuits_Init', 'Component '//i2s(Comp % ComponentId)//' is not a coil. &
          Checking if it has a component type.', Level=7)
        Comp % ComponentType = GetString(CompParams, 'Component Type', Found)
        IF (.NOT. Found) CALL Fatal ('Circuits_Init', 'Component Type not found!')
      ELSE
        Comp % ComponentType = 'coil'
      END IF
      
      Comp % i_multiplier_re = GetConstReal(CompParams, 'Current Multiplier re', Found)
      IF (.NOT. Found) Comp % i_multiplier_re = 0._dp
      Comp % i_multiplier_im = GetConstReal(CompParams, 'Current Multiplier im', Found)
      IF (.NOT. Found) Comp % i_multiplier_im = 0._dp

      Comp % VoltageFactor = GetConstReal(CompParams, 'Circuit Equation Voltage Factor', Found)
      IF (.NOT. Found) Comp % VoltageFactor = 1._dp

      Comp % ElBoundaries => ListGetIntegerArray(CompParams, 'Electrode Boundaries', Found)
      
      ! This is a feature intended to make it easier to extruded meshes internally with
      ! ElmerSolver. The idea is that the code knows which are the BCs that were created
      ! from extruding this 2D body. 
      ExtMaster = 0
      IF(.NOT. Found ) THEN
        IF( ListGetLogical( CurrentModel % Solver % Values,'Extruded Child BC Electrode', Found ) ) THEN          

          IF( ListGetLogical( CurrentModel % Simulation,"Extruded BCs Collect",Found ) ) THEN
            CALL Fatal('Circuits_init',&
                'Conflicting keywords: "Extruded Child BC Electrode" vs. "Extruded BCs Collect"')
          END IF

          CALL Info('Circuits_init','Setting "Extruded Child BCs"',Level=10)
          BLOCK
            INTEGER :: body_id
            INTEGER, POINTER :: pIntArray(:) => NULL()
            pIntArray => ListGetIntegerArray(CompParams, 'Body', Found )
            IF(.NOT. Found) pIntArray => ListGetIntegerArray(CompParams, 'Master Bodies', Found )
            IF( Found ) THEN
              IF(SIZE(pIntArray)==1) THEN
                body_id = pIntArray(1)
                NULLIFY(pIntArray)
                pIntArray => ListGetIntegerArray(CurrentModel % Bodies(body_id) % Values,&
                    'Extruded Child BCs',Found )
                IF(Found) THEN
                  CALL Info('Circuits_init','Setting Component '//I2S(CompInd)//' "Electrode Boundaries" to '&
                      //I2S(pIntArray(1))//' '//I2S(pIntArray(2)),Level=10)
                  Comp % ElBoundaries => pIntArray
                  ExtMaster = body_id
                END IF
              END IF
            END IF
          END BLOCK
        END IF
      END IF

      ! Must precede the coil type init, which already reads the
      ! direction field (ComputeFoilSheetSign).
      IF (ListGetLogical(CompParams, 'Coil Closed', Found)) THEN
        SELECT CASE (Comp % CoilType) 
        CASE ('foil winding', 'flat wire', 'foil sheet')
          CALL ComputeCoilCirculation(CompParams, CompInd)
        END SELECT
      END IF

      IF (Comp % ComponentType == 'resistor') THEN
       Comp % ivar % dofs = 1
        Comp % vvar % dofs = 1
        Comp % ivar % pdofs = 0
        Comp % vvar % pdofs = 0
      ELSE
        SELECT CASE (Comp % CoilType) 
        CASE ('stranded')
          
          Comp % nofturns = GetConstReal(CompParams, 'Number of Turns', Found)
          IF (.NOT. Found) CALL Fatal('Circuits_Init','Number of Turns not found!')
          
          Comp % ElArea = GetConstReal(CompParams, 'Electrode Area', Found)
          IF (.NOT. Found) THEN
            CALL ComputeElectrodeArea(Comp, CompParams, ExtMaster )
            WRITE(Message,'(A,ES12.5)') 'Component '//I2S(CompInd)//' "Electrode Area" is ',Comp % ElArea
            CALL Info('Circuits_Init',Message,Level=10)
          END IF
            
          Comp % CoilThickness = GetConstReal(CompParams, 'Coil Thickness', Found)
          IF (.NOT. Found) Comp % CoilThickness = 1._dp

          Comp % SymmetryCoeff = GetConstReal(CompParams, 'Symmetry Coefficient', Found)
          IF (.NOT. Found) Comp % SymmetryCoeff = 1.0_dp
          
          Comp % N_j = Comp % CoilThickness * Comp % nofturns / Comp % ElArea
          
          ! Stranded coil has current and voltage 
          ! variables (which both have a dof):
          ! ------------------------------------
          Comp % ivar % dofs = 1
          Comp % vvar % dofs = 1
          Comp % ivar % pdofs = 0
          Comp % vvar % pdofs = 0

        CASE ('massive')
          ! Massive coil has current and voltage 
          ! variables (which both have a dof):
          ! ------------------------------------
          Comp % ivar % dofs = 1
          Comp % vvar % dofs = 1
          Comp % ivar % pdofs = 0
          Comp % vvar % pdofs = 0

        CASE ('foil winding')
          Comp % polord = GetInteger(CompParams, 'Foil Winding Voltage Polynomial Order', Found)
          IF (.NOT. Found) Comp % polord = 2

          ! Foil winding has current and voltage 
          ! variables. Current has one dof and 
          ! voltage has a polynom for describing the 
          ! global voltage. The polynom has 1+"polynom order"
          ! dofs. Thus voltage variable has 1+1+"polynom order"
          ! dofs (V=V0+V1*alpha+V2*alpha^2+..):
          ! dofs:
          ! V, V0, V1, V2, ...
          ! ------------------------------------
          Comp % ivar % dofs = 1
          Comp % ivar % pdofs = 0
          Comp % vvar % dofs = Comp % polord + 2
          ! polynom dofs:
          ! -------------
          Comp % vvar % pdofs = Comp % polord + 1

          Comp % coilthickness = GetConstReal(CompParams, 'Coil Thickness', Found)
          IF (.NOT. Found) CALL Fatal('Circuits_Init','Coil Thickness not found!')

          Comp % nofturns = GetConstReal(CompParams, 'Number of Turns', Found)
          IF (.NOT. Found) CALL Fatal('Circuits_Init','Number of Turns not found!')

          Comp % ElArea = GetConstReal(CompParams, 'Electrode Area', Found)
          IF (.NOT. Found) THEN
            CALL ComputeElectrodeArea(Comp, CompParams )
            WRITE(Message,'(A,ES12.5)') 'Component '//I2S(CompInd)//' "Electrode Area" is ',Comp % ElArea
            CALL Info('Circuits_Init',Message,Level=10)
          END IF
          
          Comp % N_j = Comp % nofturns / Comp % ElArea

        CASE ('flat wire')
          CALL InitFlatWireComponent(Comp, CompParams, CompInd, ExtMaster)

        CASE ('foil sheet')
          CALL InitFoilSheetComponent(Comp, CompParams, CompInd, ExtMaster)
        END SELECT
      END IF

      CALL AddVariableToCircuit(Circuit, Comp % ivar, CId)
      CALL AddVariableToCircuit(Circuit, Comp % vvar, CId)

    END DO
      
!------------------------------------------------------------------------------
  END SUBROUTINE ReadComponents
!------------------------------------------------------------------------------

!-------------------------------------------------------------------
 SUBROUTINE ComputeElectrodeArea(Comp, CompParams, ExtMaster )
!-------------------------------------------------------------------
  USE ElementUtils
  IMPLICIT NONE
  TYPE(Component_t), POINTER :: Comp
  TYPE(ValueList_t), POINTER :: CompParams
  INTEGER, OPTIONAL :: ExtMaster 

  TYPE(ValueList_t), POINTER :: BC
  TYPE(Element_t), POINTER :: Element
  TYPE(Mesh_t), POINTER :: Mesh
  INTEGER :: t, n, BCid, NoSlices
  LOGICAL :: Found
  LOGICAL :: Parallel 
  
  Mesh => CurrentModel % Mesh
  Comp % ElArea = 0._dp

  Parallel = ( ParEnv % PEs > 1 )
  IF( Parallel ) THEN
    ! If we have single mesh then we have either parallel times or parallel slices.
    ! In both cases let us not do a parallel sum. 
    IF( Mesh % SingleMesh ) Parallel = .FALSE.
  END IF
    
  IF (CoordinateSystemDimension() == 2) THEN
    DO t=1,GetNOFActive()
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes() 
      IF (ElAssocToComp(Element, Comp)) THEN
        Comp % ElArea = Comp % ElArea + ElementAreaNoAxisTreatment(Mesh, Element, n) 
      END IF
    END DO
    
    IF( Parallel ) THEN
      Comp % ElArea = ParallelReduction(Comp % ElArea)
    END IF

    ! Add this to list since no need to compute this twice
    CALL ListAddConstReal(CompParams,'Electrode Area',Comp % ElArea )        
  ELSE
    BCid = 0

    ! This is a special case for extruded meshes where the area is computed and stored
    ! in the master body that the extruded boundary comes from. This is treated only
    ! when the extruded master is given as a parameter.
    IF( PRESENT( ExtMaster ) ) BCid = ExtMaster
    IF( BCid > 0 ) THEN
      BC => CurrentModel % Bodies(BCid) % Values
      IF (.NOT. ASSOCIATED(BC) ) CALL Fatal('ComputeElectrodeArea', 'Master body not found!')            
      Comp % ElArea = GetConstReal(BC, 'Extruded Child Area', Found ) 
      IF (.NOT. Found) CALL Fatal('ComputeElectrodeArea', '"Extruded Child Area" not found!')
    ELSE
      IF (.NOT. ASSOCIATED(Comp % ElBoundaries)) &
          CALL Fatal('ComputeElectrodeArea','Electrode Boundaries not found')      
      BCid = Comp % ElBoundaries(1)
      IF( BCid < 1 .OR. BCid > CurrentModel % NumberOfBCs ) &     
          CALL Fatal('ComputeElectrodeArea', 'BCid is beyond range: '//I2S(BCid))

      BC => CurrentModel % BCs(BCid) % Values
      IF (.NOT. ASSOCIATED(BC) ) CALL Fatal('ComputeElectrodeArea', 'Boundary not found!')
      Comp % ElArea = GetConstReal(BC, 'Area', Found)
      IF (.NOT. Found) CALL Fatal('ComputeElectrodeArea', '"Area" not found!')
    END IF
  END IF    
  
!-------------------------------------------------------------------
 END SUBROUTINE ComputeElectrodeArea
!-------------------------------------------------------------------

! This function is originally from ElementUtils. However, there is 
! some kind of treatment regarding axisymmetric cases which fails 
! here since we don't want that.
!------------------------------------------------------------------------------
   FUNCTION ElementAreaNoAxisTreatment( Mesh,Element,N ) RESULT(A)
!------------------------------------------------------------------------------
     TYPE(Mesh_t), POINTER :: Mesh
     INTEGER :: N
     TYPE(Element_t) :: Element
!------------------------------------------------------------------------------

     REAL(KIND=dp), TARGET :: NX(N),NY(N),NZ(N)

     REAL(KIND=dp) :: A

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     INTEGER :: N_Integ,t

     REAL(KIND=dp) :: Metric(3,3),Symb(3,3,3),dSymb(3,3,3,3), &
              SqrtMetric,SqrtElementMetric

     TYPE(Nodes_t) :: Nodes

     LOGICAL :: stat

     REAL(KIND=dp) :: Basis(n),u,v,w,x,y,z
     REAL(KIND=dp) :: dBasisdx(n,3)

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
!------------------------------------------------------------------------------
 
     Nodes % x => NX
     Nodes % y => NY
     Nodes % z => NZ

     Nodes % x = Mesh % Nodes % x(Element % NodeIndexes)
     Nodes % y = Mesh % Nodes % y(Element % NodeIndexes)
     Nodes % z = Mesh % Nodes % z(Element % NodeIndexes)

     IntegStuff = GaussPoints( element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ  = IntegStuff % n
!
!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
!
       A = 0.0
       DO t=1,N_Integ
!
!        Integration stuff
!
         u = U_Integ(t)
         v = V_Integ(t)
         w = W_Integ(t)
!
!------------------------------------------------------------------------------
!        Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
         stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
                    Basis,dBasisdx )
!------------------------------------------------------------------------------
!        Coordinatesystem dependent info
!------------------------------------------------------------------------------
           A =  A + SqrtElementMetric * S_Integ(t)
       END DO
!------------------------------------------------------------------------------
   END FUNCTION ElementAreaNoAxisTreatment
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
  SUBROUTINE AddVariableToCircuit(Circuit, Variable, k)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t) :: Circuit
    TYPE(CircuitVariable_t) :: Variable
    INTEGER :: Owner=-1, k
    INTEGER, POINTER :: circuit_tot_n => Null()

    CALL Info('AddVariableToCircuit','Adding variables count to circuit!',Level=20)
    
    Circuit_tot_n => CurrentModel % Circuit_tot_n
    
    IF(k==1) THEN
      IF(Owner<=0) Owner = MAX(Parenv % PEs/2,1)
      Owner = Owner - 1
      Variable % Owner = Owner
      variable % owner = 0
    ELSE
      IF(Owner<=ParEnv % PEs/2) Owner = ParEnv % PEs
      Owner = Owner - 1
      Variable % Owner = Owner
      variable % owner = ParEnv % PEs-1
    END IF

    IF (Circuit % Harmonic) THEN
      IF (Circuit % UsePerm) THEN
        Variable % valueId = Circuit % Perm(Circuit_tot_n + 1)
        Variable % ImValueId = Variable % valueId + 1
      ELSE
        Variable % valueId = Circuit_tot_n + 1
        Variable % ImValueId = Circuit_tot_n + 2
      END IF
    
      Circuit_tot_n = Circuit_tot_n + 2*Variable % dofs
    ELSE
      IF (Circuit % UsePerm) THEN
        Variable % valueId = Circuit % Perm(Circuit_tot_n + 1)
      ELSE
        Variable % valueId = Circuit_tot_n + 1
      END IF
      
      Circuit_tot_n = Circuit_tot_n + Variable % dofs
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE AddVariableToCircuit
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE AddComponentValuesToLists(CId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(Component_t), POINTER :: Comp
    TYPE(Valuelist_t), POINTER :: CompParams
    INTEGER :: CId, CompInd
    
    Circuit => CurrentModel % Circuits(CId)

    CALL Info('AddComponentValuesToLists','Adding "Circuit Voltage Variable *" keywords for '&
        //I2S(Circuit % n_comp)//' components in Circuit '//I2S(CId),Level=20)
    
    DO CompInd=1,Circuit % n_comp
 
      Comp => Circuit % Components(CompInd)   

      CompParams => CurrentModel % Components (Comp % ComponentId) % Values
      IF (.NOT. ASSOCIATED(CompParams)) CALL Fatal ('Circuits_Init', 'Component Parameters not found!')

      CALL listAddInteger(CompParams, 'Circuit Voltage Variable Id', Comp % vvar % valueId)
      CALL listAddInteger(CompParams, 'Circuit Voltage Variable dofs', Comp % vvar % dofs)
      CALL listAddInteger(CompParams, 'Circuit Current Variable Id', Comp % ivar % valueId)
      CALL listAddInteger(CompParams, 'Circuit Current Variable dofs', Comp % ivar % dofs)
      CALL listAddConstReal(CompParams, 'Stranded Coil N_j', Comp % N_j)
      CurrentModel % Components (Comp % ComponentId) % Values => CompParams
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE AddComponentValuesToLists
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE AddBareCircuitVariables(CId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(CircuitVariable_t), POINTER :: CVar
    INTEGER :: CId, i
    
    Circuit => CurrentModel % Circuits(CId)
    ! add variables that are not associated to components
    DO i=1,Circuit % n
      Cvar => Circuit % CircuitVariables(i)
      IF (Cvar % isIvar .OR. Cvar % isVvar) CYCLE
      CALL AddVariableToCircuit(Circuit, Circuit % CircuitVariables(i), CId)
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE AddBareCircuitVariables
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReadCoefficientMatrices(CId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: CId,n
    TYPE(Circuit_t), POINTER :: Circuit

    Circuit => CurrentModel % Circuits(CId)
    n = Circuit % n

    ! Read in the coefficient matrices for the circuit equations:
    ! Ax' + Bx = source:
    ! ------------------------------------------------------------

    CALL matc_get_array('C.'//i2s(CId)//'.A'//CHAR(0),Circuit % A,n,n)
    CALL matc_get_array('C.'//i2s(CId)//'.B'//CHAR(0),Circuit % B,n,n)

    IF (Circuit % Harmonic) THEN
      ! Complex multiplier matrix is used for:
      ! B = times(M,B), where B times is the element-wise product
      ! ---------------------------------------------------------
      CALL matc_get_array('C.'//i2s(CId)//'.Mre'//CHAR(0),Circuit % Mre,n,n)
      CALL matc_get_array('C.'//i2s(CId)//'.Mim'//CHAR(0),Circuit % Mim,n,n)
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE ReadCoefficientMatrices
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReadPermutationVector(CId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: CId,n,slen,i
    TYPE(Circuit_t), POINTER :: Circuit

    Circuit => CurrentModel % Circuits(CId)
    n = Circuit % n

    DO i=1,n
      Circuit % Perm(i) = NINT(GetMatcReal('C.'//i2s(CId)//'.perm('//i2s(i-1)//')'))
    END DO
    IF(ANY(Circuit % Perm /= 0)) THEN 
      Circuit % UsePerm = .TRUE.
      CALL Info( 'IHarmonic2D','Found Permutation vector for circuit '//i2s(CId), Level=4 )
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE ReadPermutationVector
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReadCircuitSources(CId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: CId,n,slen,i
    TYPE(Circuit_t), POINTER :: Circuit
    CHARACTER(:), ALLOCATABLE :: cmd
    CHARACTER(LEN=MAX_NAME_LEN) :: name

    Circuit => CurrentModel % Circuits(CId)
    n = Circuit % n
    DO i=1,n
      ! Names of the source functions, these functions should be found
      ! in the "Body Force 1" block of the .sif file.
      ! (nc: is for 'no check' e.g. don't abort if the MATC variable is not found!)
      ! ---------------------------------------------------------------------------
      slen = Matc('nc:C.'//i2s(CId)//'.source.'//i2s(i),name)
      Circuit % Source(i) = name(1:slen)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ReadCircuitSources
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE WriteCoeffVectorsForCircVariables(CId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: CId,n,slen,i,j,RowId
    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(CircuitVariable_t), POINTER :: Cvar
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)
  
    Circuit => CurrentModel % Circuits(CId)
    n = Circuit % n

    DO i=1,n
      Cvar => Circuit % CircuitVariables(i)
      RowId = Cvar % ValueId

      ALLOCATE(Cvar % A(Circuit % n), &
               Cvar % B(Circuit % n), &
               Cvar % Mre(Circuit % n), &
               Cvar % Mim(Circuit % n), &
               Cvar % SourceRe(Circuit % n), &
               Cvar % SourceIm(Circuit % n))
      Cvar % A = 0._dp
      Cvar % B = 0._dp
      Cvar % Mre = 0._dp
      Cvar % Mim = 0._dp
      Cvar % SourceRe = 0._dp
      Cvar % SourceIm = 0._dp

      DO j=1,Circuit % n
        IF (Circuit % A(i,j)/=0) Cvar % A(j) = Circuit % A(i,j)
        IF (Circuit % B(i,j)/=0) Cvar % B(j) = Circuit % B(i,j)
        IF (Circuit % Mre(i,j)/=0 .OR. Circuit % Mim(i,j)/=0) THEN
          Cvar % Mre(j) = Circuit % Mre(i,j) 
          Cvar % Mim(j) = Circuit % Mim(i,j)
        END IF
      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE WriteCoeffVectorsForCircVariables
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   FUNCTION IdInList(Id, List) RESULT (T)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     INTEGER :: List(:), Id
     LOGICAL :: T
     T = .FALSE.
     IF (ANY(List == Id)) T = .TRUE.
!------------------------------------------------------------------------------
   END FUNCTION IdInList
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   FUNCTION ElAssocToComp(Element, Component) RESULT (T)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     TYPE(Component_t), POINTER :: Component
     TYPE(Element_t), POINTER :: Element
     INTEGER :: k
     LOGICAL :: T, Found

     k = GetInteger(GetBC(Element), 'Component', Found)
     
     IF (Found) THEN
       T = (k .eq. Component % ComponentId)
     ELSE IF (ASSOCIATED(Component % BodyIds)) THEN
       T = IdInList(Element % BodyId, Component % BodyIds)
     ELSE
       T = .False.
     END IF

!------------------------------------------------------------------------------
   END FUNCTION ElAssocToComp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   FUNCTION ElAssocToCvar(Element, Cvar) RESULT (T)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     TYPE(CircuitVariable_t), POINTER :: Cvar
     TYPE(Element_t), POINTER :: Element
     LOGICAL :: T
     T = .FALSE.
     IF (ASSOCIATED(Cvar % Component)) THEN
       IF (ASSOCIATED(Cvar % Component % BodyIds)) &
       T = IdInList(Element % BodyId, Cvar % Component % BodyIds)
     END IF
!------------------------------------------------------------------------------
   END FUNCTION ElAssocToCvar
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION AddIndex(Ind, Harmonic)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    Integer :: Ind, AddIndex
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm
    
    IF (.NOT. PRESENT(Harmonic)) THEN
      harm = CurrentModel % HarmonicCircuits
    ELSE
      harm = Harmonic
    END IF
 
    IF (harm) THEN
      AddIndex = 2 * Ind
    ELSE
      AddIndex = Ind
    END IF
!------------------------------------------------------------------------------
  END FUNCTION AddIndex 
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION AddImIndex(Ind)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: Ind
    Integer :: AddImIndex
    IF ( .NOT. CurrentModel % HarmonicCircuits ) CALL Fatal ('AddImIndex','Model is not of harmonic type!')
    
    AddImIndex = 2 * Ind + 1
!------------------------------------------------------------------------------
  END FUNCTION AddImIndex 
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION ReIndex(Ind, Harmonic)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: Ind, ReIndex
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm
    
    IF (.NOT. PRESENT(Harmonic)) THEN
      harm = CurrentModel % HarmonicCircuits
    ELSE
      harm = Harmonic
    END IF
 
    IF (harm) THEN
      ReIndex = 2 * Ind - 1
    ELSE
      ReIndex = Ind
    END IF
!------------------------------------------------------------------------------
  END FUNCTION ReIndex 
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  FUNCTION ImIndex(Ind)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    Integer :: Ind, ImIndex

    ImIndex = 2 * Ind
!------------------------------------------------------------------------------
  END FUNCTION ImIndex
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   FUNCTION HasSupport(Element, nn) RESULT(support)
!------------------------------------------------------------------------------
    USE CircuitUtils
    IMPLICIT NONE
    INTEGER :: nn, dim
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: support, First=.TRUE., Found
    REAL(KIND=dp) :: wBase(nn)
    TYPE(ValueList_t), POINTER :: CompParams
    SAVE dim, First

    IF (First) THEN
      First = .FALSE.
      dim = CoordinateSystemDimension()
    END IF
    
    support = .TRUE. 
    IF (dim == 3) THEN
      ! A closed coil has no electrode potential W at all; its
      ! direction field is the CoilSolver cut potential, which is defined on
      ! every element of the coil.
      CompParams => GetComponentParams(Element)
      IF (ASSOCIATED(CompParams)) THEN
        IF (ListGetLogical(CompParams,'Coil Closed',Found)) RETURN
      END IF

      support = .FALSE.
      CALL GetLocalSolution(Wbase,'W')
      IF ( ANY(Wbase .ne. 0d0) ) support = .TRUE.
    END IF
!------------------------------------------------------------------------------
   END FUNCTION HasSupport
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Create a standard variable associated to the mesh that may be use for dependencies.
!------------------------------------------------------------------------------
  SUBROUTINE Circuits_ToMeshVariable(Solver,crt)
    
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: Crt(:)

    TYPE(Circuit_t), POINTER :: Circuit
    TYPE(CircuitVariable_t), POINTER :: CVar
    TYPE(Variable_t), POINTER :: Var, VarIm
    INTEGER :: p,i,n,nv,ni,m,iv,nsize
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: Found 
    CHARACTER(:), ALLOCATABLE :: CrtName,VarName,VarnameIm
    
    IF( .NOT. ListGetLogical( Solver % Values,'Export Circuit Variables',Found ) ) RETURN

    CALL Info('Circuit_ToMeshVariable','Adding circuit variables to be mesh variables')
    
    Mesh => Solver % Mesh
            
    DO p=1,CurrentModel % n_Circuits
      CALL Info('Circuit_ToMeshVariable','Adding circuit: '//I2S(p),Level=12)

      Circuit => CurrentModel % Circuits(p)

      n = Circuit % n
      
      IF( CurrentModel % n_Circuits == 1) THEN
        crtName = 'crt'
      ELSE
        crtName = 'crt '//I2S(p)
      END IF
    
      ! Count the v and i variables of the circuit.
      nv = 0; ni = 0
      DO i=1,n
        Cvar => Circuit % CircuitVariables(i)
        IF(Cvar % isIvar ) nv = nv + 1
        IF(Cvar % isVvar) ni = ni + 1
      END DO
      
      IF( nv + ni == 0 ) THEN
        CALL Warn('Circuits_ToMeshVariable','No voltage or current variables exists!')
        CYCLE
      END IF
      

      ! Go first through currents and then through voltages    
      DO iv=1,2      
        IF( Circuit % Harmonic ) THEN
          IF(iv==1) THEN
            varname =  crtname//' i re'
            varnameim =  crtname//' i im'
          ELSE
            varname = crtname//' v re'
            varnameim = crtname//' v im'
          END IF
        ELSE
          IF(iv==1) THEN
            varname =  crtname//' i'
          ELSE
            varname = crtname//' v'
          END IF
        END IF
        
        IF(iv==1) THEN
          nsize = ni
        ELSE
          nsize = nv
        END IF
                
        ! Get variable, if variable does not exist then we create here on-the-fly
        Var => VariableGet( Mesh % Variables,varname)
        IF(.NOT. ASSOCIATED( Var ) ) THEN
          CALL Info('Circuits_ToMeshVariable','Creating variable: '//TRIM(varname))
          CALL VariableAddVector( Mesh % Variables, Mesh, Solver,&
              varname,dofs=nsize,global=.TRUE.)
          Var => VariableGet( Mesh % Variables,varname)
        END IF
        CALL Info('Circuts_toMeshVariable','Filling variable: '//TRIM(VarName),Level=20)

        IF( Circuit % Harmonic ) THEN
          VarIm => VariableGet( Mesh % Variables,varnameim)
          IF(.NOT. ASSOCIATED( VarIm ) ) THEN
            CALL Info('Circuits_ToMeshVariable','Creating variable: '//TRIM(varnameim))
            CALL VariableAddVector( Mesh % Variables, Mesh, Solver,&
                varnameim,dofs=nsize,global=.TRUE.)
            VarIm => VariableGet( Mesh % Variables,varnameim)
          END IF
          CALL Info('Circuts_toMeshVariable','Filling variable: '//TRIM(VarNameim),Level=20)
       END IF
          
        
        ! Fill the currents or voltages
        m = 0
        DO i=1,n
          Cvar => Circuit % CircuitVariables(i)
          
          IF(iv==1 .AND. .NOT. CVar % isIvar ) CYCLE          
          IF(iv==2 .AND. .NOT. Cvar % isVvar) CYCLE
          
          CALL Info('Circuts_toMeshVariable','Inserting variable '//I2S(CVar % ValueId)//': '&
              //TRIM(Circuit % names(i)),Level=20)
                    
          m = m + 1
          Var % Values(m) = crt(Cvar % ValueId)          
          IF(Circuit % Harmonic) THEN
            VarIm % Values(m) = crt(Cvar % ImValueId)
          END IF
        END DO
      END DO
    END DO
          
  END SUBROUTINE Circuits_ToMeshVariable

   
!------------------------------------------------------------------------------
!> Flat wire winding: the coil block is split into turn cells by the normalized
!> direction fields Alpha and Beta. The stacking direction (default Beta) is
!> divided into nStack cells and the other in-plane direction into nAcross
!> cells. Every cell carries the full component current and has its own
!> voltage dof, so the component voltage is V = sum_k V_k.
!>
!> Component keywords:
!>   Number of Turns        total turns, must be an integer (= number of cells)
!>   Stacking Direction     beta (default) or alpha: direction field along the stack
!>   Turns Across           cells across the stack (default 1), divides Number of Turns
!>   Electrode Boundaries   as for foil, used for the electrode area / DC resistance
!>   Fill Factor            conductor volume fraction along the wire (default 1), or
!>   Conductor Thickness    [+ Conductor Width] to compute it from the block lengths
!>
!> Circuit dofs: vvar = V, then V_1..V_n (AddIndex(k) offsets; harmonic uses 2k).
!> Cell k of a Gauss point: FlatWireCellIndex(nStack, nAcross, sStack, sAcross).
!------------------------------------------------------------------------------
  SUBROUTINE InitFlatWireComponent(Comp, CompParams, CompInd, ExtMaster)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Component_t), POINTER :: Comp
    TYPE(ValueList_t), POINTER :: CompParams
    INTEGER :: CompInd, ExtMaster
    INTEGER :: nturns
    LOGICAL :: Found
    CHARACTER(LEN=MAX_NAME_LEN) :: str

    IF (CoordinateSystemDimension() /= 3) &
        CALL Fatal('Circuits_Init','Flat wire coil type is implemented only in 3D!')

    Comp % nofturns = GetConstReal(CompParams, 'Number of Turns', Found)
    IF (.NOT. Found) CALL Fatal('Circuits_Init','Flat wire: Number of Turns not found!')
    nturns = NINT(Comp % nofturns)
    IF (nturns < 1 .OR. ABS(Comp % nofturns - nturns) > 1.0d-8) &
        CALL Fatal('Circuits_Init','Flat wire: Number of Turns must be a positive integer!')

    Comp % nAcross = GetInteger(CompParams, 'Turns Across', Found)
    IF (.NOT. Found) Comp % nAcross = 1
    IF (Comp % nAcross < 1 .OR. MOD(nturns, Comp % nAcross) /= 0) &
        CALL Fatal('Circuits_Init','Flat wire: Number of Turns must be divisible by Turns Across!')
    Comp % nStack = nturns / Comp % nAcross

    str = GetString(CompParams, 'Stacking Direction', Found)
    IF (.NOT. Found) str = 'beta'
    SELECT CASE (str)
    CASE ('alpha')
      Comp % StackAlongAlpha = .TRUE.
    CASE ('beta')
      Comp % StackAlongAlpha = .FALSE.
    CASE DEFAULT
      CALL Fatal('Circuits_Init','Flat wire: Stacking Direction must be alpha or beta!')
    END SELECT

    IF (.NOT. ASSOCIATED(VariableGet(CurrentModel % Mesh % Variables, 'Alpha'))) &
        CALL Fatal('Circuits_Init','Flat wire needs the direction field "Alpha"!')
    IF (.NOT. ASSOCIATED(VariableGet(CurrentModel % Mesh % Variables, 'Beta'))) &
        CALL Fatal('Circuits_Init','Flat wire needs the direction field "Beta"!')

    ! Current has one dof, voltage has the total and one dof per turn cell:
    ! dofs: V, V_1, V_2, ..., V_nturns
    Comp % ivar % dofs = 1
    Comp % ivar % pdofs = 0
    Comp % vvar % dofs = nturns + 1
    Comp % vvar % pdofs = nturns

    Comp % coilthickness = 1._dp

    Comp % ElArea = GetConstReal(CompParams, 'Electrode Area', Found)
    IF (.NOT. Found) THEN
      CALL ComputeElectrodeArea(Comp, CompParams, ExtMaster)
      WRITE(Message,'(A,ES12.5)') 'Component '//I2S(CompInd)//' "Electrode Area" is ',Comp % ElArea
      CALL Info('Circuits_Init',Message,Level=10)
    END IF
    Comp % N_j = Comp % nofturns / Comp % ElArea

    Comp % FillFactor = GetConstReal(CompParams, 'Fill Factor', Found)
    IF (.NOT. Found) THEN
      IF (ListCheckPresent(CompParams, 'Conductor Thickness')) THEN
        CALL ComputeFlatWireFillFactor(Comp, CompParams)
      ELSE
        Comp % FillFactor = 1._dp
      END IF
    END IF
    WRITE(Message,'(A,ES12.5)') 'Component '//I2S(CompInd)//' flat wire fill factor is ',Comp % FillFactor
    CALL Info('Circuits_Init',Message,Level=6)

    CALL ListAddConstReal(CompParams, 'Flat Wire Fill Factor', Comp % FillFactor)
    CALL ListAddInteger(CompParams, 'Flat Wire Stack Cells', Comp % nStack)
    CALL ListAddInteger(CompParams, 'Flat Wire Across Cells', Comp % nAcross)
    CALL ListAddLogical(CompParams, 'Flat Wire Stack Along Alpha', Comp % StackAlongAlpha)
!------------------------------------------------------------------------------
  END SUBROUTINE InitFlatWireComponent
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Foil sheet winding: the turn stack is modeled as conducting SHEETS, not as a
!> conducting continuum. The block is split into nCells cells along the stacking
!> direction (the turn normal) and each cell into nSegments strands across it
!> (the turn width). Cell k lumps foilsPerCell = N/nCells turns; strand (k,j)
!> carries a uniform current density c_kj*grad(W). The intra-turn skin effect is
!> in the complex sheet conductivity (Sigma 33 / Sigma 33 im) and the intra-turn
!> proximity effect in the complex reluctivity (Nu 11/22/33), so the block
!> itself carries no volumetric eddy current.
!>
!> 'Stacking Direction' picks which direction field is the stacking normal:
!> alpha (the default) is a foil stack or a flatwise flat wire, beta is an
!> edgewise flat wire. Everything that distinguishes the stacking from the
!> across direction follows it - the cell and strand layout, the strand
!> direction, and which reluctivity components carry the complex stack
!> permeability.
!>
!> Component keywords:
!>   Number of Turns        number of foils N, must be an integer
!>   Stacking Direction     alpha (default, foils and flatwise flat wire) or
!>                          beta (edgewise flat wire): the field across the stack
!>   Sheet Cells            cells along the stacking direction (default N), must divide N
!>   Sheet Segments         strands across the stack per cell (default 16)
!>   Sheet Sublayers        strand layers through the thickness of one turn
!>                          (default 1 = one uniform current layer per turn,
!>                          0 = chosen from t/delta at the run frequency)
!>   Electrode Area         or Electrode Boundaries, as for foil winding
!>   Sigma 33 [im]          complex sheet conductivity (harmonic)
!>   Homogenization Model + Nu 11/22/33 [im]   complex reluctivity via RotM
!>
!> Circuit dofs: vvar = V, V_1..V_nCells, then c_kj cell by cell.
!------------------------------------------------------------------------------
!> True when every independent variable a keyword depends on can be resolved.
!> ListGetReal Fatals on a missing one, and 'Sheet Conductivity' is written as a
!> function of Temperature by a SIF that may also be run without any temperature
!> field, so the dependency must be checked before the first evaluation. The
!> accepted forms mirror ListParseStrToVars: a variable, 'coordinate',
!> 'prev <var>', a plain number, or another keyword of the same list.
!------------------------------------------------------------------------------
  FUNCTION FoilSheetDepsResolved(ptr, CompParams, Missing) RESULT(Ok)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueListEntry_t), POINTER :: ptr
    TYPE(ValueList_t), POINTER :: CompParams
    CHARACTER(LEN=MAX_NAME_LEN) :: Missing
    LOGICAL :: Ok
    TYPE(ValueListEntry_t), POINTER :: kptr
    TYPE(Variable_t), POINTER :: Var
    LOGICAL :: Found
    INTEGER :: l0, l1, slen

    Ok = .TRUE.
    Missing = ' '
    slen = ptr % DepNameLen
    IF (slen <= 0) RETURN

    l0 = 1
    DO WHILE (.TRUE.)
      DO WHILE (ptr % DependName(l0:l0) == ' ')
        l0 = l0 + 1
        IF (l0 > slen) EXIT
      END DO
      IF (l0 > slen) EXIT

      l1 = INDEX(ptr % DependName(l0:slen), ',')
      IF (l1 > 0) THEN
        l1 = l0 + l1 - 2
      ELSE
        l1 = slen
      END IF

      IF (ptr % DependName(l0:l1) /= 'coordinate') THEN
        Var => VariableGet(CurrentModel % Variables, TRIM(ptr % DependName(l0:l1)))
        IF (.NOT. ASSOCIATED(Var) .AND. l1-l0 > 5) THEN
          IF (ptr % DependName(l0:l0+4) == 'prev ') &
              Var => VariableGet(CurrentModel % Variables, TRIM(ptr % DependName(l0+5:l1)))
        END IF
        IF (.NOT. ASSOCIATED(Var)) THEN
          IF (VERIFY(ptr % DependName(l0:l1),'-.0123456789eE') /= 0) THEN
            kptr => ListFind(CompParams, ptr % DependName(l0:l1), Found)
            IF (.NOT. Found) THEN
              Ok = .FALSE.
              Missing = ptr % DependName(l0:l1)
              RETURN
            END IF
          END IF
        END IF
      END IF

      l0 = l1 + 2
      IF (l0 > slen) EXIT
    END DO
!------------------------------------------------------------------------------
  END FUNCTION FoilSheetDepsResolved
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> 'Sheet Conductivity' as one number for the whole component. A constant takes
!> the plain path and stays bit identical. A Real function -- TRAFOLO writes the
!> winding conductivity as sigma(Temperature) and re-solves the electromagnetic
!> problem inside a thermal iteration -- is evaluated at the nodes of the block
!> elements and averaged over the block volume, MPI reduced so that every
!> partition derives exactly the same material. 'Varies' reports the function
!> case, which is what makes the value worth recomputing on every solver call.
!------------------------------------------------------------------------------
  FUNCTION FoilSheetBlockConductivity(CompParams, Found, Varies) RESULT(sgm)
!------------------------------------------------------------------------------
    USE CircuitUtils
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: CompParams
    LOGICAL :: Found, Varies
    REAL(KIND=dp) :: sgm
    TYPE(ValueListEntry_t), POINTER :: ptr
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), sNod(:)
    REAL(KIND=dp), POINTER :: sPtr(:)
    REAL(KIND=dp) :: sInt, vol, detJ, smin, smax
    CHARACTER(LEN=MAX_NAME_LEN) :: Missing
    INTEGER :: e, n, gp, nmax
    LOGICAL :: stat, Found2

    sgm = 0._dp
    Varies = .FALSE.
    ptr => ListFind(CompParams, 'Sheet Conductivity', Found)
    IF (.NOT. Found) RETURN

    IF (ptr % TYPE == LIST_TYPE_CONSTANT_SCALAR .OR. &
        ptr % TYPE == LIST_TYPE_CONSTANT_SCALAR_STR) THEN
      sgm = ListGetConstReal(CompParams, 'Sheet Conductivity', Found)
      RETURN
    END IF

    ! The keyword is a function. Without its argument there is nothing to
    ! evaluate, so say so and let the caller keep the value it already has
    ! rather than Fatal deep inside the list machinery.
    IF (.NOT. FoilSheetDepsResolved(ptr, CompParams, Missing)) THEN
      CALL Warn('Circuits_Init','Foil sheet "Sheet Conductivity" depends on ['// &
          TRIM(Missing)//'], which does not exist; keeping the previous value')
      Found = .FALSE.
      RETURN
    END IF

    Varies = .TRUE.
    Mesh => CurrentModel % Mesh
    nmax = Mesh % MaxElementNodes
    ALLOCATE(Basis(nmax), dBasisdx(nmax,3), sNod(nmax))
    sInt = 0._dp
    vol = 0._dp
    smin = HUGE(1._dp)
    smax = -HUGE(1._dp)

    DO e = 1, GetNOFActive()
      Element => GetActiveElement(e)
      IF (.NOT. ASSOCIATED(GetComponentParams(Element), CompParams)) CYCLE
      n = GetElementNOFNodes(Element)
      sPtr => GetReal(CompParams, 'Sheet Conductivity', Found2, Element)
      IF (.NOT. Found2) CYCLE
      sNod(1:n) = sPtr(1:n)
      smin = MIN(smin, MINVAL(sNod(1:n)))
      smax = MAX(smax, MAXVAL(sNod(1:n)))
      CALL GetElementNodes(Nodes, Element)
      IP = GaussPoints(Element)
      DO gp = 1, IP % n
        stat = ElementInfo(Element, Nodes, IP % U(gp), IP % V(gp), IP % W(gp), &
            detJ, Basis, dBasisdx)
        sInt = sInt + IP % s(gp) * detJ * SUM(sNod(1:n)*Basis(1:n))
        vol  = vol  + IP % s(gp) * detJ
      END DO
    END DO
    DEALLOCATE(Basis, dBasisdx, sNod)

    sInt = ParallelReduction(sInt)
    vol  = ParallelReduction(vol)
    smin = ParallelReduction(smin, 1)
    smax = ParallelReduction(smax, 2)
    IF (vol <= 0._dp) THEN
      CALL Warn('Circuits_Init','Foil sheet block has no volume, cannot average "Sheet Conductivity"')
      Found = .FALSE.
      RETURN
    END IF

    ! A function that happens to be uniform over the block -- a sigma(T) written
    ! with a zero temperature coefficient, or a uniform temperature -- must give
    ! back exactly the number it evaluates to, so that it is bit identical to
    ! the same value written as a constant. sInt/vol is only equal to it to
    ! rounding.
    IF (smax == smin) THEN
      sgm = smin
    ELSE
      sgm = sInt / vol
    END IF
!------------------------------------------------------------------------------
  END FUNCTION FoilSheetBlockConductivity
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> How many strand layers one turn needs at this frequency, and whether the mesh
!> can carry them. One uniform current layer per turn is exact while the turn is
!> thin against the skin depth; past that the current redistributes through the
!> thickness and the single layer under-predicts Rac, so the layer count follows
!> the thickness in skin depths.
!>
!> Sub-layers only pay off when the FE field resolves the flux between the
!> bands, which is what drives the current from one band to the next. The
!> criterion is therefore the element size against the skin depth, not the band
!> thickness against the element: on a coarser mesh the layers are vetoed and
!> the single-layer model, with its own known error, is the better of the two.
!------------------------------------------------------------------------------
  FUNCTION FoilSheetAutoSublayers(tfoil, sgm, elemH) RESULT(m)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: tfoil, sgm, elemH
    INTEGER :: m
    REAL(KIND=dp) :: omega, delta, ratio, mu0
    REAL(KIND=dp), POINTER :: fptr(:,:)
    LOGICAL :: Found

    mu0 = 4.0d-7 * PI

    ! Transient has no single frequency to size the layers with, and the strand
    ! skin ladder and the reluctivity ladder are derived for one layer per turn.
    IF (TRIM(ListGetString(CurrentModel % Simulation,'Simulation Type',Found)) == 'transient') THEN
      m = 1
      CALL Info('Circuits_Init', &
          'Foil sheet transient: automatic sub-layers resolve to one layer per turn', Level=3)
      RETURN
    END IF

    ! The count fixes the circuit dof layout, so it is resolved once at init.
    ! A SIF that scans several frequencies would keep the first one's choice.
    fptr => ListGetConstRealArray(CurrentModel % Simulation, 'Frequency', Found)
    IF (Found) THEN
      IF (SIZE(fptr) > 1) CALL Warn('Circuits_Init', &
          'Foil sheet: "Sheet Sublayers = 0" is resolved once at init, but this '// &
          'simulation carries several frequencies; give "Sheet Sublayers" explicitly.')
    END IF
    omega = GetAngularFrequency(Found = Found)
    IF (.NOT. Found) omega = 0._dp
    ratio = 0._dp
    delta = HUGE(delta)
    IF (omega > 0._dp .AND. sgm > 0._dp) THEN
      delta = SQRT(2._dp / (omega * mu0 * sgm))
      ratio = tfoil / delta
    END IF

    IF (ratio <= 0.5_dp) THEN
      m = 1
    ELSE IF (ratio <= 2._dp) THEN
      m = 3
    ELSE
      m = 5
    END IF
    WRITE(Message,'(A,F8.3,A,I0)') 'Foil sheet automatic sub-layers: t/delta = ', ratio, &
        ' gives Sheet Sublayers = ', m
    CALL Info('Circuits_Init', Message, Level=3)

    IF (m > 1 .AND. elemH > delta) THEN
      WRITE(Message,'(A,F8.3,A)') 'Foil sheet: the coil mesh is too coarse for sub-layers, ' // &
          'element / skin depth = ', elemH / delta, &
          '; using one layer per turn and its own skin factor'
      CALL Info('Circuits_Init', Message, Level=3)
      m = 1
    END IF
!------------------------------------------------------------------------------
  END FUNCTION FoilSheetAutoSublayers
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Harmonic sheet material from the physics: the complex sheet conductivity and
!> the complex in-plane reluctivity at the current angular frequency, written
!> only where the SIF did not give them (the ownership is decided once, in
!> InitFoilSheetMaterial, before anything is written). Shared by the init and by
!> the per-call refresh of a temperature dependent 'Sheet Conductivity'.
!------------------------------------------------------------------------------
  SUBROUTINE SetFoilSheetHarmonicMaterial(CompParams, tfoil, ff, sgm, StackAlongAlpha, nSub)
!------------------------------------------------------------------------------
    USE CircuitUtils
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: CompParams
    REAL(KIND=dp) :: tfoil, ff, sgm
    LOGICAL :: StackAlongAlpha
    INTEGER :: nSub
    REAL(KIND=dp) :: tau0, sdc, omega, mu0, tsub
    COMPLEX(KIND=dp) :: u, th, sigs, mue, nue
    LOGICAL :: Homog, Found, FoundFreq
    INTEGER :: dStack, dPlane(2), d
    COMPLEX(KIND=dp), PARAMETER :: im = (0._dp,1._dp)

    ! Every strand is one sub-layer of a turn, so the 1-D skin solution the
    ! sheet coefficients come from is the one of a plate of thickness t/nSub.
    mu0 = 4.0d-7 * PI
    tsub = tfoil / nSub
    tau0 = mu0 * sgm * tsub**2 / 4._dp
    ! One band per turn smears the copper over the pitch, so its conductivity
    ! carries the fill factor. Copper only bands are the copper itself.
    IF (nSub > 1) THEN
      sdc = sgm
    ELSE
      sdc = ff * sgm
    END IF
    CALL ListAddConstReal(CompParams, 'Foil Sheet Tau0', tau0)
    CALL ListAddConstReal(CompParams, 'Foil Sheet Sigma DC', sdc)

    omega = GetAngularFrequency(Found = FoundFreq)
    IF (.NOT. FoundFreq) omega = 0._dp
    IF (omega < 0._dp) RETURN

    ! Frequency = 0 is an ordinary case: the harmonic solver is driven there for
    ! the DC resistance and inductance points. FoilSheetTanhOverU is 1 at omega
    ! = 0, which leaves the real DC sheet conductivity ff*sigma and the plain air
    ! reluctivity, the exact DC limit.
    u  = SQRT(im * omega * tau0)
    th = FoilSheetTanhOverU(u)

    IF (GetLogical(CompParams, 'Foil Sheet Derive Sigma 33', Found)) THEN
      sigs = sdc * th
      CALL ListAddConstReal(CompParams, 'Sigma 33', REAL(sigs, KIND=dp))
      CALL ListAddConstReal(CompParams, 'Sigma 33 im', AIMAG(sigs))
      WRITE(Message,'(A,ES12.5,SP,ES12.5,A)') 'Foil sheet Sigma 33 = ', &
          REAL(sigs, KIND=dp), AIMAG(sigs), ' i'
      CALL Info('Circuits_Init', Message, Level=3)
    END IF

    Homog = GetLogical(CompParams, 'Homogenization Model', Found)
    IF (.NOT. Found) Homog = .FALSE.
    IF (Homog .AND. GetLogical(CompParams, 'Foil Sheet Derive Nu', Found)) THEN
      mue = mu0 * ((1._dp - ff) + ff * th)
      nue = 1._dp / mue
      CALL FoilSheetNuDirections(StackAlongAlpha, dStack, dPlane)
      ! A thin layer drives no eddy loop for a field along its own normal, so
      ! only the two in-plane directions get the complex stack permeability.
      CALL ListAddConstReal(CompParams, FoilSheetNuKey(dStack), 1._dp / mu0)
      DO d = 1, 2
        CALL ListAddConstReal(CompParams, FoilSheetNuKey(dPlane(d)), REAL(nue, KIND=dp))
        CALL ListAddConstReal(CompParams, FoilSheetNuKey(dPlane(d))//' im', AIMAG(nue))
      END DO
      WRITE(Message,'(A,ES12.5,SP,ES12.5,A)') 'Foil sheet '//FoilSheetNuKey(dPlane(1))// &
          ' = '//FoilSheetNuKey(dPlane(2))//' = ', REAL(nue, KIND=dp), AIMAG(nue), ' i'
      CALL Info('Circuits_Init', Message, Level=3)
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE SetFoilSheetHarmonicMaterial
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Refresh the harmonic sheet material from a temperature dependent
!> 'Sheet Conductivity'. Called on every harmonic circuits solver call, so a
!> thermal iteration that updates the temperature field also updates the coil
!> material. A constant keyword, an explicit 'Sigma 33' or a transient run all
!> leave this a no-op.
!------------------------------------------------------------------------------
  SUBROUTINE UpdateFoilSheetMaterial(CompParams)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: CompParams
    REAL(KIND=dp) :: tfoil, ff, sgm
    INTEGER :: nSub
    LOGICAL :: FoundT, FoundF, FoundS, Found, Varies, StackAlongAlpha

    IF (.NOT. GetLogical(CompParams, 'Foil Sheet Sigma Varies', Found)) RETURN

    tfoil = GetConstReal(CompParams, 'Foil Thickness', FoundT)
    ff    = GetConstReal(CompParams, 'Fill Factor', FoundF)
    sgm   = FoilSheetBlockConductivity(CompParams, FoundS, Varies)
    IF (.NOT. (FoundT .AND. FoundF .AND. FoundS)) RETURN
    StackAlongAlpha = GetLogical(CompParams, 'Foil Sheet Stack Along Alpha', Found)
    IF (.NOT. Found) StackAlongAlpha = .TRUE.
    nSub = GetInteger(CompParams, 'Foil Sheet Sublayers', Found)
    IF (.NOT. Found .OR. nSub < 1) nSub = 1

    WRITE(Message,'(A,ES12.5)') 'Foil sheet mean sheet conductivity = ', sgm
    CALL Info('Circuits_Init', Message, Level=3)
    CALL SetFoilSheetHarmonicMaterial(CompParams, tfoil, ff, sgm, StackAlongAlpha, nSub)
!------------------------------------------------------------------------------
  END SUBROUTINE UpdateFoilSheetMaterial
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Derive everything the foil sheet needs from the three physical keywords
!> 'Foil Thickness', 'Fill Factor' and 'Sheet Conductivity', so that the SIF
!> writer supplies physics and Elmer owns the homogenization formulas. tau0 and
!> the DC sheet conductivity are stored for the transient ladders. In harmonic
!> mode the complex sheet conductivity and in-plane reluctivity are filled in
!> here, but only where the SIF did not give them explicitly: an explicit
!> 'Sigma 33' or 'Nu 22' always wins, which keeps hand written SIFs valid.
!------------------------------------------------------------------------------
  SUBROUTINE InitFoilSheetMaterial(CompParams, StackAlongAlpha, nSub)
!------------------------------------------------------------------------------
    USE CircuitUtils
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: CompParams
    LOGICAL :: StackAlongAlpha
    INTEGER :: nSub
    REAL(KIND=dp) :: tfoil, ff, sgm, tau0, sdc, mu0
    REAL(KIND=dp) :: nuinf, rr(6), tt(6), arr(6,1)
    INTEGER :: nlad, nlk, dStack, dPlane(2), d
    LOGICAL :: FoundT, FoundF, FoundS, HavePhys, Homog, Found, Transient
    LOGICAL :: SigmaVaries

    mu0 = 4.0d-7 * PI
    CALL FoilSheetNuDirections(StackAlongAlpha, dStack, dPlane)
    tfoil = GetConstReal(CompParams, 'Foil Thickness', FoundT)
    ff    = GetConstReal(CompParams, 'Fill Factor', FoundF)
    sgm   = FoilSheetBlockConductivity(CompParams, FoundS, SigmaVaries)
    HavePhys = FoundT .AND. FoundF .AND. FoundS

    ! Who owns the derived material has to be decided here, before anything is
    ! written: an explicit 'Sigma 33' or 'Nu 33' in the SIF always wins, and
    ! after the first derivation both keywords are present whoever wrote them.
    ! 'Nu 33' is the marker because Gamma is in the turn plane for either
    ! stacking direction, so it is derived in both cases.
    CALL ListAddLogical(CompParams, 'Foil Sheet Derive Sigma 33', &
        HavePhys .AND. .NOT. ListCheckPresent(CompParams,'Sigma 33'))
    CALL ListAddLogical(CompParams, 'Foil Sheet Derive Nu', &
        HavePhys .AND. .NOT. ListCheckPresent(CompParams,'Nu 33'))
    ! Only worth refreshing when the derivation actually uses it.
    CALL ListAddLogical(CompParams, 'Foil Sheet Sigma Varies', SigmaVaries .AND. HavePhys)

    IF (HavePhys) THEN
      IF (tfoil <= 0._dp) CALL Fatal('Circuits_Init','Foil sheet: "Foil Thickness" must be positive!')
      IF (ff <= 0._dp .OR. ff > 1._dp) CALL Fatal('Circuits_Init','Foil sheet: "Fill Factor" must be in (0,1]!')
      IF (sgm <= 0._dp) CALL Fatal('Circuits_Init','Foil sheet: "Sheet Conductivity" must be positive!')
      tau0 = mu0 * sgm * tfoil**2 / 4._dp
      ! SetFoilSheetHarmonicMaterial writes the values the run actually uses,
      ! with the sub-layer thickness and, for copper bands, no ff smearing.
      ! These are the one-band-per-turn defaults the transient ladder needs.
      sdc  = ff * sgm
      CALL ListAddConstReal(CompParams, 'Foil Sheet Tau0', tau0)
      CALL ListAddConstReal(CompParams, 'Foil Sheet Sigma DC', sdc)
      WRITE(Message,'(A,ES12.5,A,ES12.5,A,ES12.5)') 'Foil sheet physics: t = ', tfoil, &
          ', ff = ', ff, ', sigma = ', sgm
      CALL Info('Circuits_Init', Message, Level=3)
      IF (nSub > 1) THEN
        WRITE(Message,'(A,I0,A,ES12.5,A)') 'Foil sheet sub-layers: ', nSub, &
            ' per turn, the sheet coefficients use t/nSub = ', tfoil / nSub, ' m'
        CALL Info('Circuits_Init', Message, Level=3)
      END IF
      IF (nSub == 1) THEN
        WRITE(Message,'(A,ES12.5,A,ES12.5)') 'Foil sheet derived: tau0 = ', tau0, &
            ' s, DC sheet conductivity = ', sdc
        CALL Info('Circuits_Init', Message, Level=3)
      END IF
    END IF

    Transient = (TRIM(ListGetString(CurrentModel % Simulation,'Simulation Type',Found)) == 'transient')

    IF (Transient) THEN
      IF (nSub /= 1) CALL Fatal('Circuits_Init', &
          'Foil sheet transient supports "Sheet Sublayers = 1" only; the strand '// &
          'skin ladder and the reluctivity ladder are derived for one layer per turn.')
      IF (.NOT. HavePhys) CALL Fatal('Circuits_Init', &
          'Foil sheet transient needs "Foil Thickness", "Fill Factor" and "Sheet Conductivity"!')

      ! In transient the in-plane reluctivity becomes a Foster ladder. Derive it
      ! here from the same physics the harmonic closed form uses and publish it
      ! as the ordinary ladder keywords, so the AV solver reads one uniform form
      ! whoever wrote it. An explicit 'Nu 22 Residues' in the SIF wins.
      ! Post processing wants a conductivity; in transient that is the DC one,
      ! since the skin correction is carried by the strand ladder.
      IF (.NOT. ListCheckPresent(CompParams,'Sigma 33')) &
          CALL ListAddConstReal(CompParams, 'Sigma 33', sdc)

      Homog = GetLogical(CompParams, 'Homogenization Model', Found)
      IF (.NOT. Found) Homog = .FALSE.
      IF (Homog .AND. .NOT. ListCheckPresent(CompParams,'Nu 33 Residues')) THEN
        nlad = GetInteger(CompParams, 'Homogenization Ladder Order', Found)
        IF (.NOT. Found) nlad = 4
        CALL FoilSheetNuFoster(tau0, ff, nlad, nuinf, rr(1:nlad), tt(1:nlad))

        CALL ListAddConstReal(CompParams, FoilSheetNuKey(dStack)//' y0', 1._dp/mu0)
        ! Post processing still wants a plain reluctivity; give it the DC value
        ! of the same ladder, nu(0) = nu_inf + sum r_k, which for a non magnetic
        ! foil is 1/mu0 whatever the fill factor.
        CALL ListAddConstReal(CompParams, FoilSheetNuKey(dStack), 1._dp/mu0)
        DO d = 1, 2
          CALL ListAddConstReal(CompParams, FoilSheetNuKey(dPlane(d))//' y0', nuinf)
          arr = 0._dp
          arr(1:nlad,1) = rr(1:nlad)
          CALL ListAddConstRealArray(CompParams, FoilSheetNuKey(dPlane(d))//' Residues', &
              nlad, 1, arr(1:nlad,1:1))
          arr(1:nlad,1) = tt(1:nlad)
          CALL ListAddConstRealArray(CompParams, FoilSheetNuKey(dPlane(d))//' Taus', &
              nlad, 1, arr(1:nlad,1:1))
          CALL ListAddConstReal(CompParams, FoilSheetNuKey(dPlane(d)), nuinf + SUM(rr(1:nlad)))
        END DO

        WRITE(Message,'(A,I0,A,ES13.6,A,ES13.6)') 'Foil sheet Nu ladder: order ', nlad, &
            ', nu_inf = ', nuinf, ', nu(0) = ', nuinf + SUM(rr(1:nlad))
        CALL Info('Circuits_Init', Message, Level=5)
        DO nlk = 1, nlad
          WRITE(Message,'(A,I0,A,ES13.6,A,ES13.6,A)') '  pole ', nlk, ': r = ', rr(nlk), &
              ', T = ', tt(nlk), ' s'
          CALL Info('Circuits_Init', Message, Level=5)
        END DO
      END IF
      RETURN
    END IF

    IF (.NOT. HavePhys) RETURN
    CALL SetFoilSheetHarmonicMaterial(CompParams, tfoil, ff, sgm, StackAlongAlpha, nSub)
!------------------------------------------------------------------------------
  END SUBROUTINE InitFoilSheetMaterial
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE InitFoilSheetComponent(Comp, CompParams, CompInd, ExtMaster)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Component_t), POINTER :: Comp
    TYPE(ValueList_t), POINTER :: CompParams
    INTEGER :: CompInd, ExtMaster
    INTEGER :: nfoils, nLayers
    REAL(KIND=dp) :: tfoil, sgm, elemH, tsub
    LOGICAL :: Found, FoundS, Varies
    CHARACTER(LEN=MAX_NAME_LEN) :: str

    IF (CoordinateSystemDimension() /= 3) &
        CALL Fatal('Circuits_Init','Foil sheet coil type is implemented only in 3D!')

    Comp % nofturns = GetConstReal(CompParams, 'Number of Turns', Found)
    IF (.NOT. Found) CALL Fatal('Circuits_Init','Foil sheet: Number of Turns not found!')
    nfoils = NINT(Comp % nofturns)
    IF (nfoils < 1 .OR. ABS(Comp % nofturns - nfoils) > 1.0d-8) &
        CALL Fatal('Circuits_Init','Foil sheet: Number of Turns must be a positive integer!')

    ! Foils and flatwise flat wire stack along Alpha, edgewise flat wire along
    ! Beta. The default reproduces the foil layout, so existing inputs are
    ! untouched. The flat wire type defaults the same keyword to beta, which is
    ! why it is spelled out in both places rather than shared.
    str = GetString(CompParams, 'Stacking Direction', Found)
    IF (.NOT. Found) str = 'alpha'
    SELECT CASE (str)
    CASE ('alpha')
      Comp % StackAlongAlpha = .TRUE.
    CASE ('beta')
      Comp % StackAlongAlpha = .FALSE.
    CASE DEFAULT
      CALL Fatal('Circuits_Init','Foil sheet: Stacking Direction must be alpha or beta!')
    END SELECT
    CALL ListAddLogical(CompParams, 'Foil Sheet Stack Along Alpha', Comp % StackAlongAlpha)

    ! 0 asks the kernel to pick the layout from the block geometry and the mesh.
    Comp % nCells = GetInteger(CompParams, 'Sheet Cells', Found)
    IF (.NOT. Found) Comp % nCells = nfoils
    Comp % nSegments = GetInteger(CompParams, 'Sheet Segments', Found)
    IF (.NOT. Found) Comp % nSegments = 16

    IF (.NOT. ASSOCIATED(VariableGet(CurrentModel % Mesh % Variables, 'Alpha'))) &
        CALL Fatal('Circuits_Init','Foil sheet needs the direction field "Alpha"!')
    IF (.NOT. ASSOCIATED(VariableGet(CurrentModel % Mesh % Variables, 'Beta'))) &
        CALL Fatal('Circuits_Init','Foil sheet needs the direction field "Beta"!')

    ! 0 asks for the sub-layer count to be chosen from the skin depth; it needs
    ! the turn thickness, so it is resolved once the geometry is known.
    Comp % nSublayers = GetInteger(CompParams, 'Sheet Sublayers', Found)
    IF (.NOT. Found) Comp % nSublayers = 1

    ! The block geometry is also what a missing 'Foil Thickness' is derived
    ! from, so measure it whenever anything is left to the kernel. This has to
    ! run before InitFoilSheetMaterial, which turns the thickness into the
    ! homogenization coefficients.
    IF (Comp % nCells <= 0 .OR. Comp % nSegments <= 0 .OR. Comp % nSublayers <= 0 .OR. &
        .NOT. ListCheckPresent(CompParams, 'Foil Thickness')) &
        CALL FoilSheetAutoLayout(Comp, CompParams, nfoils)

    IF (Comp % nCells < 1 .OR. MOD(nfoils, Comp % nCells) /= 0) &
        CALL Fatal('Circuits_Init','Foil sheet: Number of Turns must be divisible by Sheet Cells!')
    Comp % foilsPerCell = nfoils / Comp % nCells

    IF (Comp % nSegments < 1) &
        CALL Fatal('Circuits_Init','Foil sheet: Sheet Segments must be positive!')

    IF (Comp % nSublayers <= 0) THEN
      IF (Comp % foilsPerCell > 1) THEN
        ! The copper only bands of a sub-layered cell are the conductor of ONE
        ! turn: they sit in the central fill factor of the cell and are spaced
        ! t/m. A cell that lumps several turns has that copper in several
        ! separate turns instead, with insulation between them, so the bands
        ! would be laid over the wrong geometry.
        Comp % nSublayers = 1
        CALL Info('Circuits_Init', &
            'Foil sheet: cells lump several turns, sub-layers off', Level=3)
      ELSE
        tfoil = GetConstReal(CompParams, 'Foil Thickness', Found)
        sgm = FoilSheetBlockConductivity(CompParams, FoundS, Varies)
        IF (.NOT. (Found .AND. FoundS)) CALL Fatal('Circuits_Init', &
            'Foil sheet: "Sheet Sublayers = 0" needs the thickness and "Sheet Conductivity"!')
        elemH = GetConstReal(CompParams, 'Foil Sheet Element Size', Found)
        IF (.NOT. Found) elemH = 0._dp
        Comp % nSublayers = FoilSheetAutoSublayers(tfoil, sgm, elemH)
      END IF
    END IF
    IF (Comp % nSublayers < 1) &
        CALL Fatal('Circuits_Init','Foil sheet: Sheet Sublayers must be positive!')
    IF (Comp % nSublayers > 1 .AND. Comp % foilsPerCell > 1) &
        CALL Fatal('Circuits_Init','Foil sheet: cells lump several turns, sub-layers off; '// &
            'raise "Sheet Cells" to one cell per turn or set "Sheet Sublayers = 1"!')
    CALL ListAddInteger(CompParams, 'Foil Sheet Sublayers', Comp % nSublayers)

    Comp % FillFactor = GetConstReal(CompParams, 'Fill Factor', Found)
    IF (.NOT. Found) Comp % FillFactor = 1._dp
    IF (Comp % nSublayers > 1 .AND. (Comp % FillFactor <= 0._dp .OR. Comp % FillFactor > 1._dp)) &
        CALL Fatal('Circuits_Init', &
            'Foil sheet: "Sheet Sublayers" > 1 needs "Fill Factor" in (0,1]!')

    CALL InitFoilSheetMaterial(CompParams, Comp % StackAlongAlpha, Comp % nSublayers)

    ! Sub-layers subdivide the stack further but stay one conductor: they share
    ! the voltage dof and the current constraint of their turn cell, so only the
    ! strand dofs multiply.
    nLayers = Comp % nCells * Comp % nSublayers

    ! dofs: V, V_1..V_nCells, then one current dof per (sub-layer, segment)
    Comp % ivar % dofs = 1
    Comp % ivar % pdofs = 0
    Comp % vvar % dofs = 1 + Comp % nCells + nLayers * Comp % nSegments
    Comp % vvar % pdofs = Comp % nCells + nLayers * Comp % nSegments

    Comp % coilthickness = 1._dp

    Comp % ElArea = GetConstReal(CompParams, 'Electrode Area', Found)
    IF (.NOT. Found) THEN
      CALL ComputeElectrodeArea(Comp, CompParams, ExtMaster)
      WRITE(Message,'(A,ES12.5)') 'Component '//I2S(CompInd)//' "Electrode Area" is ',Comp % ElArea
      CALL Info('Circuits_Init',Message,Level=10)
    END IF
    Comp % N_j = Comp % nofturns / Comp % ElArea

    IF (ALLOCATED(Comp % StrandWeight)) DEALLOCATE(Comp % StrandWeight)
    ALLOCATE(Comp % StrandWeight(nLayers * Comp % nSegments))
    Comp % StrandWeight = 0._dp

    ! Scale of the strand dofs: with c_kj = SigmaRef * y_kj the unknown y_kj is a
    ! voltage (y_kj = f V_k at DC), so the circuit block is as well conditioned
    ! as the flat wire one. Without it the strand dofs are sigma times larger
    ! than the cell voltages and the coupled solve stagnates around 1e-6.
    ! The DC sheet conductivity scales the strand dofs. In transient it comes
    ! from the physics keywords; 'Sigma 33' is not read there.
    Comp % SigmaRef = GetConstReal(CompParams, 'Foil Sheet Sigma DC', Found)
    IF (.NOT. Found) Comp % SigmaRef = GetConstReal(CompParams, 'Sigma 33', Found)
    IF (.NOT. Found .OR. Comp % SigmaRef <= 0._dp) Comp % SigmaRef = 1._dp
    CALL ListAddConstReal(CompParams, 'Foil Sheet Sigma Ref', Comp % SigmaRef)

    CALL ComputeFoilSheetSign(Comp, CompParams)

    CALL ListAddInteger(CompParams, 'Foil Sheet Cells', Comp % nCells)
    CALL ListAddInteger(CompParams, 'Foil Sheet Segments', Comp % nSegments)
    CALL ListAddInteger(CompParams, 'Foil Sheet Foils Per Cell', Comp % foilsPerCell)

    WRITE(Message,'(A,I0,A,I0,A,I0,A,I0,A)') 'Component '//I2S(CompInd)//' foil sheet: ', &
        Comp % nCells,' cells x ',Comp % nSublayers,' sub-layers x ',Comp % nSegments, &
        ' segments (',Comp % foilsPerCell,' turns per cell)'
    CALL Info('Circuits_Init',Message,Level=6)

    ! A sub-layer thinner than an element is still integrated exactly, because
    ! FoilSheetPieces clips it, but the field it drives cannot be resolved that
    ! finely. Say so rather than let it look like a converged layout.
    elemH = GetConstReal(CompParams, 'Foil Sheet Element Size', Found)
    IF (Found .AND. elemH > 0._dp) THEN
      tsub = GetConstReal(CompParams, 'Foil Sheet Stack Extent', Found) / MAX(1, nLayers)
      IF (Found .AND. tsub < elemH) THEN
        WRITE(Message,'(A,ES10.3,A,ES10.3,A)') 'Foil sheet sub-layer pitch ', tsub, &
            ' m is below the mean element size ', elemH, &
            ' m: the strand currents resolve the turn, the field does not'
        CALL Info('Circuits_Init', Message, Level=4)
      END IF
    END IF
    IF (Comp % StackAlongAlpha) THEN
      CALL Info('Circuits_Init','Component '//I2S(CompInd)//' foil sheet stacks along Alpha',Level=6)
    ELSE
      CALL Info('Circuits_Init','Component '//I2S(CompInd)//' foil sheet stacks along Beta',Level=6)
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE InitFoilSheetComponent
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Orientation of the Euler potential strand direction: grad(stack) x
!> grad(across) points along grad(W) or against it depending on how Alpha and
!> Beta are laid out and on which of the two is the stacking direction, so fix
!> the sign once per component from int grad(W) . (gStack x gAcross) dV.
!------------------------------------------------------------------------------
!> Choose the strand layout from the block geometry and the mesh, for
!> 'Sheet Cells' or 'Sheet Segments' given as 0, and derive a missing
!> 'Foil Thickness'. The SIF writer knows the coil but not the element size,
!> while the kernel can measure both: the stacking and across fields run from 0
!> to 1 over the block, so the block extents are the inverses of the mean
!> magnitudes of their gradients, and the mean element size is the cube root of
!> the mean element volume. A strand narrower than about one and a half elements
!> cannot be resolved, which sets the cap.
!------------------------------------------------------------------------------
  SUBROUTINE FoilSheetAutoLayout(Comp, CompParams, nfoils)
!------------------------------------------------------------------------------
    USE CircuitUtils
    IMPLICIT NONE
    TYPE(Component_t), POINTER :: Comp
    TYPE(ValueList_t), POINTER :: CompParams
    INTEGER :: nfoils

    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), sStack(:), sAcross(:)
    REAL(KIND=dp) :: detJ, ga(3), gb(3), wgp, vol, sga, sgb, sh, ve, blkT, blkH, elemH
    REAL(KIND=dp) :: ff, tfoil
    INTEGER :: e, n, gp, nmax, nel, nCellAuto, nSegAuto
    LOGICAL :: stat, Found

    nmax = CurrentModel % Mesh % MaxElementNodes
    ALLOCATE(Basis(nmax), dBasisdx(nmax,3), sStack(nmax), sAcross(nmax))
    vol = 0._dp; sga = 0._dp; sgb = 0._dp; sh = 0._dp; nel = 0

    DO e = 1, GetNOFActive()
      Element => GetActiveElement(e)
      IF (.NOT. ASSOCIATED(GetComponentParams(Element), CompParams)) CYCLE
      n = GetElementNOFNodes(Element)
      CALL GetElementNodes(Nodes, Element)
      CALL GetFlatWireLocalFields(Comp % StackAlongAlpha, Element, n, sStack, sAcross)
      IP = GaussPoints(Element)
      ve = 0._dp
      DO gp = 1, IP % n
        stat = ElementInfo(Element, Nodes, IP % U(gp), IP % V(gp), IP % W(gp), &
            detJ, Basis, dBasisdx)
        wgp = IP % s(gp) * detJ
        ga = MATMUL(sStack(1:n), dBasisdx(1:n,:))
        gb = MATMUL(sAcross(1:n), dBasisdx(1:n,:))
        sga = sga + wgp * SQRT(SUM(ga*ga))
        sgb = sgb + wgp * SQRT(SUM(gb*gb))
        ve = ve + wgp
      END DO
      vol = vol + ve
      sh = sh + ve**(1._dp/3._dp)
      nel = nel + 1
    END DO

    vol = ParallelReduction(vol)
    sga = ParallelReduction(sga)
    sgb = ParallelReduction(sgb)
    sh  = ParallelReduction(sh)
    nel = ParallelReduction(nel)

    IF (vol <= 0._dp .OR. sga <= 0._dp .OR. sgb <= 0._dp .OR. nel < 1) &
        CALL Fatal('Circuits_Init','Foil sheet automatic layout: the coil block is empty!')

    blkT = vol / sga
    blkH = vol / sgb
    elemH = sh / nel

    ! The cap is on the TURN cells only. Sub-layers subdivide a turn further on
    ! purpose, so they are never traded away against the mesh here.
    IF (Comp % nCells <= 0) THEN
      nCellAuto = MIN(nfoils, MAX(1, FLOOR(blkT / (1.5_dp * elemH))))
      DO WHILE (nCellAuto > 1 .AND. MOD(nfoils, nCellAuto) /= 0)
        nCellAuto = nCellAuto - 1
      END DO
      Comp % nCells = nCellAuto
    END IF

    ! V_e^(1/3) is about half a tetrahedron's edge length, so one strand per
    ! V_e^(1/3) is about half an element edge. Strands narrower than an element
    ! are legitimate: FoilSheetPieces clips them exactly, so each still carries
    ! its own volume and its own current. Resolving the across direction matters
    ! for an edgewise flat wire, whose across direction is the wide side of the
    ! turn and carries the current crowding into the inner edge; for a foil
    ! stack, which has no such gradient across the width, it is nearly free.
    IF (Comp % nSegments <= 0) THEN
      nSegAuto = NINT(blkH / elemH)
      Comp % nSegments = MIN(40, MAX(4, nSegAuto))
    END IF

    ! 'Foil Thickness' is optional: a winding that gives its fill factor and its
    ! turn count already fixes the conductor thickness through the stack extent,
    ! t = ff * L_stack / N. An explicit keyword always wins.
    IF (.NOT. ListCheckPresent(CompParams, 'Foil Thickness')) THEN
      ff = GetConstReal(CompParams, 'Fill Factor', Found)
      IF (Found .AND. ff > 0._dp) THEN
        tfoil = ff * blkT / nfoils
        CALL ListAddConstReal(CompParams, 'Foil Thickness', tfoil)
        WRITE(Message,'(A,ES12.5,A,ES11.4,A,F7.4,A,I0)') &
            'Foil sheet derived "Foil Thickness" = ', tfoil, ' m from stack extent ', blkT, &
            ', fill factor ', ff, ', turns ', nfoils
        CALL Info('Circuits_Init', Message, Level=3)
      END IF
    END IF

    CALL ListAddConstReal(CompParams, 'Foil Sheet Stack Extent', blkT)
    CALL ListAddConstReal(CompParams, 'Foil Sheet Element Size', elemH)
    WRITE(Message,'(A,ES11.4,A,ES11.4,A,ES11.4)') 'Foil sheet block: stack extent ', blkT, &
        ', across extent ', blkH, ', mean element size ', elemH
    CALL Info('Circuits_Init', Message, Level=3)
    WRITE(Message,'(A,I0,A,I0,A,I0,A)') 'Foil sheet automatic layout: Sheet Cells = ', &
        Comp % nCells, ', Sheet Segments = ', Comp % nSegments, ' (', nel, ' block elements)'
    CALL Info('Circuits_Init', Message, Level=3)

    DEALLOCATE(Basis, dBasisdx, sStack, sAcross)
!------------------------------------------------------------------------------
  END SUBROUTINE FoilSheetAutoLayout
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Mean circulation of a closed coil's direction field over its current loops.
! With t = grad(Alpha) x grad(Beta) the volume integral of u.t is the
! double integral of the loop circulation of u over (Alpha, Beta), both of which
! span [0,1], so it is the mean circulation: exactly 1 for the electrode
! potential W of an open coil. The CoilSolver cut potential of a closed coil is
! not: its Dirichlet cut pins every node of the element layer it crosses, so the
! potential ramps over the loop minus that layer and the field is too large by
! that angular fraction. Measure it once and let GetCoilWBase divide it out.
!------------------------------------------------------------------------------
  SUBROUTINE ComputeCoilCirculation(CompParams, CompInd)
!------------------------------------------------------------------------------
    USE CircuitUtils
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: CompParams
    INTEGER :: CompInd
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), Wloc(:), Aloc(:), Bloc(:)
    REAL(KIND=dp) :: detJ, gw(3), ga(3), gb(3), tv(3), circ
    INTEGER :: e, n, gp, nmax
    LOGICAL :: stat
    REAL(KIND=dp), PARAMETER :: MinCirculation = 0.5_dp, MaxCirculation = 2.0_dp
!------------------------------------------------------------------------------

    ! Measuring it again on the already normalized field would just return 1.
    IF (ListCheckPresent(CompParams, 'Coil Circulation')) RETURN

    nmax = CurrentModel % Mesh % MaxElementNodes
    ALLOCATE(Basis(nmax), dBasisdx(nmax,3), Wloc(nmax), Aloc(nmax), Bloc(nmax))

    circ = 0._dp
    DO e = 1, GetNOFActive()
      Element => GetActiveElement(e)
      IF (.NOT. ASSOCIATED(GetComponentParams(Element), CompParams)) CYCLE
      n = GetElementNOFNodes(Element)
      CALL GetElementNodes(Nodes, Element)
      CALL GetCoilWBase(Element, n, CompParams, Wloc)
      CALL GetScalarLocalSolution(Aloc, 'Alpha', UElement=Element)
      CALL GetScalarLocalSolution(Bloc, 'Beta', UElement=Element)
      IP = GaussPoints(Element)
      DO gp = 1, IP % n
        stat = ElementInfo(Element, Nodes, IP % U(gp), IP % V(gp), IP % W(gp), &
            detJ, Basis, dBasisdx)
        gw = MATMUL(Wloc(1:n), dBasisdx(1:n,:))
        ga = MATMUL(Aloc(1:n), dBasisdx(1:n,:))
        gb = MATMUL(Bloc(1:n), dBasisdx(1:n,:))
        tv = CrossProduct(ga, gb)
        circ = circ + IP % s(gp) * detJ * SUM(gw*tv)
      END DO
    END DO

    ! The sign is the coil's winding sense, which ComputeFoilSheetSign owns.
    circ = ABS(ParallelReduction(circ))

    WRITE(Message,'(A,ES13.6)') 'Component '//I2S(CompInd)// &
        ' closed coil mean circulation: ', circ
    CALL Info('Circuits_Init', Message, Level=5)

    IF (circ < MinCirculation .OR. circ > MaxCirculation) THEN
      WRITE(Message,'(A,ES13.6,A)') 'Component '//I2S(CompInd)// &
          ' closed coil circulation is ', circ, ', expected near 1. The CoilSolver '// &
          'cut potential is not a valid direction field: check "Coil Normal", '// &
          '"Coil Center" and that the coil really is a closed loop.'
      CALL Fatal('Circuits_Init', Message)
    END IF

    CALL ListAddConstReal(CompParams, 'Coil Circulation', circ)

    DEALLOCATE(Basis, dBasisdx, Wloc, Aloc, Bloc)
!------------------------------------------------------------------------------
  END SUBROUTINE ComputeCoilCirculation
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ComputeFoilSheetSign(Comp, CompParams)
!------------------------------------------------------------------------------
    USE CircuitUtils
    IMPLICIT NONE
    TYPE(Component_t), POINTER :: Comp
    TYPE(ValueList_t), POINTER :: CompParams
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), sStack(:), sAcross(:), Wloc(:)
    REAL(KIND=dp), ALLOCATABLE :: dNode(:), aNode(:)
    LOGICAL, ALLOCATABLE :: inBlk(:), skipn(:)
    REAL(KIND=dp) :: detJ, gw(3), tv(3), flux, sgn
    INTEGER :: e, n, gp, nmax, nno, i, v, gnode, ind, nPiece
    INTEGER :: pCell(MaxFoilSheetPieces), pSeg(MaxFoilSheetPieces)
    REAL(KIND=dp) :: pVol(MaxFoilSheetPieces), pBary(4,MaxFoilSheetPieces)
    REAL(KIND=dp) :: volerr, dmax, amax, cval, contr, Vpar
    REAL(KIND=dp) :: scov, tcov
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: stat

    Mesh => CurrentModel % Mesh
    nmax = Mesh % MaxElementNodes
    nno = Mesh % NumberOfNodes
    ALLOCATE(Basis(nmax), dBasisdx(nmax,3), sStack(nmax), sAcross(nmax), Wloc(nmax))
    ALLOCATE(dNode(nno), aNode(nno), inBlk(nno), skipn(nno))
    dNode = 0._dp; aNode = 0._dp; inBlk = .FALSE.; skipn = .FALSE.
    flux = 0._dp
    volerr = 0._dp
    scov = 0._dp
    tcov = 0._dp

    ! A node is interior to the coil block when every element that carries its
    ! basis function is a block element, it is not on a mesh boundary, and it is
    ! not shared with another partition. Only there must the discrete divergence
    ! of the strand source vanish; on the electrodes and the block surface it
    ! carries the terminal current and is nonzero by design.
    DO e = 1, Mesh % NumberOfBulkElements
      Element => Mesh % Elements(e)
      n = Element % TYPE % NumberOfNodes
      IF (ASSOCIATED(GetComponentParams(Element), CompParams)) THEN
        inBlk(Element % NodeIndexes(1:n)) = .TRUE.
      ELSE
        skipn(Element % NodeIndexes(1:n)) = .TRUE.
      END IF
    END DO
    DO e = Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(e)
      n = Element % TYPE % NumberOfNodes
      skipn(Element % NodeIndexes(1:n)) = .TRUE.
    END DO
    IF (ASSOCIATED(Mesh % ParallelInfo % GInterface)) THEN
      WHERE (Mesh % ParallelInfo % GInterface(1:nno)) skipn = .TRUE.
    END IF

    DO e = 1, GetNOFActive()
      Element => GetActiveElement(e)
      IF (.NOT. ASSOCIATED(GetComponentParams(Element), CompParams)) CYCLE
      n = GetElementNOFNodes(Element)
      CALL GetElementNodes(Nodes, Element)
      CALL GetFlatWireLocalFields(Comp % StackAlongAlpha, Element, n, sStack, sAcross)
      CALL GetCoilWBase(Element, n, CompParams, Wloc)
      IP = GaussPoints(Element)
      DO gp = 1, IP % n
        stat = ElementInfo(Element, Nodes, IP % U(gp), IP % V(gp), IP % W(gp), &
            detJ, Basis, dBasisdx)
        gw = MATMUL(Wloc(1:n), dBasisdx(1:n,:))
        tv = FoilSheetDirection(sStack, sAcross, dBasisdx, n, 1._dp)
        flux = flux + IP % s(gp) * detJ * SUM(gw*tv)
      END DO

      ! Clipping checks (a) and (c): the pieces must fill the element exactly,
      ! and a strand-wise constant source built on them must have zero discrete
      ! divergence at every interior node, for an arbitrary set of strand
      ! currents. The currents are a fixed function of the strand index so that
      ! every partition uses the same ones.
      IF (n /= 4 .OR. Element % TYPE % ElementCode /= 504) CYCLE
      CALL FoilSheetPieces(Comp % nCells, Comp % nSegments, sStack(1:4), sAcross(1:4), &
          nPiece, pCell, pSeg, pVol, pBary, Comp % nSublayers, Comp % FillFactor)
      ! One band per turn tiles the element; copper only bands cover the fill
      ! factor of it, so the global coverage is the check there.
      IF (Comp % nSublayers == 1) volerr = MAX(volerr, ABS(SUM(pVol(1:nPiece)) - 1._dp))
      scov = scov + SUM(pVol(1:nPiece)) * detJ
      tcov = tcov + detJ
      stat = ElementInfo(Element, Nodes, 0.25_dp, 0.25_dp, 0.25_dp, detJ, Basis, dBasisdx)
      tv = FoilSheetDirection(sStack, sAcross, dBasisdx, n, 1._dp)
      Vpar = detJ / 6._dp
      DO i = 1, nPiece
        ind = (pCell(i)-1) * Comp % nSegments + pSeg(i)  ! pCell is the sub-layer here
        cval = SIN(1.234_dp*ind) + 0.5_dp*COS(2.345_dp*ind)
        DO v = 1, 4
          gnode = Element % NodeIndexes(v)
          contr = cval * SUM(tv*dBasisdx(v,:)) * pVol(i) * Vpar
          dNode(gnode) = dNode(gnode) + contr
          aNode(gnode) = aNode(gnode) + ABS(contr)
        END DO
      END DO
    END DO

    dmax = 0._dp; amax = 0._dp
    DO i = 1, nno
      IF (.NOT. inBlk(i) .OR. skipn(i)) CYCLE
      dmax = MAX(dmax, ABS(dNode(i)))
      amax = MAX(amax, aNode(i))
    END DO
    dmax = ParallelReduction(dmax, 2)
    amax = ParallelReduction(amax, 2)
    volerr = ParallelReduction(volerr, 2)

    flux = ParallelReduction(flux)
    sgn = 1._dp
    IF (flux < 0._dp) sgn = -1._dp
    CALL ListAddConstReal(CompParams, 'Foil Sheet Direction Sign', sgn)
    WRITE(Message,'(A,F5.1,A,ES12.5)') 'Foil sheet strand direction sign ', sgn, &
        ', int gradW . (gA x gB) = ', flux
    CALL Info('Circuits_Init', Message, Level=5)

    scov = ParallelReduction(scov)
    tcov = ParallelReduction(tcov)
    IF (Comp % nSublayers == 1) THEN
      WRITE(Message,'(A,ES10.3)') 'Foil sheet clipping, max |sum(piece vol)/V - 1| over block: ', volerr
    ELSE
      WRITE(Message,'(A,F10.7,A,F10.7)') 'Foil sheet clipping, strand volume / block volume: ', &
          scov / MAX(tcov, TINY(tcov)), ', fill factor ', Comp % FillFactor
    END IF
    CALL Info('Circuits_Init', Message, Level=5)
    WRITE(Message,'(A,ES10.3,A,ES10.3,A,ES10.3)') &
        'Foil sheet source divergence at interior nodes: max ', dmax, &
        ' of scale ', amax, ', relative ', dmax / MAX(amax, TINY(amax))
    CALL Info('Circuits_Init', Message, Level=5)

    DEALLOCATE(Basis, dBasisdx, sStack, sAcross, Wloc, dNode, aNode, inBlk, skipn)
!------------------------------------------------------------------------------
  END SUBROUTINE ComputeFoilSheetSign
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Geometric fill factor f = nStack*t_c/L_stack [* nAcross*w_c/L_across]. The
!> block lengths follow from the normalized direction fields, L = V/int|grad s|.
!------------------------------------------------------------------------------
  SUBROUTINE ComputeFlatWireFillFactor(Comp, CompParams)
!------------------------------------------------------------------------------
    USE CircuitUtils
    IMPLICIT NONE
    TYPE(Component_t), POINTER :: Comp
    TYPE(ValueList_t), POINTER :: CompParams
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), sStack(:), sAcross(:)
    REAL(KIND=dp) :: Vol, GradStack, GradAcross, detJ, Lstack, Lacross, tc, wc, f
    INTEGER :: t, n, i, nmax
    LOGICAL :: stat, Found, Parallel

    nmax = CurrentModel % Mesh % MaxElementNodes
    ALLOCATE(Basis(nmax), dBasisdx(nmax,3), sStack(nmax), sAcross(nmax))

    Vol = 0._dp
    GradStack = 0._dp
    GradAcross = 0._dp
    DO t=1,GetNOFActive()
      Element => GetActiveElement(t)
      IF (.NOT. ElAssocToComp(Element, Comp)) CYCLE
      n = GetElementNOFNodes(Element)
      CALL GetElementNodes(Nodes, Element)
      CALL GetFlatWireLocalFields(Comp % StackAlongAlpha, Element, n, sStack, sAcross)
      IP = GaussPoints(Element)
      DO i=1,IP % n
        stat = ElementInfo(Element, Nodes, IP % U(i), IP % V(i), IP % W(i), detJ, Basis, dBasisdx)
        Vol = Vol + IP % s(i)*detJ
        GradStack = GradStack + IP % s(i)*detJ*SQRT(SUM(MATMUL(sStack(1:n), dBasisdx(1:n,:))**2))
        GradAcross = GradAcross + IP % s(i)*detJ*SQRT(SUM(MATMUL(sAcross(1:n), dBasisdx(1:n,:))**2))
      END DO
    END DO

    Parallel = (ParEnv % PEs > 1) .AND. .NOT. CurrentModel % Mesh % SingleMesh
    IF (Parallel) THEN
      Vol = ParallelReduction(Vol)
      GradStack = ParallelReduction(GradStack)
      GradAcross = ParallelReduction(GradAcross)
    END IF
    IF (GradStack <= 0._dp .OR. GradAcross <= 0._dp) &
        CALL Fatal('ComputeFlatWireFillFactor','Direction fields are constant over the component!')

    Lstack = Vol / GradStack
    Lacross = Vol / GradAcross

    tc = GetConstReal(CompParams, 'Conductor Thickness', Found)
    f = Comp % nStack * tc / Lstack
    wc = GetConstReal(CompParams, 'Conductor Width', Found)
    IF (Found) f = f * Comp % nAcross * wc / Lacross

    WRITE(Message,'(A,ES12.5,A,ES12.5)') 'Flat wire block length along stack ',Lstack,' across ',Lacross
    CALL Info('ComputeFlatWireFillFactor',Message,Level=6)
    IF (f <= 0._dp .OR. f > 1._dp) THEN
      WRITE(Message,'(A,ES12.5)') 'Flat wire fill factor out of (0,1]: ',f
      CALL Fatal('ComputeFlatWireFillFactor',Message)
    END IF
    Comp % FillFactor = f
!------------------------------------------------------------------------------
  END SUBROUTINE ComputeFlatWireFillFactor
!------------------------------------------------------------------------------

END MODULE CircuitsMod

MODULE CircMatInitMod

  USE CircuitsMod
  IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
   SUBROUTINE SetCircuitsParallelInfo()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Matrix_t), POINTER :: CM
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(Solver_t), POINTER :: ASolver
    TYPE(Circuit_t), POINTER :: Circuits(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: i, nm, Circuit_tot_n, p, j, &
               cnt(Parenv % PEs), r_cnt(ParEnv % PEs), &
               RowId, nn, l, k, n_Circuits
    
    CM => CurrentModel % CircuitMatrix
    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('SetCircuitsParallelInfo','ASolver not found!')
    nm = ASolver % Matrix % NumberOfRows
    Circuit_tot_n = CurrentModel % Circuit_tot_n
    Circuits => CurrentModel % Circuits
    n_Circuits = CurrentModel % n_Circuits
    
    IF(.NOT.ASSOCIATED(CM % ParallelInfo)) THEN
      ALLOCATE(CM % ParallelInfo)
      ALLOCATE(CM % ParallelInfo % NeighbourList(nm+Circuit_tot_n))
      DO i=1,nm+Circuit_tot_n
        CM % ParallelInfo % NeighbourList(i) % Neighbours => Null()
      END DO
    END IF

    DO p = 1,n_Circuits
      DO i=1,Circuits(p) % n
        cnt  = 0
        Cvar => Circuits(p) % CircuitVariables(i)
        IF(ASSOCIATED(CVar%Component)) THEN
          DO j=1,GetNOFACtive()
             Element => GetActiveElement(j)
               IF(ElAssocToCvar(Element, Cvar)) THEN
                 cnt(ParEnv % mype+1)=cnt(ParEnv % mype+1)+1
               END IF
          END DO
        END IF
        CALL MPI_ALLREDUCE(cnt,r_cnt,ParEnv % PEs, MPI_INTEGER, &
                MPI_MAX,ASolver % Matrix % Comm,j)


        RowId = Cvar % ValueId + nm

        nn = COUNT(r_cnt>0)
        IF( r_cnt(CVar % Owner+1)<=0 ) THEN
          r_cnt(CVar % Owner+1) = 1
          nn=nn + 1
        END IF

        IF(r_cnt(Parenv % myPE+1)<=0) THEN
          r_cnt(parenv % mype+1) = 1
          nn = nn + 1
        END IF

        r_cnt  = 1; nn=Parenv % PEs ! for now

        IF (Circuits(p) % Harmonic) THEN
          DO j=1,Cvar % Dofs
            IF(.NOT.ASSOCIATED(CM % ParallelInfo % NeighbourList(RowId+AddIndex(j-1))%Neighbours)) THEN
              ALLOCATE(CM % ParallelInfo % NeighbourList(RowId+AddIndex(j-1)) % Neighbours(nn))
              ALLOCATE(CM % ParallelInfo % NeighbourList(RowId+AddImIndex(j-1)) % Neighbours(nn))
            END IF
            CM % ParallelInfo % NeighbourList(RowId+AddIndex(j-1)) % Neighbours(1)   = CVar % Owner
            CM % ParallelInfo % NeighbourList(RowId+AddImIndex(j-1)) % Neighbours(1) = Cvar % Owner
            l = 1
            DO k=0,ParEnv % PEs-1
              IF(k==CVar % Owner) CYCLE
              IF(r_cnt(k+1)>0) THEN
                l = l + 1
                CM % ParallelInfo % NeighbourList(RowId+AddIndex(j-1)) % Neighbours(l) = k
                CM % ParallelInfo % NeighbourList(RowId+AddImIndex(j-1)) % Neighbours(l) = k
              END IF
            END DO
            CM % RowOwner(RowId + AddIndex(j-1))   = Cvar % Owner
            CM % RowOwner(RowId + AddImIndex(j-1)) = Cvar % Owner
          END DO
        ELSE
          DO j=1,Cvar % Dofs
            IF(.NOT.ASSOCIATED(CM % ParallelInfo % NeighbourList(RowId+j-1)%Neighbours)) THEN
              ALLOCATE(CM % ParallelInfo % NeighbourList(RowId+j-1) % Neighbours(nn))
            END IF
            CM % ParallelInfo % NeighbourList(RowId+j-1) % Neighbours(1) = CVar % Owner
            l = 1
            DO k=0,ParEnv % PEs-1
              IF(k==CVar % Owner) CYCLE
              IF(r_cnt(k+1)>0) THEN
                l = l + 1
                CM % ParallelInfo % NeighbourList(RowId+j-1) % Neighbours(l) = k
              END IF
            END DO
            CM % RowOwner(RowId+j-1) = Cvar % Owner
          END DO
        END IF
      END DO
    END DO


    IF ( parenv % mype==0 ) THEN
      DO i=1,parenv % pes
        CALL INFO('SetCircuitsParallelInfo','owners: '//i2s(i)//' '//i2s(count(cm % rowowner==i-1)), Level=9)
      END DO
    END IF
!------------------------------------------------------------------------------
   END SUBROUTINE SetCircuitsParallelInfo
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE CountCmplxMatElement(Rows, Cnts, RowId, dofs)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: Rows(:), Cnts(:)
    INTEGER :: RowId, dofs

    ! Matrix element structure:
    !
    ! Re -Im
    ! Im Re
    !
    ! First do Re -Im:
    ! ----------------
    Cnts(RowId) = Cnts(RowId) + 2 * dofs

    ! Then do Im Re:
    ! --------------
    Cnts(RowId+1) = Cnts(RowId+1) + 2 * dofs

!------------------------------------------------------------------------------
   END SUBROUTINE CountCmplxMatElement
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CountMatElement(Rows, Cnts, RowId, dofs, Harmonic)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: Rows(:), Cnts(:)
    INTEGER :: RowId, dofs
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm
    
    IF (.NOT. PRESENT(Harmonic)) THEN
      harm = CurrentModel % HarmonicCircuits
    ELSE
      harm = Harmonic
    END IF
    
    IF (harm) THEN
      CALL CountCmplxMatElement(Rows, Cnts, RowId, dofs)
    ELSE
      Cnts(RowId) = Cnts(RowId) + dofs
    END IF

!------------------------------------------------------------------------------
   END SUBROUTINE CountMatElement
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CreateCmplxMatElement(Rows, Cols, Cnts, RowId, ColId)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: Rows(:), Cols(:), Cnts(:)
    INTEGER :: RowId, ColId

    ! Matrix element structure:
    !
    ! Re -Im
    ! Im Re
    !
    ! First do Re (0,0):
    ! ------------------
    Cols(Rows(RowId) + Cnts(RowId)) = ColId
    Cnts(RowId) = Cnts(RowId) + 1

    ! Then do -Im (0,1):
    ! ------------------
    Cols(Rows(RowId) + Cnts(RowId)) = ColId + 1
    Cnts(RowId) = Cnts(RowId) + 1

    ! Then do Re (1,0):
    ! -----------------
    Cols(Rows(RowId+1) + Cnts(RowId+1)) = ColId
    Cnts(RowId+1) = Cnts(RowId+1) + 1

    ! Then do Im (1,1):
    ! -----------------
    Cols(Rows(RowId+1) + Cnts(RowId+1)) = ColId + 1
    Cnts(RowId+1) = Cnts(RowId+1) + 1
    
!------------------------------------------------------------------------------
   END SUBROUTINE CreateCmplxMatElement
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CreateMatElement(Rows, Cols, Cnts, RowId, ColId, Harmonic)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: Rows(:), Cols(:), Cnts(:)
    INTEGER :: RowId, ColId
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm
    
    IF (.NOT. PRESENT(Harmonic)) THEN
      harm = CurrentModel % HarmonicCircuits
    ELSE
      harm = Harmonic
    END IF
    
    IF (harm) THEN
      CALL CreateCmplxMatElement(Rows, Cols, Cnts, RowId, ColId)
    ELSE
      Cols(Rows(RowId) + Cnts(RowId)) = ColId
      Cnts(RowId) = Cnts(RowId) + 1
    END IF

!------------------------------------------------------------------------------
   END SUBROUTINE CreateMatElement
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE CountBasicCircuitEquations(Rows, Cnts)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuits(:)
    TYPE(CircuitVariable_t), POINTER :: Cvar
    INTEGER :: i, j, p, nm, RowId, n_Circuits
    INTEGER, POINTER :: Rows(:), Cnts(:)
    
    Circuits => CurrentModel % Circuits
    n_Circuits = CurrentModel % n_Circuits
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    
    ! Basic circuit equations...
    ! ---------------------------
    DO p = 1,n_Circuits
      DO i=1,Circuits(p) % n
        Cvar => Circuits(p) % CircuitVariables(i)
        IF(CVar % Owner /= ParEnv % myPE) CYCLE

        RowId = Cvar % ValueId + nm
        DO j=1,Circuits(p) % n
          IF(Cvar % A(j)/=0._dp.OR.Cvar % B(j)/=0._dp) &
             CALL CountMatElement(Rows, Cnts, RowId, 1)
        END DO
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountBasicCircuitEquations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CreateBasicCircuitEquations(Rows, Cols, Cnts)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuits(:)
    TYPE(CircuitVariable_t), POINTER :: Cvar
    INTEGER :: i, j, p, nm, RowId, ColId, n_Circuits
    INTEGER, POINTER :: Rows(:), Cols(:), Cnts(:)
    
    Circuits => CurrentModel % Circuits
    n_Circuits = CurrentModel % n_Circuits
    nm = CurrentModel % Asolver % Matrix % NumberOfRows
    
    ! Basic circuit equations...
    ! ---------------------------
    DO p = 1,n_Circuits
      DO i=1,Circuits(p) % n
        Cvar => Circuits(p) % CircuitVariables(i)
        IF(Cvar % Owner /= ParEnv % myPE) CYCLE

        RowId = Cvar % ValueId + nm
        DO j=1,Circuits(p) % n
          IF(Cvar % A(j)/=0._dp .OR. Cvar % B(j)/=0._dp) THEN
            ColId = Circuits(p) % CircuitVariables(j) % ValueId + nm
            CALL CreateMatElement(Rows, Cols, Cnts, RowId, ColId)
          END IF
        END DO
      END DO
    END DO

!------------------------------------------------------------------------------
   END SUBROUTINE CreateBasicCircuitEquations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CountComponentEquations(Rows, Cnts, Done, dofsdone)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuits(:)
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(Solver_t), POINTER :: ASolver
    TYPE(Element_t), POINTER :: Element
    TYPE(Component_t), POINTER :: Comp
    INTEGER :: i, j, p, nm, nn, nd, &
               RowId, ColId, n_Circuits, &
               CompInd, q
    INTEGER, POINTER :: Rows(:), Cnts(:)
    LOGICAL :: dofsdone
    LOGICAL*1 :: Done(:)
    
    Circuits => CurrentModel % Circuits
    n_Circuits = CurrentModel % n_Circuits
    Asolver => CurrentModel % Asolver
    nm = Asolver % Matrix % NumberOfRows
    DO p=1,n_Circuits
      DO CompInd=1,Circuits(p) % n_comp
        Done = .FALSE.
        Comp => Circuits(p) % Components(CompInd)
        Cvar => Comp % vvar
        RowId = Cvar % ValueId + nm
        ColId = Cvar % ValueId + nm
        IF (Comp % ComponentType == 'resistor') THEN
            CALL CountMatElement(Rows, Cnts, RowId, 1)
            CALL CountMatElement(Rows, Cnts, RowId, 1)
            CYCLE
        ELSE
          SELECT CASE (Comp % CoilType)
          CASE('stranded')
             CALL CountMatElement(Rows, Cnts, RowId, 1)
             CALL CountMatElement(Rows, Cnts, RowId, 1)
          CASE('massive')
             CALL CountMatElement(Rows, Cnts, RowId, 1)
             CALL CountMatElement(Rows, Cnts, RowId, 1)
          CASE('foil winding')
            ! V = V0 + V1*alpha + V2*alpha^2 + ...
            CALL CountMatElement(Rows, Cnts, RowId, Cvar % dofs)

            ! Circuit eqns for the pdofs:
            ! I(Vj) - I = 0
            ! ------------------------------------
            DO j=1, Cvar % pdofs
              CALL CountMatElement(Rows, Cnts, RowId + AddIndex(j), Cvar % dofs)
            END DO
          CASE('flat wire')
            ! V - sum_k V_k = 0
            CALL CountMatElement(Rows, Cnts, RowId, Cvar % dofs)
            ! cell rows: (V_k, V_k) and (V_k, I)
            DO j=1, Cvar % pdofs
              CALL CountMatElement(Rows, Cnts, RowId + AddIndex(j), 2)
            END DO
          CASE('foil sheet')
            ! V - sum_k m_k V_k = 0
            CALL CountMatElement(Rows, Cnts, RowId, 1 + Comp % nCells)
            DO j=1, Comp % nCells
              ! cell row: (V_k, I) and one (V_k, c) per strand of the turn
              CALL CountMatElement(Rows, Cnts, RowId + AddIndex(j), &
                  1 + Comp % nSublayers * Comp % nSegments)
            END DO
            DO j=Comp % nCells + 1, Cvar % pdofs
              ! strand row: (c_kj, c_kj) and (c_kj, V_k)
              CALL CountMatElement(Rows, Cnts, RowId + AddIndex(j), 2)
            END DO
          END SELECT
        END IF

!        temp = SUM(Cnts)
!print *, "Active elements", ParEnv % Mype, ":", GetNOFActive()
        DO q=GetNOFActive(),1,-1
          Element => GetActiveElement(q)
          CALL CountComponentElements(Element, Comp, RowId, Rows, Cnts, Done, dofsdone)
        END DO

        DO q=GetNOFBoundaryElements(),1,-1
          Element => GetBoundaryElement(q)
          CALL CountComponentElements(Element, Comp, RowId, Rows, Cnts, Done, dofsdone)
        END DO
!        Comp % nofcnts = SUM(Cnts) - temp
!        print *, ParEnv % Mype, "CompInd:", CompInd, "Comp % nofcnts", Comp % nofcnts
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountComponentEquations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CreateComponentEquations(Rows, Cols, Cnts, Done, dofsdone)
!------------------------------------------------------------------------------
    USE CircuitUtils
    IMPLICIT NONE
    TYPE(Circuit_t), POINTER :: Circuits(:)
    TYPE(CircuitVariable_t), POINTER :: Cvar
    TYPE(Solver_t), POINTER :: ASolver
    TYPE(Element_t), POINTER :: Element
    TYPE(Component_t), POINTER :: Comp
    INTEGER :: i, j, jj, ll, p, nm, nn, nd, &
               VvarId, IvarId, n_Circuits, &
               CompInd, q
    INTEGER, POINTER :: Rows(:), Cols(:), Cnts(:)
    LOGICAL :: dofsdone
    LOGICAL*1 :: Done(:)
    
    Circuits => CurrentModel % Circuits
    n_Circuits = CurrentModel % n_Circuits
    Asolver => CurrentModel % Asolver
    nm = Asolver % Matrix % NumberOfRows

    DO p=1,n_Circuits
      DO CompInd = 1, Circuits(p) % n_comp
        Done = .FALSE.
        Comp => Circuits(p) % Components(CompInd)
        Cvar => Comp % vvar
        VvarId = Comp % vvar % ValueId + nm
        IvarId = Comp % ivar % ValueId + nm

        IF (Comp % ComponentType == 'resistor') THEN
            CALL CreateMatElement(Rows, Cols, Cnts, VvarId, IvarId)
            CALL CreateMatElement(Rows, Cols, Cnts, VvarId, VvarId)
            CYCLE
        ELSE
          SELECT CASE (Comp % CoilType)
          CASE('stranded')
            CALL CreateMatElement(Rows, Cols, Cnts, VvarId, IvarId)
            CALL CreateMatElement(Rows, Cols, Cnts, VvarId, VvarId)
          CASE('massive')
            CALL CreateMatElement(Rows, Cols, Cnts, VvarId, IvarId)
            CALL CreateMatElement(Rows, Cols, Cnts, VvarId, VvarId)
          CASE('foil winding')
            DO j=0, Cvar % pdofs
              ! V = V0 + V1*alpha + V2*alpha^2 + ...
              CALL CreateMatElement(Rows, Cols, Cnts, VvarId, VvarId + AddIndex(j))
              IF (j/=0) THEN
                ! Circuit eqns for the pdofs:
                ! I(Vi) - I = 0
                ! ------------------------------------
                CALL CreateMatElement(Rows, Cols, Cnts, VvarId + AddIndex(j), IvarId)
                DO jj = 1, Cvar % pdofs
                    CALL CreateMatElement(Rows, Cols, Cnts, VvarId + AddIndex(j), VvarId + AddIndex(j))
                END DO
              END IF
            END DO
          CASE('flat wire')
            CALL CreateMatElement(Rows, Cols, Cnts, VvarId, VvarId)
            DO j=1, Cvar % pdofs
              CALL CreateMatElement(Rows, Cols, Cnts, VvarId, VvarId + AddIndex(j))
              CALL CreateMatElement(Rows, Cols, Cnts, VvarId + AddIndex(j), IvarId)
              CALL CreateMatElement(Rows, Cols, Cnts, VvarId + AddIndex(j), VvarId + AddIndex(j))
            END DO
          CASE('foil sheet')
            CALL CreateMatElement(Rows, Cols, Cnts, VvarId, VvarId)
            DO j=1, Comp % nCells
              CALL CreateMatElement(Rows, Cols, Cnts, VvarId, VvarId + AddIndex(j))
              CALL CreateMatElement(Rows, Cols, Cnts, VvarId + AddIndex(j), IvarId)
              DO ll=(j-1)*Comp % nSublayers + 1, j*Comp % nSublayers
                DO jj=1, Comp % nSegments
                  i = FoilSheetStrandDof(Comp % nCells, Comp % nSegments, ll, jj)
                  ! cell current balance and the strand equation
                  CALL CreateMatElement(Rows, Cols, Cnts, VvarId + AddIndex(j), VvarId + AddIndex(i))
                  CALL CreateMatElement(Rows, Cols, Cnts, VvarId + AddIndex(i), VvarId + AddIndex(j))
                  CALL CreateMatElement(Rows, Cols, Cnts, VvarId + AddIndex(i), VvarId + AddIndex(i))
                END DO
              END DO
            END DO
          END SELECT
        END IF

!        temp = SUM(Cnts)
!print *, "Active elements ", ParEnv % Mype, ":", GetNOFActive()
        DO q=GetNOFActive(),1,-1
          Element => GetActiveElement(q)
          CALL CreateComponentElements(Element, Comp, VvarId, IvarId, Rows, Cols, Cnts, Done, dofsdone)
        END DO

        DO q=GetNOFBoundaryElements(),1,-1
          Element => GetBoundaryElement(q)
          CALL CreateComponentElements(Element, Comp, VvarId, IvarId, Rows, Cols, Cnts, Done, dofsdone)
        END DO
!        Comp % nofcnts = SUM(Cnts) - temp
!        print *, ParEnv % Mype, "CompInd:", CompInd, "Coil Type:", Comp % CoilType, &
!                 "Comp % BodyId:", Comp % BodyId, "Comp % nofcnts", Comp % nofcnts

      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CreateComponentEquations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CountComponentElements(Element, Comp, RowId, Rows, Cnts, Done, dofsdone)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    TYPE(Component_t), POINTER :: Comp
    INTEGER :: nn, nd, RowId
    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: Rows(:), Cnts(:)
    LOGICAL*1 :: Done(:)
    LOGICAL :: dofsdone

    IF (ElAssocToComp(Element, Comp)) THEN
      Asolver => CurrentModel % Asolver
      nn = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element,ASolver)
      SELECT CASE (Comp % CoilType)
      CASE('stranded')           
        CALL CountAndCreateStranded(Element,nn,nd,RowId,Cnts,Done,Rows)
      CASE('massive')
        IF (HasSupport(Element,nn)) THEN
          CALL CountAndCreateMassive(Element,nn,nd,RowId,Cnts,Done,Rows)
        END IF
     CASE('foil winding')
        IF (HasSupport(Element,nn)) THEN
          CALL CountAndCreateFoilWinding(Element,nn,nd,Comp,Cnts,Done,Rows)
        END IF
     CASE('flat wire')
        IF (HasSupport(Element,nn)) THEN
          CALL CountAndCreateFlatWire(Element,nn,nd,Comp,Cnts,Done,Rows)
        END IF
     CASE('foil sheet')
        IF (HasSupport(Element,nn)) THEN
          CALL CountAndCreateFoilSheet(Element,nn,nd,Comp,Cnts,Done,Rows)
        END IF
      END SELECT
    END IF
!------------------------------------------------------------------------------
   END SUBROUTINE CountComponentElements
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CreateComponentElements(Element, Comp, VvarId, IvarId, Rows, Cols, Cnts, Done, dofsdone)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    TYPE(Component_t), POINTER :: Comp
    TYPE(Solver_t), POINTER :: ASolver
    INTEGER :: nn, nd, VvarId, IvarId
    INTEGER, POINTER :: Rows(:), Cols(:), Cnts(:)
    LOGICAL*1 :: Done(:)
    LOGICAL :: dofsdone
    
    IF (ElAssocToComp(Element, Comp)) THEN
      Asolver => CurrentModel % Asolver
      nn = GetElementNOFNodes(Element)
      nd = GetElementNOFDOFs(Element,ASolver)
      SELECT CASE (Comp % CoilType)
      CASE('stranded')
        CALL CountAndCreateStranded(Element,nn,nd,VvarId,Cnts,Done,Rows,Cols,IvarId)
      CASE('massive')
        IF (HasSupport(Element,nn)) THEN
          CALL CountAndCreateMassive(Element,nn,nd,VvarId,Cnts,Done,Rows,Cols=Cols)
        END IF
     CASE('foil winding')
        IF (HasSupport(Element,nn)) THEN
          CALL CountAndCreateFoilWinding(Element,nn,nd,Comp,Cnts,Done,Rows,Cols=Cols)
        END IF
     CASE('flat wire')
        IF (HasSupport(Element,nn)) THEN
          CALL CountAndCreateFlatWire(Element,nn,nd,Comp,Cnts,Done,Rows,Cols=Cols)
        END IF
     CASE('foil sheet')
        IF (HasSupport(Element,nn)) THEN
          CALL CountAndCreateFoilSheet(Element,nn,nd,Comp,Cnts,Done,Rows,Cols=Cols)
        END IF
      END SELECT
    END IF
!------------------------------------------------------------------------------
   END SUBROUTINE CreateComponentElements
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CountAndCreateStranded(Element,nn,nd,i,Cnts,Done,Rows,Cols,Jsind,Harmonic)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    INTEGER :: nn, nd, ncdofs1, ncdofs2, dim
    OPTIONAL :: Cols
    INTEGER :: Rows(:), Cols(:), Cnts(:)
    INTEGER :: p,i,j,k,Indexes(nd)
    INTEGER, OPTIONAL :: Jsind
    INTEGER, POINTER :: PS(:)
    LOGICAL*1 :: Done(:)
    LOGICAL :: First=.TRUE.
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm
    SAVE dim, First

    IF (First) THEN
      First = .FALSE.
      dim = CoordinateSystemDimension()
    END IF

    IF (.NOT. PRESENT(Harmonic)) THEN
      harm = CurrentModel % HarmonicCircuits
    ELSE
      harm = Harmonic
    END IF

    
    IF (.NOT. ASSOCIATED(CurrentModel % ASolver) ) CALL Fatal ('CountAndCreateStranded','ASolver not found!')
    PS => CurrentModel % Asolver % Variable % Perm

    nd = GetElementDOFs(Indexes,Element,CurrentModel % ASolver)
    IF(dim==2) THEN
      ncdofs1=1
      ncdofs2=nd
    ELSE IF(dim==3) THEN
      ncdofs1=nn+1
      ncdofs2=nd
    END IF

    DO p=ncdofs1,ncdofs2
      j = Indexes(p)

      IF( ASSOCIATED( CurrentModel % Mesh % PeriodicPerm ) ) THEN
        ! If we have periodicity eliminated only flag the master in Done
        k = CurrentModel % Mesh % PeriodicPerm(j)
        IF( k > 0 ) j = k
      END IF

      IF(.NOT.Done(j)) THEN
        Done(j) = .TRUE.
        j = PS(j)
        IF (harm) j = ReIndex(j)
        IF(PRESENT(Cols)) THEN
          CALL CreateMatElement(Rows, Cols, Cnts, i, j, harm) 
          CALL CreateMatElement(Rows, Cols, Cnts, j, Jsind, harm)
!         CALL CreateMatElement(Rows, Cols, Cnts, j, Jsind)
        ELSE
          CALL CountMatElement(Rows, Cnts, i, 1, harm)
          CALL CountMatElement(Rows, Cnts, j, 1, harm)
!         CALL CountMatElement(Rows, Cnts, j, 1)
        END IF
      END IF
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountAndCreateStranded
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE CountAndCreateMassive(Element,nn,nd,i,Cnts,Done,Rows,Cols,Harmonic)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    INTEGER :: nn, nd, ncdofs1, ncdofs2, dim
    OPTIONAL :: Cols
    INTEGER :: Rows(:), Cols(:), Cnts(:)
    INTEGER :: p,i,j,k,Indexes(nd)
    INTEGER, POINTER :: PS(:)
    LOGICAL*1 :: Done(:)
    LOGICAL :: First=.TRUE.
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm
    SAVE dim, First

    IF (First) THEN
      First = .FALSE.
      dim = CoordinateSystemDimension()
    END IF
    
    IF (.NOT. PRESENT(Harmonic)) THEN
      harm = CurrentModel % HarmonicCircuits
    ELSE
      harm = Harmonic
    END IF
    
    IF (.NOT. ASSOCIATED(CurrentModel % ASolver) ) CALL Fatal ('CountAndCreateMassive','ASolver not found!')
    PS => CurrentModel % Asolver % Variable % Perm
    nd = GetElementDOFs(Indexes,Element,CurrentModel % ASolver)
    IF(dim==2) THEN
      ncdofs1=1
      ncdofs2=nd
    ELSE IF(dim==3) THEN
      ncdofs1=nn+1
      ncdofs2=nd
    END IF
    DO p=ncdofs1,ncdofs2
      j = Indexes(p)

      IF( ASSOCIATED( CurrentModel % Mesh % PeriodicPerm ) ) THEN
        ! If we have periodicity eliminated only flag the master in Done
        k = CurrentModel % Mesh % PeriodicPerm(j)
        IF( k > 0 ) j = k
      END IF

      IF(.NOT.Done(j)) THEN
        Done(j) = .TRUE.
        j = PS(j)
        IF (harm) j = ReIndex(j)
        IF(PRESENT(Cols)) THEN
          CALL CreateMatElement(Rows, Cols, Cnts, i, j, harm)
          CALL CreateMatElement(Rows, Cols, Cnts, j, i, harm)
        ELSE
          CALL CountMatElement(Rows, Cnts, i, 1, harm)
          CALL CountMatElement(Rows, Cnts, j, 1, harm)
        END IF
      END IF
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountAndCreateMassive
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
    SUBROUTINE CountAndCreateFoilWinding(Element,nn,nd,Comp,Cnts,Done,Rows,Cols,Harmonic)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Element_t) :: Element
    TYPE(Component_t), POINTER :: Comp
    INTEGER :: nn, nd, ncdofs, dim
    OPTIONAL :: Cols
    INTEGER :: Rows(:), Cols(:), Cnts(:)
    INTEGER :: Indexes(nd)
    INTEGER :: p,j,q,vpolord,vpolordtest,vpolord_tot,&
      dofId,dofIdtest,vvarId, nm, ni, qn
    LOGICAL :: dofsdone, First=.TRUE.
    INTEGER, POINTER :: PS(:)
    LOGICAL*1 :: Done(:)
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm, ConstraintActive, Found
    SAVE dim, First

    IF (First) THEN
      First = .FALSE.
      dim = CoordinateSystemDimension()
    END IF
    
    IF (.NOT. PRESENT(Harmonic)) THEN
      harm = CurrentModel % HarmonicCircuits
    ELSE
      harm = Harmonic
    END IF
    
    IF (.NOT. ASSOCIATED(CurrentModel % ASolver) ) CALL Fatal ('CountAndCreateFoilWinding','ASolver not found!')
    PS => CurrentModel % Asolver % Variable % Perm
    nd = GetElementDOFs(Indexes,Element,CurrentModel % ASolver)
    nm = CurrentModel % ASolver % Matrix % NumberOfRows

    ncdofs=nd
    IF (dim == 3) ncdofs=nd-nn

    vvarId = Comp % vvar % ValueId
    vpolord_tot = Comp % vvar % pdofs - 1

    ! DEV-1491: reserve room for the nodal (scalar potential) couplings added by
    ! Add_foil_winding when the component activates the constraint equation.
    ! Why: without them the circuit source lives in the edge rows only and the
    ! DC (omega=0) system is inconsistent -> the linear solver floors. Before:
    ! only edge dofs were counted/created here. Limitation: gated exactly like
    ! the assembly (same keyword, 3D only), otherwise the counting pass and the
    ! creation pass disagree and Circuits_MatrixInit aborts.
    ConstraintActive = .FALSE.
    IF (dim == 3) ConstraintActive = ListGetLogical( &
        CurrentModel % Components(Comp % ComponentId) % Values, 'Activate Constraint', Found)

    DO vpolordtest=0,vpolord_tot ! V'(alpha)
      dofIdtest = AddIndex(vpolordtest + 1) + vvarId
      DO vpolord = 0, vpolord_tot ! V(alpha)
        dofId = AddIndex(vpolord + 1) + vvarId
        IF (PRESENT(Cols)) THEN  
          CALL CreateMatElement(Rows, Cols, Cnts, dofIdtest+nm, dofId+nm, harm)
        ELSE
          CALL CountMatElement(Rows, Cnts, dofIdtest+nm, 1, harm)
        END IF
      END DO

      DO j=1,ncdofs
        q=j                        
        IF (dim == 3) q=q+nn
        IF (PRESENT(Cols)) THEN  
          q = PS(Indexes(q))
          IF (harm) q = ReIndex(q)
          CALL CreateMatElement(Rows, Cols, Cnts, dofIdtest+nm, q, harm)
        ELSE
          CALL CountMatElement(Rows, Cnts, dofIdtest+nm, 1, harm)
        END IF
      END DO

      IF (ConstraintActive) THEN
        DO ni=1,nn
          IF (PRESENT(Cols)) THEN
            qn = PS(Indexes(ni))
            IF (harm) qn = ReIndex(qn)
            CALL CreateMatElement(Rows, Cols, Cnts, dofIdtest+nm, qn, harm)
          ELSE
            CALL CountMatElement(Rows, Cnts, dofIdtest+nm, 1, harm)
          END IF
        END DO
      END IF
    END DO

    DO vpolord = 0, vpolord_tot ! V(alpha)
      dofId = AddIndex(vpolord + 1) + vvarId
      DO j=1,ncdofs
        q=j
        IF (dim == 3) q=q+nn
        q = PS(Indexes(q))
        IF (harm) q = ReIndex(q)
        IF (PRESENT(Cols)) THEN  
          CALL CreateMatElement(Rows, Cols, Cnts, q, dofId+nm, harm)
        ELSE
          CALL CountMatElement(Rows, Cnts, q, 1, harm)
        END IF
      END DO

      IF (ConstraintActive) THEN
        DO ni=1,nn
          qn = PS(Indexes(ni))
          IF (harm) qn = ReIndex(qn)
          IF (PRESENT(Cols)) THEN
            CALL CreateMatElement(Rows, Cols, Cnts, qn, dofId+nm, harm)
          ELSE
            CALL CountMatElement(Rows, Cnts, qn, 1, harm)
          END IF
        END DO
      END IF
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountAndCreateFoilWinding
!------------------------------------------------------------------------------





!------------------------------------------------------------------------------
!> Matrix structure of the flat wire couplings: every cell voltage row that the
!> element can touch couples with the edge dofs of the element and vice versa.
!------------------------------------------------------------------------------
  SUBROUTINE CountAndCreateFlatWire(Element,nn,nd,Comp,Cnts,Done,Rows,Cols,Harmonic)
!------------------------------------------------------------------------------
    USE CircuitUtils
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    TYPE(Component_t), POINTER :: Comp
    INTEGER :: nn, nd
    OPTIONAL :: Cols
    INTEGER :: Rows(:), Cols(:), Cnts(:)
    INTEGER :: Indexes(nd)
    INTEGER :: j, q, k, ks, ka, ks1, ks2, ka1, ka2, dofId, vvarId, nm, ncdofs, ni, qn
    INTEGER, POINTER :: PS(:)
    LOGICAL*1 :: Done(:)
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm, ConstraintActive, Found
    REAL(KIND=dp) :: sStack(nn), sAcross(nn)

    IF (.NOT. PRESENT(Harmonic)) THEN
      harm = CurrentModel % HarmonicCircuits
    ELSE
      harm = Harmonic
    END IF

    IF (.NOT. ASSOCIATED(CurrentModel % ASolver) ) CALL Fatal ('CountAndCreateFlatWire','ASolver not found!')
    IF (CoordinateSystemDimension() /= 3) CALL Fatal('CountAndCreateFlatWire','Flat wire is implemented only in 3D!')
    PS => CurrentModel % Asolver % Variable % Perm
    nd = GetElementDOFs(Indexes,Element,CurrentModel % ASolver)
    nm = CurrentModel % ASolver % Matrix % NumberOfRows
    ncdofs = nd - nn
    vvarId = Comp % vvar % ValueId

    ! DEV-1491: reserve room for the nodal (scalar potential) couplings added by
    ! Add_flat_wire when the component activates the constraint equation. See the
    ! comment in CountAndCreateFoilWinding; the gate must match the assembly
    ! exactly or Circuits_MatrixInit aborts on an element-count mismatch.
    ConstraintActive = ListGetLogical( &
        CurrentModel % Components(Comp % ComponentId) % Values, 'Activate Constraint', Found)

    CALL GetFlatWireLocalFields(Comp % StackAlongAlpha, Element, nn, sStack, sAcross)

    ! Cells the element can touch: the nodal range plus one cell of margin
    ! for higher order elements and cell edges cutting through elements.
    ks1 = MAX(1, FLOOR(MINVAL(sStack) * Comp % nStack))
    ks2 = MIN(Comp % nStack, FLOOR(MAXVAL(sStack) * Comp % nStack) + 2)
    ka1 = MAX(1, FLOOR(MINVAL(sAcross) * Comp % nAcross))
    ka2 = MIN(Comp % nAcross, FLOOR(MAXVAL(sAcross) * Comp % nAcross) + 2)

    DO ks = ks1, ks2
      DO ka = ka1, ka2
        k = (ks-1) * Comp % nAcross + ka
        dofId = AddIndex(k, harm) + vvarId
        DO j=1,ncdofs
          q = PS(Indexes(j+nn))
          IF (harm) q = ReIndex(q)
          IF (PRESENT(Cols)) THEN
            CALL CreateMatElement(Rows, Cols, Cnts, dofId+nm, q, harm)
            CALL CreateMatElement(Rows, Cols, Cnts, q, dofId+nm, harm)
          ELSE
            CALL CountMatElement(Rows, Cnts, dofId+nm, 1, harm)
            CALL CountMatElement(Rows, Cnts, q, 1, harm)
          END IF
        END DO

        IF (ConstraintActive) THEN
          DO ni=1,nn
            qn = PS(Indexes(ni))
            IF (harm) qn = ReIndex(qn)
            IF (PRESENT(Cols)) THEN
              CALL CreateMatElement(Rows, Cols, Cnts, dofId+nm, qn, harm)
              CALL CreateMatElement(Rows, Cols, Cnts, qn, dofId+nm, harm)
            ELSE
              CALL CountMatElement(Rows, Cnts, dofId+nm, 1, harm)
              CALL CountMatElement(Rows, Cnts, qn, 1, harm)
            END IF
          END DO
        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountAndCreateFlatWire
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Matrix structure of the foil sheet couplings: every strand current row c_kj
!> the element can touch couples with the edge dofs of the element and vice
!> versa. The circuit-internal couplings (c_kj with V_k and with itself) are
!> created once per component in CreateComponentEquations.
!------------------------------------------------------------------------------
  SUBROUTINE CountAndCreateFoilSheet(Element,nn,nd,Comp,Cnts,Done,Rows,Cols,Harmonic)
!------------------------------------------------------------------------------
    USE CircuitUtils
    IMPLICIT NONE
    TYPE(Element_t), POINTER :: Element
    TYPE(Component_t), POINTER :: Comp
    INTEGER :: nn, nd
    OPTIONAL :: Cols
    INTEGER :: Rows(:), Cols(:), Cnts(:)
    INTEGER :: Indexes(nd)
    INTEGER :: j, q, k, ks, ka, ks1, ks2, ka1, ka2, dofId, vvarId, nm, ncdofs, nLayers
    INTEGER, POINTER :: PS(:)
    LOGICAL*1 :: Done(:)
    LOGICAL, OPTIONAL :: Harmonic
    LOGICAL :: harm
    REAL(KIND=dp) :: sStack(nn), sAcross(nn)

    IF (.NOT. PRESENT(Harmonic)) THEN
      harm = CurrentModel % HarmonicCircuits
    ELSE
      harm = Harmonic
    END IF

    IF (.NOT. ASSOCIATED(CurrentModel % ASolver) ) CALL Fatal ('CountAndCreateFoilSheet','ASolver not found!')
    IF (CoordinateSystemDimension() /= 3) CALL Fatal('CountAndCreateFoilSheet','Foil sheet is implemented only in 3D!')
    PS => CurrentModel % Asolver % Variable % Perm
    nd = GetElementDOFs(Indexes,Element,CurrentModel % ASolver)
    nm = CurrentModel % ASolver % Matrix % NumberOfRows
    ncdofs = nd - nn
    vvarId = Comp % vvar % ValueId

    CALL GetFlatWireLocalFields(Comp % StackAlongAlpha, Element, nn, sStack, sAcross)

    ! Strands the element can touch: the nodal range plus one cell of margin
    ! for higher order elements and strand borders cutting through elements.
    nLayers = Comp % nCells * Comp % nSublayers
    ks1 = MAX(1, FLOOR(MINVAL(sStack) * nLayers))
    ks2 = MIN(nLayers, FLOOR(MAXVAL(sStack) * nLayers) + 2)
    ka1 = MAX(1, FLOOR(MINVAL(sAcross) * Comp % nSegments))
    ka2 = MIN(Comp % nSegments, FLOOR(MAXVAL(sAcross) * Comp % nSegments) + 2)

    DO ks = ks1, ks2
      DO ka = ka1, ka2
        k = FoilSheetStrandDof(Comp % nCells, Comp % nSegments, ks, ka)
        dofId = AddIndex(k, harm) + vvarId
        DO j=1,ncdofs
          q = PS(Indexes(j+nn))
          IF (harm) q = ReIndex(q)
          IF (PRESENT(Cols)) THEN
            CALL CreateMatElement(Rows, Cols, Cnts, dofId+nm, q, harm)
            CALL CreateMatElement(Rows, Cols, Cnts, q, dofId+nm, harm)
          ELSE
            CALL CountMatElement(Rows, Cnts, dofId+nm, 1, harm)
            CALL CountMatElement(Rows, Cnts, q, 1, harm)
          END IF
        END DO
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountAndCreateFoilSheet
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE Circuits_MatrixInit()
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Matrix_t), POINTER :: CM
    TYPE(Solver_t), POINTER :: ASolver
    INTEGER, POINTER :: PS(:), Cnts(:)
    INTEGER, POINTER CONTIG :: Rows(:), Cols(:)
    INTEGER :: nm, Circuit_tot_n, n, i
    LOGICAL :: dofsdone
    LOGICAL*1, ALLOCATABLE :: Done(:)
    REAL(KIND=dp), POINTER CONTIG :: Values(:)
    LOGICAL :: Parallel, Found
    
    ASolver => CurrentModel % Asolver
    IF (.NOT.ASSOCIATED(ASolver)) CALL Fatal('Circuits_MatrixInit','ASolver not found!')
    Circuit_tot_n = CurrentModel % Circuit_tot_n
    
    ! Initialize Circuit matrix:
    ! -----------------------------
    PS => Asolver % Variable % Perm
    nm =  Asolver % Matrix % NumberOfRows

    CM => AllocateMatrix()
    CurrentModel % CircuitMatrix=>CM
    
    CM % Format = MATRIX_CRS
    Asolver % Matrix % AddMatrix => CM
    ALLOCATE(CM % RHS(nm + Circuit_tot_n)); CM % RHS=0._dp

    CM % NumberOfRows = nm + Circuit_tot_n

    n = CM % NumberOfRows

    ALLOCATE(Rows(n+1), Cnts(n)); Rows=0; Cnts=0
    ALLOCATE(Done(SIZE(PS)), CM % RowOwner(n)); Cm % RowOwner=-1


    Parallel = CurrentModel % Solver % Parallel      
    IF( Parallel ) CALL SetCircuitsParallelInfo()

    ! COUNT SIZES:
    ! ============
    dofsdone = .FALSE.
    
    CALL CountBasicCircuitEquations(Rows, Cnts)
    CALL CountComponentEquations(Rows, Cnts, Done, dofsdone)

    ! ALLOCATE CRS STRUCTURES (if need be):
    ! =====================================

    n = SUM(Cnts)

    IF (n<=0) THEN
      CM % NUmberOfRows = 0
      DEALLOCATE(Rows,Cnts,Done,CM); CM=>Null()
      Asolver %  Matrix % AddMatrix => CM
      CurrentModel % CircuitMatrix=>CM
      RETURN 
    END IF

    ALLOCATE(Cols(n+1), Values(n+1))
    Cols = 0; Values = 0._dp

    ! CREATE ROW POINTERS:
    ! ====================

    CM % NumberOfRows = nm + Circuit_tot_n
    Rows(1) = 1
    DO i=2,CM % NumberOfRows+1
      Rows(i) = Rows(i-1) + Cnts(i-1)
    END DO

    Cnts = 0

    ! CREATE COLUMNS:
    ! ===============

    CALL CreateBasicCircuitEquations(Rows, Cols, Cnts)
    CALL CreateComponentEquations(Rows, Cols, Cnts, Done, dofsdone)
    
    IF (n /= SUM(Cnts)) THEN
      CALL Fatal('Circuits_MatrixInit', &
                 'Inconsistent number of matrix elements: '//I2S(n)//' vs. '//I2S(SUM(CNTs)))
    END IF

    DEALLOCATE( Cnts, Done )
    CM % Rows => Rows
    CM % Cols => Cols
    CM % Values => Values
    CALL CRS_SortMatrix(CM)
    
    Asolver %  Matrix % AddMatrix => CM
!------------------------------------------------------------------------------
  END SUBROUTINE Circuits_MatrixInit
!------------------------------------------------------------------------------

      
END MODULE CircMatInitMod


