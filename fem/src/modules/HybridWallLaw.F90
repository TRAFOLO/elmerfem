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
! *  Authors: TRAFOLO Engineering
! *  Email:   info@trafolo.eu
! *  Web:     https://trafolo.eu
! *  Address: Lastādijas Street 12 k-3,
! *           Riga, LV-1050, Latvia
! *
! *  Original Date: November 22, 2026
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Wall heat transfer of a conjugate heat transfer case whose air mesh is too
!> coarse to resolve the wall layer: the solid and the air exchange heat through
!> a wall formula evaluated against a reference air temperature near the wall.
!>
!> The solid side gets a heat transfer coefficient against the reference air,
!> the air side receives the same heat as a flux:
!>
!>    solid: Heat Transfer Coefficient = Variable Tsolid, Tref, Coordinate 2
!>      Real Procedure "HybridWallLaw" "WallHeatTransferCoefficient"
!>    solid: External Temperature = Equals Tref
!>    air: Heat Flux = Variable Tsolid, Tref, Coordinate 2
!>      Real Procedure "HybridWallLaw" "WallHeatFlux"
!>
!> The temperature-independent part of each formula is set up per face by the
!> case writer and given in the Boundary Condition, each keyword starting with
!> "Hybrid Wall Law":
!>    (no suffix, Integer)  the formula, one of the WALL_LAW_* codes below
!>    A, B                  h = A + B dT^(1/4); for the run laws B multiplies
!>                          (dT/x)^(1/4), x the height above Run Start
!>    Channel               e1 of the channel factor r = 1/sqrt(1 + 200.5/(e1 dT)^1.5),
!>                          0 on an open face, negative for a closed slot (r = 0)
!>    Forced                forced-flow coefficient hf, blended as (hn^3 + hf^3)^(1/3)
!>    Outside               coefficient added after the blend
!>    Run Start             height of the foot of the wall [m]
!>    Strip dT, Strip Exponent
!>                          the B term of a strip face is scaled by
!>                          clip(dT / Strip dT, 0.5, 2)^Strip Exponent; Strip dT is 50 K if not given
!> dT = Tsolid - Tref. The coefficient takes dT >= 0.02 K on both sides; the
!> flux is that coefficient times dT with its sign: where a wall is colder than
!> the air, the air gives up the heat the wall takes.
!> For HeatSolve. HeatSolveVec hands over values interpolated to integration points,
!> where a node without Tref is no longer recognisable; do not use it with these laws.
!> \ingroup UDFs
!------------------------------------------------------------------------------
MODULE HybridWallLawCore
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: WallCoefficient, WALL_LAW_NATURAL, WALL_LAW_NATURAL_STRIP, WALL_LAW_RUN, &
      WALL_LAW_RUN_STRIP, WALL_LAW_CHANNEL, WALL_LAW_MIXED, WALL_LAW_MIXED_OUTSIDE, &
      WallLawParameters_t, DEFAULT_STRIP_DT_K, MissingInput

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)

  INTEGER, PARAMETER :: WALL_LAW_NATURAL       = 1
  INTEGER, PARAMETER :: WALL_LAW_NATURAL_STRIP = 2
  INTEGER, PARAMETER :: WALL_LAW_RUN           = 3
  INTEGER, PARAMETER :: WALL_LAW_RUN_STRIP     = 4
  INTEGER, PARAMETER :: WALL_LAW_CHANNEL       = 5
  INTEGER, PARAMETER :: WALL_LAW_MIXED         = 6
  INTEGER, PARAMETER :: WALL_LAW_MIXED_OUTSIDE = 7

  REAL(KIND=dp), PARAMETER :: MIN_DIFFERENCE_K = 0.02_dp
  REAL(KIND=dp), PARAMETER :: MIN_RUN_M = 0.003_dp
  REAL(KIND=dp), PARAMETER :: DEFAULT_STRIP_DT_K = 50.0_dp

  TYPE WallLawParameters_t
    INTEGER :: Kind = 0
    REAL(KIND=dp) :: A = 0.0_dp, B = 0.0_dp, Channel = 0.0_dp, Forced = 0.0_dp, Outside = 0.0_dp
    REAL(KIND=dp) :: RunStart = 0.0_dp, StripDT = DEFAULT_STRIP_DT_K, StripExponent = 0.0_dp
  END TYPE WallLawParameters_t

CONTAINS

  !> Elmer marks a dependency it could not find at a node with HUGE; no heat passes there.
  PURE LOGICAL FUNCTION MissingInput(x)
    REAL(KIND=dp), INTENT(IN) :: x(3)
    MissingInput = ANY(ABS(x) >= 0.5_dp * HUGE(1.0_dp))
  END FUNCTION MissingInput

  PURE FUNCTION StripFactor(dt, p) RESULT(s)
    REAL(KIND=dp), INTENT(IN) :: dt
    TYPE(WallLawParameters_t), INTENT(IN) :: p
    REAL(KIND=dp) :: s
    s = MIN(MAX(dt / p % StripDT, 0.5_dp), 2.0_dp) ** p % StripExponent
  END FUNCTION StripFactor

  PURE FUNCTION ChannelFactor(dt, e1) RESULT(r)
    REAL(KIND=dp), INTENT(IN) :: dt, e1
    REAL(KIND=dp) :: r
    r = 1.0_dp
    IF (e1 > 0.0_dp) r = 1.0_dp / SQRT(1.0_dp + 200.5_dp / (e1 * dt) ** 1.5_dp)
    IF (e1 < 0.0_dp) r = 0.0_dp
  END FUNCTION ChannelFactor

  !> Coefficient of the wall formula at temperature difference dt and height y.
  !> Ok is false for an unknown formula.
  PURE SUBROUTINE WallCoefficient(p, dt, y, h, Ok)
    TYPE(WallLawParameters_t), INTENT(IN) :: p
    REAL(KIND=dp), INTENT(IN) :: dt, y
    REAL(KIND=dp), INTENT(OUT) :: h
    LOGICAL, INTENT(OUT) :: Ok
    REAL(KIND=dp) :: x, hn

    Ok = .TRUE.
    SELECT CASE (p % Kind)
    CASE (WALL_LAW_NATURAL)
      h = p % A + p % B * dt ** 0.25_dp
    CASE (WALL_LAW_NATURAL_STRIP)
      h = p % A + p % B * StripFactor(dt, p) * dt ** 0.25_dp
    CASE (WALL_LAW_RUN, WALL_LAW_RUN_STRIP)
      x = MAX(y - p % RunStart, MIN_RUN_M)
      h = p % B * dt ** 0.25_dp / x ** 0.25_dp
      IF (p % Kind == WALL_LAW_RUN_STRIP) h = h * StripFactor(dt, p)
      h = p % A + h
    CASE (WALL_LAW_CHANNEL)
      h = (p % A + p % B * dt ** 0.25_dp) * ChannelFactor(dt, p % Channel)
    CASE (WALL_LAW_MIXED, WALL_LAW_MIXED_OUTSIDE)
      ! natural and forced convection combined; the coefficient part never sees dT below the floor
      hn = (p % A + p % B * MAX(dt, MIN_DIFFERENCE_K) ** 0.25_dp) * &
          ChannelFactor(MAX(dt, MIN_DIFFERENCE_K), p % Channel)
      h = (hn ** 3 + p % Forced ** 3) ** (1.0_dp / 3.0_dp)
      IF (p % Kind == WALL_LAW_MIXED_OUTSIDE) h = h + p % Outside
    CASE DEFAULT
      h = 0.0_dp
      Ok = .FALSE.
    END SELECT
  END SUBROUTINE WallCoefficient

END MODULE HybridWallLawCore


!> Reads the formula of the current boundary element, see HybridWallLawCore.
SUBROUTINE HybridWallLawParameters(p)
  USE DefUtils
  USE HybridWallLawCore
  IMPLICIT NONE
  TYPE(WallLawParameters_t), INTENT(OUT) :: p
  TYPE(ValueList_t), POINTER :: BC
  LOGICAL :: Found

  BC => GetBC()
  IF (.NOT. ASSOCIATED(BC)) CALL Fatal('HybridWallLaw', 'Called outside a boundary condition.')
  p % Kind = GetInteger(BC, 'Hybrid Wall Law', Found)
  IF (.NOT. Found) CALL Fatal('HybridWallLaw', 'Boundary condition '//I2S(GetBCId())//' has no "Hybrid Wall Law".')
  p % A = GetConstReal(BC, 'Hybrid Wall Law A', Found)
  p % B = GetConstReal(BC, 'Hybrid Wall Law B', Found)
  p % Channel = GetConstReal(BC, 'Hybrid Wall Law Channel', Found)
  p % Forced = GetConstReal(BC, 'Hybrid Wall Law Forced', Found)
  p % Outside = GetConstReal(BC, 'Hybrid Wall Law Outside', Found)
  p % RunStart = GetConstReal(BC, 'Hybrid Wall Law Run Start', Found)
  p % StripDT = GetConstReal(BC, 'Hybrid Wall Law Strip dT', Found)
  IF (.NOT. Found) p % StripDT = DEFAULT_STRIP_DT_K
  IF (p % StripDT <= 0.0_dp .AND. (p % Kind == WALL_LAW_NATURAL_STRIP .OR. p % Kind == WALL_LAW_RUN_STRIP)) &
      CALL Fatal('HybridWallLaw', 'Boundary condition '//I2S(GetBCId())// &
      ': "Hybrid Wall Law Strip dT" must be positive.')
  p % StripExponent = GetConstReal(BC, 'Hybrid Wall Law Strip Exponent', Found)
END SUBROUTINE HybridWallLawParameters


!> Heat transfer coefficient of the solid against the reference air.
!> x = (Tsolid, Tref, height).
FUNCTION WallHeatTransferCoefficient(Model, n, x) RESULT(h)
  USE DefUtils
  USE HybridWallLawCore
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: x(*), h
  TYPE(WallLawParameters_t) :: p
  LOGICAL :: Ok

  h = 0.0_dp
  IF (MissingInput(x(1:3))) RETURN
  CALL HybridWallLawParameters(p)
  CALL WallCoefficient(p, MAX(x(1) - x(2), 0.02_dp), x(3), h, Ok)
  IF (.NOT. Ok) CALL Fatal('HybridWallLaw', 'Unknown "Hybrid Wall Law" '//I2S(p % Kind)//'.')
END FUNCTION WallHeatTransferCoefficient


!> Heat flux from the wall into the air, the same heat the solid loses.
!> x = (Tsolid, Tref, height).
FUNCTION WallHeatFlux(Model, n, x) RESULT(q)
  USE DefUtils
  USE HybridWallLawCore
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: x(*), q
  TYPE(WallLawParameters_t) :: p
  REAL(KIND=dp) :: h
  LOGICAL :: Ok

  q = 0.0_dp
  IF (MissingInput(x(1:3))) RETURN
  CALL HybridWallLawParameters(p)
  CALL WallCoefficient(p, MAX(x(1) - x(2), 0.02_dp), x(3), h, Ok)
  IF (.NOT. Ok) CALL Fatal('HybridWallLaw', 'Unknown "Hybrid Wall Law" '//I2S(p % Kind)//'.')
  q = h * (x(1) - x(2))
END FUNCTION WallHeatFlux


!------------------------------------------------------------------------------
!> Heat conductivity of the air at the walls of such a case in forced flow:
!>
!>    Heat Conductivity = Variable Wall Distance
!>      Real Procedure "HybridWallLaw" "FilmConductivity"
!>
!> The wall heat enters the air as a flux that does not depend on the air at the
!> wall. Air that cannot carry it away would get hotter than the wall. Two things
!> stand for the transport the coarse air solve does not resolve:
!>  - a well-mixed film in air that stands nearly still within the film thickness
!>    of a wall (pockets, dead corners): the film conductivity, fading linearly to
!>    that of the air as the fastest node of the element reaches the film speed;
!>  - the cells at a wall conduct at least like the wall formula: wall coefficient
!>    times the height of the cell, which keeps the air at the wall near the wall
!>    temperature also where the cell is thicker than the thermal layer.
!> Material keywords, each starting with "Hybrid Film":
!>    Air Conductivity   conductivity of the air away from the walls
!>    Conductivity       conductivity of the film
!>    Thickness          wall distance up to which the film is full; it is gone at 1.3 times it
!>    Speed              air speed at which the film is gone; not given: the film does not fade
!>    Wall Coefficient   heat transfer coefficient the cells at a wall conduct at least [W/m2K]
!> x = (wall distance). Needs the field Flow Solution. For HeatSolve, which asks
!> for the conductivity node by node of the element it assembles.
!> \ingroup UDFs
!------------------------------------------------------------------------------
FUNCTION FilmConductivity(Model, n, x) RESULT(k)
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: x(*), k

  REAL(KIND=dp), PARAMETER :: FADE_RATIO = 1.3_dp, AT_WALL_M = 1.0e-6_dp
  TYPE(Element_t), POINTER :: Element, Parent
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Mesh_t), POINTER, SAVE :: KnownMesh => NULL()
  TYPE(Variable_t), POINTER, SAVE :: Distance => NULL(), Flow => NULL()
  INTEGER, SAVE :: KnownBody = -1, KnownElement = -1
  REAL(KIND=dp), SAVE :: AirK, FilmK, Thickness, StillSpeed, WallCoefficient, Still, WallCellK
  REAL(KIND=dp) :: d, NearestM, FarthestM, Fastest
  INTEGER :: i, j, node, dofs
  LOGICAL :: Found, FromFace

  Element => Model % CurrentElement
  IF (.NOT. ASSOCIATED(Element)) CALL Fatal('HybridWallLaw', 'The film was asked for outside an element.')
  ! Asked for on a face (a boundary condition that takes the conductivity of its parent):
  ! the film is that of the air element behind the face.
  FromFace = ASSOCIATED(Element % BoundaryInfo)
  IF (FromFace) THEN
    Parent => Element % BoundaryInfo % Left
    IF (.NOT. HasFilm(Parent)) Parent => Element % BoundaryInfo % Right
    IF (.NOT. HasFilm(Parent)) CALL Fatal('HybridWallLaw', 'Face without an air element that has a film.')
    Element => Parent
  END IF
  Mesh => GetMesh()
  IF (.NOT. ASSOCIATED(Mesh, KnownMesh)) THEN
    Distance => VariableGet(Mesh % Variables, 'Wall Distance')
    Flow => VariableGet(Mesh % Variables, 'Flow Solution')
    IF (.NOT. ASSOCIATED(Distance) .OR. .NOT. ASSOCIATED(Flow)) &
        CALL Fatal('HybridWallLaw', 'The film needs the fields Wall Distance and Flow Solution.')
    IF (.NOT. ASSOCIATED(Distance % Perm) .OR. .NOT. ASSOCIATED(Flow % Perm) .OR. &
        Distance % TYPE /= Variable_on_nodes .OR. Flow % TYPE /= Variable_on_nodes .OR. &
        Distance % DOFs /= 1 .OR. Flow % DOFs /= CoordinateSystemDimension() + 1) &
        CALL Fatal('HybridWallLaw', 'The film needs the nodal fields Wall Distance and Flow Solution (velocity, pressure).')
    KnownMesh => Mesh
    KnownBody = -1
  END IF
  IF (Element % BodyId /= KnownBody) THEN
    Material => GetMaterial(Element)
    AirK = GetConstReal(Material, 'Hybrid Film Air Conductivity', Found)
    IF (.NOT. Found .OR. AirK <= 0.0_dp) &
        CALL Fatal('HybridWallLaw', '"Hybrid Film Air Conductivity" must be given and positive.')
    FilmK = GetConstReal(Material, 'Hybrid Film Conductivity', Found)
    IF (.NOT. Found .OR. FilmK <= 0.0_dp) &
        CALL Fatal('HybridWallLaw', '"Hybrid Film Conductivity" must be given and positive.')
    Thickness = GetConstReal(Material, 'Hybrid Film Thickness', Found)
    IF (.NOT. Found .OR. Thickness <= 0.0_dp) &
        CALL Fatal('HybridWallLaw', '"Hybrid Film Thickness" must be given and positive.')
    StillSpeed = GetConstReal(Material, 'Hybrid Film Speed', Found)
    IF (.NOT. Found) StillSpeed = 0.0_dp
    WallCoefficient = GetConstReal(Material, 'Hybrid Film Wall Coefficient', Found)
    IF (.NOT. Found) WallCoefficient = 0.0_dp
    KnownBody = Element % BodyId
    KnownElement = -1
  END IF

  ! Once per element: HeatSolve starts every element with its first node. A face asks for
  ! its own nodes only, so its parent is taken anew each time.
  IF (FromFace .OR. Element % ElementIndex /= KnownElement .OR. n == Element % NodeIndexes(1)) THEN
    NearestM = HUGE(NearestM)
    FarthestM = 0.0_dp
    Fastest = 0.0_dp
    dofs = Flow % DOFs
    DO i = 1, Element % TYPE % NumberOfNodes
      node = Element % NodeIndexes(i)
      j = Distance % Perm(node)
      IF (j <= 0 .OR. Flow % Perm(node) <= 0) &
          CALL Fatal('HybridWallLaw', 'The film needs Wall Distance and Flow Solution at every node of its elements.')
      NearestM = MIN(NearestM, Distance % Values(j))
      FarthestM = MAX(FarthestM, Distance % Values(j))
      j = Flow % Perm(node)
      Fastest = MAX(Fastest, SQRT(SUM(Flow % Values(dofs*(j-1)+1:dofs*j-1)**2)))
    END DO
    Still = 1.0_dp
    IF (StillSpeed > 0.0_dp) Still = MIN(MAX(1.0_dp - Fastest / StillSpeed, 0.0_dp), 1.0_dp)
    WallCellK = 0.0_dp
    IF (NearestM <= AT_WALL_M) WallCellK = WallCoefficient * FarthestM
    KnownElement = Element % ElementIndex
  END IF

  d = MAX(x(1), 0.0_dp)
  k = AirK + (FilmK - AirK) * Still * &
      MIN(MAX((FADE_RATIO * Thickness - d) / ((FADE_RATIO - 1.0_dp) * Thickness), 0.0_dp), 1.0_dp)
  k = MAX(k, WallCellK)

CONTAINS

  FUNCTION HasFilm(Candidate) RESULT(Has)
    TYPE(Element_t), POINTER :: Candidate
    LOGICAL :: Has
    TYPE(ValueList_t), POINTER :: Values

    Has = .FALSE.
    IF (.NOT. ASSOCIATED(Candidate)) RETURN
    IF (Candidate % BodyId <= 0) RETURN
    Values => GetMaterial(Candidate)
    IF (ASSOCIATED(Values)) Has = ListCheckPresent(Values, 'Hybrid Film Air Conductivity')
  END FUNCTION HasFilm
END FUNCTION FilmConductivity
