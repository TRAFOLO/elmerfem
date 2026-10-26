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
! *  Original Date: June 9, 2026
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
!> Temperature-dependent convective heat-transfer coefficient h(T_surface) for a
!> boundary face. Wired into the HeatSolver, evaluated at every nonlinear
!> iteration (steady state) or timestep (transient), via:
!>
!>    Heat Transfer Coefficient = Variable Temperature
!>      Real Procedure "htc_udf" "htc_evaluate"
!>
!> The temperature-INDEPENDENT per-face data (orientation codes + characteristic
!> lengths) is classified offline by TRAFOLO and written into each Boundary
!> Condition block as keywords; this routine does the temperature-DEPENDENT
!> physics live: film temperature -> fluid properties -> Ra/Re -> Nu -> h_nat,
!> h_forced -> Churchill-Usagi blend -> h_final.
!>
!> Temperatures are degrees Celsius at the SIF boundary (matching Elmer's
!> Temperature variable and External Temperature); converted to Kelvin
!> internally for the property correlations. SI elsewhere (m, m/s, W/m^2/K).
!>
!> Per-BC keywords:
!>    NatOrient (Integer)  1=HU 2=HD 3=V 4=VC 5=INC
!>    FcdOrient (Integer)  0=none 1=WW 2=LW 3=PAR
!>    Lc_nat, Lc_forced, theta_deg, channel_gap, channel_length (Real, m / deg)
!>    T_ambient, velocity_mag (Real) -- per-BC, or once in the Constants block
!>    Fluid Id (Integer)   1=air (default)
!> \ingroup UDFs
!------------------------------------------------------------------------------

! Physical constants, real-kind, and integer regime codes. No Elmer dependency.
MODULE htc_convection_constants
  IMPLICIT NONE
  PUBLIC

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)

  REAL(KIND=dp), PARAMETER :: STANDARD_GRAVITY              = 9.80665_dp
  REAL(KIND=dp), PARAMETER :: R_AIR                         = 287.058_dp
  REAL(KIND=dp), PARAMETER :: P_ATM                         = 101325.0_dp
  REAL(KIND=dp), PARAMETER :: STILL_AIR_VELOCITY_THRESHOLD  = 0.1_dp
  REAL(KIND=dp), PARAMETER :: DEG_TO_RAD                    = 0.017453292519943295_dp

  ! Natural-convection orientation codes (SIF: NatOrient)
  INTEGER, PARAMETER :: HORIZONTAL_UP    = 1
  INTEGER, PARAMETER :: HORIZONTAL_DOWN  = 2
  INTEGER, PARAMETER :: VERTICAL         = 3
  INTEGER, PARAMETER :: VERTICAL_CHANNEL = 4
  INTEGER, PARAMETER :: INCLINED         = 5

  ! Forced-convection orientation codes (SIF: FcdOrient)
  INTEGER, PARAMETER :: FCD_NONE     = 0
  INTEGER, PARAMETER :: WINDWARD     = 1
  INTEGER, PARAMETER :: LEEWARD      = 2
  INTEGER, PARAMETER :: PARALLEL_FCD = 3

  ! Working-fluid codes (SIF: Fluid Id). Only air is implemented today;
  ! oils / water / glycol plug into compute_fluid_properties below.
  INTEGER, PARAMETER :: FLUID_AIR = 1

  ! Range of the air property fits; temperatures outside it are taken at the nearest limit.
  REAL(KIND=dp), PARAMETER :: T_MIN_K = 200.0_dp
  REAL(KIND=dp), PARAMETER :: T_MAX_K = 1500.0_dp
  REAL(KIND=dp), PARAMETER :: T_AMBIENT_FALLBACK_K = 293.15_dp
  ! Lengths and speeds are kept in a range where Ra, Re and h stay finite.
  REAL(KIND=dp), PARAMETER :: LENGTH_MIN = 1.0e-6_dp
  REAL(KIND=dp), PARAMETER :: LENGTH_MAX = 1.0e3_dp
  REAL(KIND=dp), PARAMETER :: LENGTH_FALLBACK = 0.05_dp
  REAL(KIND=dp), PARAMETER :: VELOCITY_MAX = 1.0e3_dp
  ! Returned when the correlations still give no finite coefficient (W/m^2/K).
  REAL(KIND=dp), PARAMETER :: H_FALLBACK = 10.0_dp
END MODULE htc_convection_constants


! Working-fluid thermophysical properties at film temperature. Input is Kelvin
! (conversion from Celsius happens at the SIF boundary). Air uses Sutherland
! viscosity + Tsilingiris polynomial fits for k, cp.
MODULE htc_fluid_properties
  USE htc_convection_constants, ONLY: dp, P_ATM, R_AIR, FLUID_AIR
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: compute_fluid_properties

CONTAINS

  ! Dispatcher on the SIF "Fluid Id". Air is the only fluid implemented; to add
  ! a coolant (transformer oil, water, glycol) give it a code in
  ! htc_convection_constants and a branch here returning the same eight fields.
  ! NOTE: liquids need temperature-dependent fits -- a single datasheet value is
  ! a 20-50% error because mu and Pr swing 2-5x over a normal operating range.
  PURE SUBROUTINE compute_fluid_properties(fluid_id, T_film_K, rho, mu, k_fluid, cp, beta, Pr, nu)
    INTEGER,       INTENT(IN)  :: fluid_id
    REAL(KIND=dp), INTENT(IN)  :: T_film_K
    REAL(KIND=dp), INTENT(OUT) :: rho, mu, k_fluid, cp, beta, Pr, nu

    SELECT CASE (fluid_id)
    ! CASE (FLUID_OIL);   CALL compute_oil_properties(T_film_K, rho, mu, k_fluid, cp, beta, Pr, nu)
    CASE DEFAULT   ! FLUID_AIR
      CALL compute_air_properties(T_film_K, rho, mu, k_fluid, cp, beta, Pr, nu)
    END SELECT
  END SUBROUTINE compute_fluid_properties

  PURE SUBROUTINE compute_air_properties(T_film_K, rho, mu, k_air, cp, beta, Pr, nu)
    REAL(KIND=dp), INTENT(IN)  :: T_film_K
    REAL(KIND=dp), INTENT(OUT) :: rho, mu, k_air, cp, beta, Pr, nu

    rho   = P_ATM / (R_AIR * T_film_K)
    mu    = 1.458e-6_dp * T_film_K**1.5_dp / (T_film_K + 110.4_dp)
    k_air = 0.02624_dp * (T_film_K / 300.0_dp)**0.8646_dp
    cp    = 1030.5_dp - 0.19975_dp * T_film_K + 3.9734e-4_dp * T_film_K**2
    beta  = 1.0_dp / T_film_K
    Pr    = mu * cp / k_air
    nu    = mu / rho
  END SUBROUTINE compute_air_properties

END MODULE htc_fluid_properties


! Nusselt-number correlations for natural and forced convection.
! Every function applies the Nu >= 1 conduction floor before returning.
MODULE htc_correlations
  USE htc_convection_constants, ONLY: dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: nu_churchill_chu, nu_horiz_up, nu_horiz_down, &
            nu_vertical_channel, nu_parallel, nu_churchill_bernstein

CONTAINS

  ! Churchill-Chu (1975) for a vertical isothermal plate.
  PURE FUNCTION nu_churchill_chu(Ra, Pr) RESULT(Nu)
    REAL(KIND=dp), INTENT(IN) :: Ra, Pr
    REAL(KIND=dp) :: Nu, denom, psi
    IF (Ra <= 1.0e9_dp) THEN
      denom = (1.0_dp + (0.492_dp / Pr)**(9.0_dp / 16.0_dp))**(4.0_dp / 9.0_dp)
      Nu = 0.68_dp + 0.67_dp * Ra**0.25_dp / denom
    ELSE
      psi = (1.0_dp + (0.492_dp / Pr)**(9.0_dp / 16.0_dp))**(-8.0_dp / 27.0_dp)
      Nu = (0.825_dp + 0.387_dp * Ra**(1.0_dp / 6.0_dp) * psi)**2
    END IF
    Nu = MAX(Nu, 1.0_dp)
  END FUNCTION nu_churchill_chu

  ! McAdams / Lloyd-Moran for hot-face-up (or cold-face-down) horizontal plate.
  PURE FUNCTION nu_horiz_up(Ra) RESULT(Nu)
    REAL(KIND=dp), INTENT(IN) :: Ra
    REAL(KIND=dp) :: Nu
    IF (Ra <= 1.0e7_dp) THEN
      Nu = 0.54_dp * Ra**0.25_dp
    ELSE
      Nu = 0.15_dp * Ra**(1.0_dp / 3.0_dp)
    END IF
    Nu = MAX(Nu, 1.0_dp)
  END FUNCTION nu_horiz_up

  ! Radziemska-Lewandowski / Incropera for hot-face-down horizontal plate.
  PURE FUNCTION nu_horiz_down(Ra) RESULT(Nu)
    REAL(KIND=dp), INTENT(IN) :: Ra
    REAL(KIND=dp) :: Nu
    Nu = 0.27_dp * Ra**0.25_dp
    Nu = MAX(Nu, 1.0_dp)
  END FUNCTION nu_horiz_down

  ! Bar-Cohen-Rohsenow Elenbaas for a vertical parallel-plate channel.
  ! Ra_s is the Rayleigh based on gap s; L is the streamwise channel length.
  PURE FUNCTION nu_vertical_channel(Ra_s, s, L) RESULT(Nu)
    REAL(KIND=dp), INTENT(IN) :: Ra_s, s, L
    REAL(KIND=dp) :: Nu, El
    El = Ra_s * (s / L)
    ! El = 0 (surface at ambient temperature) is the conduction limit; skip the division by zero.
    Nu = 1.0_dp
    IF (El > 0.0_dp) Nu = (576.0_dp / El**2 + 2.873_dp / SQRT(El))**(-0.5_dp)
    Nu = MAX(Nu, 1.0_dp)
  END FUNCTION nu_vertical_channel

  ! Pohlhausen (laminar) / Whitaker (turbulent) flat-plate parallel flow.
  PURE FUNCTION nu_parallel(Re, Pr) RESULT(Nu)
    REAL(KIND=dp), INTENT(IN) :: Re, Pr
    REAL(KIND=dp) :: Nu, Pr13
    Pr13 = Pr**(1.0_dp / 3.0_dp)
    IF (Re < 5.0e5_dp) THEN
      Nu = 0.664_dp * Re**0.5_dp * Pr13
    ELSE
      Nu = 0.037_dp * Re**0.8_dp * Pr13
    END IF
    Nu = MAX(Nu, 1.0_dp)
  END FUNCTION nu_parallel

  ! Churchill-Bernstein (1977) for bluff-body / cylinder cross-flow.
  ! Used for windward faces; leeward applies a 0.33 multiplier in the core.
  PURE FUNCTION nu_churchill_bernstein(Re, Pr) RESULT(Nu)
    REAL(KIND=dp), INTENT(IN) :: Re, Pr
    REAL(KIND=dp) :: Nu, denom1, term1, term2, Pr13
    Pr13   = Pr**(1.0_dp / 3.0_dp)
    denom1 = (1.0_dp + (0.4_dp / Pr)**(2.0_dp / 3.0_dp))**(1.0_dp / 4.0_dp)
    term1  = 0.3_dp + 0.62_dp * Re**0.5_dp * Pr13 / denom1
    term2  = (1.0_dp + (Re / 282000.0_dp)**(5.0_dp / 8.0_dp))**(4.0_dp / 5.0_dp)
    Nu = term1 * term2
    Nu = MAX(Nu, 1.0_dp)
  END FUNCTION nu_churchill_bernstein

END MODULE htc_correlations


! Pure inner function -- all physics, no Elmer types -- so it stays mathematically
! total (Nu=1 conduction floor instead of aborting on degenerate Ra/Re).
MODULE htc_core
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE, IEEE_IS_NAN
  USE htc_convection_constants
  USE htc_fluid_properties, ONLY: compute_fluid_properties
  USE htc_correlations,     ONLY: nu_churchill_chu, nu_horiz_up, nu_horiz_down, &
                                  nu_vertical_channel, nu_parallel, nu_churchill_bernstein
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: htc_evaluate_core, htc_evaluate_safe

  ! Bits of the status of htc_evaluate_safe, one per kind of input it had to correct.
  INTEGER, PARAMETER, PUBLIC :: STATUS_HOT_SURFACE  = 0
  INTEGER, PARAMETER, PUBLIC :: STATUS_COLD_SURFACE = 1
  INTEGER, PARAMETER, PUBLIC :: STATUS_NAN_SURFACE  = 2
  INTEGER, PARAMETER, PUBLIC :: STATUS_BAD_AMBIENT  = 3
  INTEGER, PARAMETER, PUBLIC :: STATUS_BAD_LENGTH   = 4
  INTEGER, PARAMETER, PUBLIC :: STATUS_BAD_VELOCITY = 5
  INTEGER, PARAMETER, PUBLIC :: STATUS_BAD_ANGLE    = 6
  INTEGER, PARAMETER, PUBLIC :: STATUS_BAD_RESULT   = 7
  INTEGER, PARAMETER, PUBLIC :: STATUS_BITS         = 8

CONTAINS

  PURE FUNCTION htc_evaluate_core(fluid_id, T_surface_K, T_ambient_K, NatOrient_in, FcdOrient_in, &
                                  Lc_nat, Lc_forced, theta_deg, velocity_mag, &
                                  channel_gap, channel_length) RESULT(h_final)
    INTEGER,       INTENT(IN) :: fluid_id
    REAL(KIND=dp), INTENT(IN) :: T_surface_K, T_ambient_K
    INTEGER,       INTENT(IN) :: NatOrient_in, FcdOrient_in
    REAL(KIND=dp), INTENT(IN) :: Lc_nat, Lc_forced, theta_deg, velocity_mag
    REAL(KIND=dp), INTENT(IN) :: channel_gap, channel_length
    REAL(KIND=dp) :: h_final

    REAL(KIND=dp) :: T_film_K, rho, mu, k_fluid, cp, beta, Pr, nu_fluid, alpha
    REAL(KIND=dp) :: delta_T, Ra_base, Ra, Ra_s, Nu_nat, Nu_fcd
    REAL(KIND=dp) :: h_nat, h_forced, Re, h_div_length

    T_film_K = 0.5_dp * (T_surface_K + T_ambient_K)
    CALL compute_fluid_properties(fluid_id, T_film_K, rho, mu, k_fluid, cp, beta, Pr, nu_fluid)
    alpha   = k_fluid / (rho * cp)
    delta_T = ABS(T_surface_K - T_ambient_K)
    Ra_base = STANDARD_GRAVITY * beta * delta_T / (nu_fluid * alpha)

    ! Natural-convection branch.
    Nu_nat       = 1.0_dp
    h_div_length = Lc_nat
    SELECT CASE (NatOrient_in)
    CASE (HORIZONTAL_UP)
      Ra     = Ra_base * Lc_nat**3
      Nu_nat = nu_horiz_up(Ra)
    CASE (HORIZONTAL_DOWN)
      Ra     = Ra_base * Lc_nat**3
      Nu_nat = nu_horiz_down(Ra)
    CASE (VERTICAL)
      Ra     = Ra_base * Lc_nat**3
      Nu_nat = nu_churchill_chu(Ra, Pr)
    CASE (INCLINED)
      ! g_eff = g * |cos(theta)| folded into Ra; slant Lc_nat baked at SIF.
      Ra     = Ra_base * Lc_nat**3 * ABS(COS(theta_deg * DEG_TO_RAD))
      Nu_nat = nu_churchill_chu(Ra, Pr)
    CASE (VERTICAL_CHANNEL)
      Ra_s         = Ra_base * channel_gap**3
      Nu_nat       = nu_vertical_channel(Ra_s, channel_gap, channel_length)
      h_div_length = channel_gap   ! Elenbaas Lc is the gap s, not Lc_nat.
    CASE DEFAULT
      Nu_nat = 1.0_dp
    END SELECT
    h_nat = Nu_nat * k_fluid / h_div_length

    ! Forced-convection branch (skip if still air or no forced orient).
    IF (FcdOrient_in == FCD_NONE .OR. velocity_mag < STILL_AIR_VELOCITY_THRESHOLD) THEN
      h_final = h_nat
      RETURN
    END IF

    Re = velocity_mag * Lc_forced / nu_fluid
    SELECT CASE (FcdOrient_in)
    CASE (WINDWARD)
      Nu_fcd = nu_churchill_bernstein(Re, Pr)
    CASE (LEEWARD)
      Nu_fcd = MAX(0.33_dp * nu_churchill_bernstein(Re, Pr), 1.0_dp)
    CASE (PARALLEL_FCD)
      Nu_fcd = nu_parallel(Re, Pr)
    CASE DEFAULT
      Nu_fcd = 1.0_dp
    END SELECT
    h_forced = Nu_fcd * k_fluid / Lc_forced

    ! Churchill-Usagi hybrid blend, n=3 (assisting / mixed convection).
    h_final = (h_nat**3 + h_forced**3)**(1.0_dp / 3.0_dp)
  END FUNCTION htc_evaluate_core

  ! htc_evaluate_core for any input, physical or not: a finite h >= 0 always. Inputs outside the range of the
  ! model are taken at the nearest valid value, and each kind of correction sets its STATUS_* bit.
  PURE SUBROUTINE htc_evaluate_safe(fluid_id, T_surface_K, T_ambient_K, NatOrient_in, FcdOrient_in, &
                                    Lc_nat, Lc_forced, theta_deg, velocity_mag, &
                                    channel_gap, channel_length, h_final, status)
    INTEGER,       INTENT(IN)  :: fluid_id, NatOrient_in, FcdOrient_in
    REAL(KIND=dp), INTENT(IN)  :: T_surface_K, T_ambient_K, Lc_nat, Lc_forced, theta_deg, velocity_mag
    REAL(KIND=dp), INTENT(IN)  :: channel_gap, channel_length
    REAL(KIND=dp), INTENT(OUT) :: h_final
    INTEGER,       INTENT(OUT) :: status

    REAL(KIND=dp) :: T_s, T_a, L_nat, L_fcd, theta, speed, gap, ch_length
    INTEGER :: nat_orient

    status = 0
    T_a = T_ambient_K
    IF (IEEE_IS_NAN(T_a)) THEN
      T_a = T_AMBIENT_FALLBACK_K
      status = IBSET(status, STATUS_BAD_AMBIENT)
    ELSE IF (T_a < T_MIN_K .OR. T_a > T_MAX_K) THEN
      T_a = MIN(MAX(T_a, T_MIN_K), T_MAX_K)
      status = IBSET(status, STATUS_BAD_AMBIENT)
    END IF

    T_s = T_surface_K
    IF (IEEE_IS_NAN(T_s)) THEN
      T_s = T_a
      status = IBSET(status, STATUS_NAN_SURFACE)
    ELSE IF (T_s > T_MAX_K) THEN
      T_s = T_MAX_K
      status = IBSET(status, STATUS_HOT_SURFACE)
    ELSE IF (T_s < T_MIN_K) THEN
      T_s = T_MIN_K
      status = IBSET(status, STATUS_COLD_SURFACE)
    END IF

    L_nat = Lc_nat
    CALL fix_length(L_nat, LENGTH_FALLBACK, status)

    nat_orient = NatOrient_in
    gap = channel_gap
    ch_length = channel_length
    IF (nat_orient == VERTICAL_CHANNEL) THEN
      IF (is_length(gap) .AND. is_length(ch_length)) THEN
        CALL fix_length(gap, LENGTH_FALLBACK, status)
        CALL fix_length(ch_length, LENGTH_FALLBACK, status)
      ELSE
        ! No usable channel: the face is taken as a vertical plate.
        nat_orient = VERTICAL
        status = IBSET(status, STATUS_BAD_LENGTH)
      END IF
    END IF

    theta = theta_deg
    IF (nat_orient == INCLINED .AND. .NOT. IEEE_IS_FINITE(theta)) THEN
      theta = 0.0_dp
      status = IBSET(status, STATUS_BAD_ANGLE)
    END IF

    speed = velocity_mag
    IF (IEEE_IS_NAN(speed)) THEN
      speed = 0.0_dp
      status = IBSET(status, STATUS_BAD_VELOCITY)
    ELSE IF (speed < 0.0_dp .OR. speed > VELOCITY_MAX) THEN
      speed = MIN(ABS(speed), VELOCITY_MAX)
      status = IBSET(status, STATUS_BAD_VELOCITY)
    END IF

    ! TRAFOLO writes Lc_forced = 0 for still air: it only has to be valid where forced convection acts.
    L_fcd = Lc_forced
    IF (FcdOrient_in /= FCD_NONE .AND. speed >= STILL_AIR_VELOCITY_THRESHOLD) &
      CALL fix_length(L_fcd, L_nat, status)

    h_final = htc_evaluate_core(fluid_id, T_s, T_a, nat_orient, FcdOrient_in, L_nat, L_fcd, theta, speed, &
                                gap, ch_length)
    IF (.NOT. IEEE_IS_FINITE(h_final) .OR. h_final < 0.0_dp) THEN
      h_final = H_FALLBACK
      status = IBSET(status, STATUS_BAD_RESULT)
    END IF
  END SUBROUTINE htc_evaluate_safe

  PURE LOGICAL FUNCTION is_length(L)
    REAL(KIND=dp), INTENT(IN) :: L
    is_length = IEEE_IS_FINITE(L) .AND. L > 0.0_dp
  END FUNCTION is_length

  ! A length that is not positive and finite becomes fallback; one outside [LENGTH_MIN, LENGTH_MAX] its limit.
  PURE SUBROUTINE fix_length(L, fallback, status)
    REAL(KIND=dp), INTENT(INOUT) :: L
    REAL(KIND=dp), INTENT(IN)    :: fallback
    INTEGER,       INTENT(INOUT) :: status
    IF (.NOT. is_length(L)) THEN
      L = fallback
    ELSE IF (L < LENGTH_MIN .OR. L > LENGTH_MAX) THEN
      L = MIN(MAX(L, LENGTH_MIN), LENGTH_MAX)
    ELSE
      RETURN
    END IF
    status = IBSET(status, STATUS_BAD_LENGTH)
  END SUBROUTINE fix_length

END MODULE htc_core


! Elmer shim -- external (non-module) function so it is callable by name from a
! SIF Procedure line. Reads BC parameters via DefUtils and defers to htc_evaluate_safe,
! which corrects any unphysical input. It never stops the run: missing parameters get
! defaults, and each kind of correction is warned about once per process.
FUNCTION htc_evaluate(Model, n, T_surface) RESULT(h_final)
  USE DefUtils
  USE htc_core,                 ONLY: htc_evaluate_safe, STATUS_HOT_SURFACE, STATUS_COLD_SURFACE, &
                                      STATUS_NAN_SURFACE, STATUS_BAD_AMBIENT, STATUS_BAD_LENGTH, &
                                      STATUS_BAD_VELOCITY, STATUS_BAD_ANGLE, STATUS_BAD_RESULT, STATUS_BITS
  USE htc_convection_constants, ONLY: FLUID_AIR, VERTICAL, FCD_NONE, T_MIN_K, T_MAX_K, H_FALLBACK
  IMPLICIT NONE
  ! dp comes from DefUtils (Elmer's Types module). It is numerically equal to
  ! SELECTED_REAL_KIND(15, 307), the kind used by htc_core / htc_convection_constants.

  TYPE(Model_t)    :: Model
  INTEGER          :: n
  REAL(KIND=dp)    :: T_surface, h_final

  TYPE(ValueList_t), POINTER :: BC, Constants
  TYPE(Variable_t),  POINTER :: htc_var
  INTEGER       :: NatOrient_in, FcdOrient_in, fluid_id, status
  REAL(KIND=dp) :: Lc_nat, Lc_forced, theta_deg, T_ambient_C, velocity_mag
  REAL(KIND=dp) :: channel_gap, channel_length
  REAL(KIND=dp), PARAMETER :: CELSIUS_TO_KELVIN = 273.15_dp
  REAL(KIND=dp), PARAMETER :: T_AMBIENT_DEFAULT_C = 20.0_dp
  ! The report of a missing parameter, after those of the STATUS_* bits.
  INTEGER, PARAMETER :: MISSING_PARAMETER = STATUS_BITS
  LOGICAL       :: Found, Missing
  LOGICAL, SAVE :: Reported(0:STATUS_BITS) = .FALSE.

  BC => GetBC()
  IF (.NOT. ASSOCIATED(BC)) THEN
    CALL Report(MISSING_PARAMETER, 'Called without a boundary condition, h = '//TRIM(Num(H_FALLBACK))// &
        ' W/m^2/K is used.')
    h_final = H_FALLBACK
    RETURN
  END IF
  Constants => GetConstants()
  Missing = .FALSE.

  NatOrient_in = GetInteger(BC, 'NatOrient', Found)
  IF (.NOT. Found) THEN
    NatOrient_in = VERTICAL
    Missing = .TRUE.
  END IF
  FcdOrient_in = GetInteger(BC, 'FcdOrient', Found)
  IF (.NOT. Found) THEN
    FcdOrient_in = FCD_NONE
    Missing = .TRUE.
  END IF
  Lc_nat = GetConstReal(BC, 'Lc_nat', Found)
  IF (.NOT. Found) THEN
    Lc_nat = 0.0_dp    ! htc_evaluate_safe takes LENGTH_FALLBACK
    Missing = .TRUE.
  END IF

  ! T_ambient, velocity_mag and Fluid Id are uniform across the case in this
  ! workflow. Read from the BC first (per-face override); fall back to the SIF
  ! Constants block, then to the BC's External Temperature. Temperatures are
  ! Celsius (matches Elmer's Temperature).
  T_ambient_C = GetConstReal(BC, 'T_ambient', Found)
  IF (.NOT. Found .AND. ASSOCIATED(Constants)) &
    T_ambient_C = GetConstReal(Constants, 'T_ambient', Found)
  IF (.NOT. Found) T_ambient_C = GetConstReal(BC, 'External Temperature', Found)
  IF (.NOT. Found) THEN
    T_ambient_C = T_AMBIENT_DEFAULT_C
    Missing = .TRUE.
  END IF
  IF (Missing) CALL Report(MISSING_PARAMETER, 'A convection parameter (NatOrient, FcdOrient, Lc_nat or '// &
      'T_ambient) is missing, defaults are used: vertical face, no forced flow, 0.05 m, 20 C.')

  velocity_mag = GetConstReal(BC, 'velocity_mag', Found)
  IF (.NOT. Found .AND. ASSOCIATED(Constants)) &
    velocity_mag = GetConstReal(Constants, 'velocity_mag', Found)
  IF (.NOT. Found) velocity_mag = 0.0_dp    ! genuine still-air case

  fluid_id = GetInteger(BC, 'Fluid Id', Found)
  IF (.NOT. Found .AND. ASSOCIATED(Constants)) fluid_id = GetInteger(Constants, 'Fluid Id', Found)
  IF (.NOT. Found) fluid_id = FLUID_AIR

  Lc_forced      = GetConstReal(BC, 'Lc_forced',      Found); IF (.NOT. Found) Lc_forced      = 0.0_dp
  theta_deg      = GetConstReal(BC, 'theta_deg',      Found); IF (.NOT. Found) theta_deg      = 0.0_dp
  channel_gap    = GetConstReal(BC, 'channel_gap',    Found); IF (.NOT. Found) channel_gap    = 0.0_dp
  channel_length = GetConstReal(BC, 'channel_length', Found); IF (.NOT. Found) channel_length = 0.0_dp

  ! SIF Temperature variable and T_ambient are both Celsius; the core works in
  ! Kelvin (air-property correlations require absolute temperature).
  CALL htc_evaluate_safe(fluid_id, T_surface + CELSIUS_TO_KELVIN, T_ambient_C + CELSIUS_TO_KELVIN, &
                         NatOrient_in, FcdOrient_in, Lc_nat, Lc_forced, theta_deg, velocity_mag, &
                         channel_gap, channel_length, h_final, status)
  IF (status /= 0) CALL ReportCorrections()

  ! Expose h as a nodal field for the full-results output. The "Htc" variable
  ! is declared on the HeatSolver only when Save Full Results is enabled, so
  ! this is a no-op (VariableGet returns null) otherwise.
  htc_var => VariableGet(Model % Variables, 'Htc', ThisOnly = .TRUE.)
  IF (ASSOCIATED(htc_var)) THEN
    IF (ASSOCIATED(htc_var % Perm)) THEN
      IF (htc_var % Perm(n) > 0) htc_var % Values(htc_var % Perm(n)) = h_final
    ELSE
      htc_var % Values(n) = h_final
    END IF
  END IF

CONTAINS

  SUBROUTINE ReportCorrections()
    IF (BTEST(status, STATUS_HOT_SURFACE)) CALL Report(STATUS_HOT_SURFACE, 'Surface temperature '// &
        TRIM(Num(T_surface))//' C is above the range of the convection model, h taken at '// &
        TRIM(Num(T_MAX_K - CELSIUS_TO_KELVIN))//' C. Check the losses and the cooling.')
    IF (BTEST(status, STATUS_COLD_SURFACE)) CALL Report(STATUS_COLD_SURFACE, 'Surface temperature '// &
        TRIM(Num(T_surface))//' C is below the range of the convection model, h taken at '// &
        TRIM(Num(T_MIN_K - CELSIUS_TO_KELVIN))//' C. Check the heat sources and the ambient temperature.')
    IF (BTEST(status, STATUS_NAN_SURFACE)) CALL Report(STATUS_NAN_SURFACE, &
        'Surface temperature is not a number, h taken at the ambient temperature.')
    IF (BTEST(status, STATUS_BAD_AMBIENT)) CALL Report(STATUS_BAD_AMBIENT, 'Ambient temperature '// &
        TRIM(Num(T_ambient_C))//' C is outside the range of the convection model, the nearest valid value is used.')
    IF (BTEST(status, STATUS_BAD_LENGTH)) CALL Report(STATUS_BAD_LENGTH, 'A characteristic length (Lc_nat '// &
        TRIM(Num(Lc_nat))//', Lc_forced '//TRIM(Num(Lc_forced))//', channel '//TRIM(Num(channel_gap))//' x '// &
        TRIM(Num(channel_length))//' m) is not a positive length up to 1000 m, a valid one is used.')
    IF (BTEST(status, STATUS_BAD_VELOCITY)) CALL Report(STATUS_BAD_VELOCITY, 'Air velocity '// &
        TRIM(Num(velocity_mag))//' m/s is not a valid speed, the nearest valid value is used.')
    IF (BTEST(status, STATUS_BAD_ANGLE)) CALL Report(STATUS_BAD_ANGLE, 'Inclination '//TRIM(Num(theta_deg))// &
        ' deg is not a number, 0 deg is used.')
    IF (BTEST(status, STATUS_BAD_RESULT)) CALL Report(STATUS_BAD_RESULT, &
        'The convection model gave no finite h, '//TRIM(Num(H_FALLBACK))//' W/m^2/K is used.')
  END SUBROUTINE ReportCorrections

  SUBROUTINE Report(kind, text)
    INTEGER,          INTENT(IN) :: kind
    CHARACTER(LEN=*), INTENT(IN) :: text
    IF (Reported(kind)) RETURN
    Reported(kind) = .TRUE.
    CALL Warn('htc_evaluate', text)
  END SUBROUTINE Report

  ! Plain for ordinary magnitudes; exponent form for the others, NaN and Infinity.
  FUNCTION Num(x) RESULT(text)
    REAL(KIND=dp), INTENT(IN) :: x
    CHARACTER(LEN=32) :: text
    IF (ABS(x) >= 0.1_dp .AND. ABS(x) < 1.0e6_dp) THEN
      WRITE(text, '(F0.1)') x
    ELSE
      WRITE(text, '(ES12.4E3)') x
    END IF
    text = ADJUSTL(text)
  END FUNCTION Num
END FUNCTION htc_evaluate
