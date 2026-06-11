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
! *           50-062 Wroclaw, Poland
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
    Nu = (576.0_dp / El**2 + 2.873_dp / SQRT(El))**(-0.5_dp)
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
  USE htc_convection_constants
  USE htc_fluid_properties, ONLY: compute_fluid_properties
  USE htc_correlations,     ONLY: nu_churchill_chu, nu_horiz_up, nu_horiz_down, &
                                  nu_vertical_channel, nu_parallel, nu_churchill_bernstein
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: htc_evaluate_core

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
      ! g_eff = g * cos(theta) folded into Ra; slant Lc_nat baked at SIF.
      Ra     = Ra_base * Lc_nat**3 * COS(theta_deg * DEG_TO_RAD)
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

END MODULE htc_core


! Elmer shim -- external (non-module) function so it is callable by name from a
! SIF Procedure line. Reads BC parameters via DefUtils, clamps T_surface,
! asserts unit sanity, then defers to the pure htc_evaluate_core.
FUNCTION htc_evaluate(Model, n, T_surface) RESULT(h_final)
  USE DefUtils
  USE htc_core,                 ONLY: htc_evaluate_core
  USE htc_convection_constants, ONLY: FLUID_AIR
  IMPLICIT NONE
  ! dp comes from DefUtils (Elmer's Types module). It is numerically equal to
  ! SELECTED_REAL_KIND(15, 307), the kind used by htc_core / htc_convection_constants.

  TYPE(Model_t)    :: Model
  INTEGER          :: n
  REAL(KIND=dp)    :: T_surface, h_final

  TYPE(ValueList_t), POINTER :: BC, Constants
  TYPE(Variable_t),  POINTER :: htc_var
  INTEGER       :: NatOrient_in, FcdOrient_in, fluid_id
  REAL(KIND=dp) :: Lc_nat, Lc_forced, theta_deg, T_ambient_C, velocity_mag
  REAL(KIND=dp) :: channel_gap, channel_length, T_surf_K, T_amb_K
  REAL(KIND=dp), PARAMETER :: CELSIUS_TO_KELVIN = 273.15_dp
  LOGICAL       :: Found

  BC => GetBC()
  IF (.NOT. ASSOCIATED(BC)) CALL Fatal('htc_evaluate', 'No active boundary condition')
  Constants => GetConstants()

  NatOrient_in = GetInteger(BC, 'NatOrient', Found)
  IF (.NOT. Found) CALL Fatal('htc_evaluate', 'Missing SIF parameter: NatOrient')
  FcdOrient_in = GetInteger(BC, 'FcdOrient', Found)
  IF (.NOT. Found) CALL Fatal('htc_evaluate', 'Missing SIF parameter: FcdOrient')

  Lc_nat = GetConstReal(BC, 'Lc_nat', Found)
  IF (.NOT. Found) CALL Fatal('htc_evaluate', 'Missing SIF parameter: Lc_nat')

  ! T_ambient, velocity_mag and Fluid Id are uniform across the case in this
  ! workflow. Read from the BC first (per-face override); fall back to the SIF
  ! Constants block. Temperatures are Celsius (matches Elmer's Temperature).
  T_ambient_C = GetConstReal(BC, 'T_ambient', Found)
  IF (.NOT. Found .AND. ASSOCIATED(Constants)) &
    T_ambient_C = GetConstReal(Constants, 'T_ambient', Found)
  IF (.NOT. Found) CALL Fatal('htc_evaluate', 'Missing T_ambient (BC or Constants)')

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

  ! SIF Temperature variable and T_ambient are both Celsius; convert to Kelvin
  ! before the core (air-property correlations require absolute temperature).
  ! Clamp T_surface here, not in core -- Elmer may call with T_surface = 0 C on
  ! the first iteration of a fresh run, which is in range; defend against a
  ! stray sub-physical value instead.
  T_amb_K  = T_ambient_C + CELSIUS_TO_KELVIN
  T_surf_K = MAX(T_surface + CELSIUS_TO_KELVIN, 200.0_dp)
  IF (T_surf_K < 100.0_dp .OR. T_surf_K > 1500.0_dp) &
    CALL Fatal('htc_evaluate', 'T_surface out of plausible range -- check SIF units')

  h_final = htc_evaluate_core(fluid_id, T_surf_K, T_amb_K, NatOrient_in, FcdOrient_in, &
                              Lc_nat, Lc_forced, theta_deg, velocity_mag, &
                              channel_gap, channel_length)

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
END FUNCTION htc_evaluate
