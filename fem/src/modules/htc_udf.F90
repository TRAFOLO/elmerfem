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
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PUBLIC
  ! Not re-exported: the blanket PUBLIC above would otherwise hand
  ! IEEE_IS_FINITE to every consumer and risk a clash with DefUtils in the
  ! shim. Use the clamp / finite_or wrappers instead.
  PRIVATE :: IEEE_IS_FINITE

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

  ! ---------------------------------------------------------------------------
  ! Numerical guard rails (DEV-1465)
  !
  ! WHY: this UDF runs inside Elmer's boundary assembly on every nonlinear
  ! iteration, so whatever it returns goes straight into the stiffness matrix.
  ! A non-finite, zero or negative h poisons the matrix and ElmerSolver then
  ! aborts with `Fatal: Norm of solution appears to be NaN`, leaving the user
  ! with no result at all. The first Picard iterate of a lossy case also
  ! routinely overshoots by orders of magnitude, because h is evaluated at the
  ! previous iterate where dT is still ~0 and h is therefore at its smallest.
  !
  ! BEFORE: no output guards existed. Degenerate-but-legal SIF input returned
  ! +Inf (Lc_nat = 0; channel_gap = 0 on a VERTICAL_CHANNEL face, which TRAFOLO
  ! emits on every face today; Lc_forced = 0 with forced convection), NaN
  ! (T_surface below 0 K), or a negative h (Lc_nat < 0); and the shim aborted
  ! the whole solve via Fatal once T_surface passed 1500 K.
  !
  ! LIMITATION: clamping keeps the solve alive and the matrix well posed, but a
  ! clamped h is an extrapolation, not validated physics. The correlations and
  ! the dry-air property fits are only meaningful up to ~1500 K; above that the
  ! returned h is h(1500 K), which understates the true value (h rises with T),
  ! so the reported temperature is on the conservative/hot side. The shim warns
  ! once per category so a clamped run is never mistaken for a clean one.
  ! ---------------------------------------------------------------------------

  REAL(KIND=dp), PARAMETER :: CELSIUS_TO_KELVIN = 273.15_dp

  ! Validity band of the Sutherland viscosity and Tsilingiris k/cp fits.
  REAL(KIND=dp), PARAMETER :: T_FILM_MIN_K  = 100.0_dp
  REAL(KIND=dp), PARAMETER :: T_FILM_MAX_K  = 1500.0_dp

  ! Surface temperature accepted from Elmer before clamping. The lower bound is
  ! sub-cryogenic and the upper bound is the correlation validity limit; both
  ! exist to keep the arithmetic finite, not to model anything.
  REAL(KIND=dp), PARAMETER :: T_SURF_MIN_K  = 100.0_dp
  REAL(KIND=dp), PARAMETER :: T_SURF_MAX_K  = 1500.0_dp

  ! Smallest characteristic length treated as real geometry: 1 um. Anything
  ! smaller (or non-positive, or non-finite, or above LENGTH_MAX) is a
  ! generator/meshing artefact rather than a small face, so it is replaced by
  ! LENGTH_FALLBACK instead of clamped to LENGTH_MIN. Clamping to 1 um would
  ! give h = Nu*k/L of ~3e4 W/m2K -- finite, but silently over-cooling the part
  ! by three orders of magnitude; the 10 mm fallback lands on an ordinary
  ! still-air h of a few W/m2K, i.e. errs hot rather than cold. The upper bound
  ! matters because Ra scales with Lc**3 and would overflow to Inf first.
  REAL(KIND=dp), PARAMETER :: LENGTH_MIN      = 1.0e-6_dp
  REAL(KIND=dp), PARAMETER :: LENGTH_MAX      = 1.0e3_dp
  REAL(KIND=dp), PARAMETER :: LENGTH_FALLBACK = 1.0e-2_dp

  REAL(KIND=dp), PARAMETER :: PR_MIN        = 1.0e-3_dp   ! keeps Pr**(1/3) finite
  REAL(KIND=dp), PARAMETER :: EL_MIN        = 1.0e-12_dp  ! keeps 576/El**2 finite
  REAL(KIND=dp), PARAMETER :: NU_FLOOR      = 1.0_dp      ! pure-conduction limit
  REAL(KIND=dp), PARAMETER :: THETA_MAX_DEG = 89.0_dp     ! g_eff -> 0 at 90 deg

  ! h band. Still air gives 2-25 W/m2K and forced air 10-500, so H_MAX only
  ! ever catches a numerical blow-up. H_MIN is strictly positive so the Robin
  ! boundary term stays positive definite.
  REAL(KIND=dp), PARAMETER :: H_MIN         = 1.0e-3_dp
  REAL(KIND=dp), PARAMETER :: H_MAX         = 1.0e5_dp

  ! Last-resort h when the whole evaluation came out non-finite: the order of
  ! magnitude of still-air natural convection. Always accompanied by a warning.
  REAL(KIND=dp), PARAMETER :: H_FALLBACK    = 5.0_dp

CONTAINS

  ! Clamp x into [lo, hi]. A non-finite x collapses to lo, so NaN can never
  ! survive a clamp (MAX/MIN alone propagate NaN on most compilers).
  PURE FUNCTION clamp(x, lo, hi) RESULT(y)
    REAL(KIND=dp), INTENT(IN) :: x, lo, hi
    REAL(KIND=dp) :: y
    IF (.NOT. IEEE_IS_FINITE(x)) THEN
      y = lo
    ELSE
      y = MIN(MAX(x, lo), hi)
    END IF
  END FUNCTION clamp

  ! x if finite, otherwise fallback. Used where a clamp would hide the scale.
  PURE FUNCTION finite_or(x, fallback) RESULT(y)
    REAL(KIND=dp), INTENT(IN) :: x, fallback
    REAL(KIND=dp) :: y
    IF (IEEE_IS_FINITE(x)) THEN
      y = x
    ELSE
      y = fallback
    END IF
  END FUNCTION finite_or

  ! .TRUE. when L can be used as a characteristic length (finite, 1 um..1 km).
  ! Kept as a predicate so the shim can report exactly which SIF keyword was
  ! rejected while the core silently substitutes LENGTH_FALLBACK.
  PURE FUNCTION length_is_valid(L) RESULT(ok)
    REAL(KIND=dp), INTENT(IN) :: L
    LOGICAL :: ok
    ok = IEEE_IS_FINITE(L) .AND. (L >= LENGTH_MIN) .AND. (L <= LENGTH_MAX)
  END FUNCTION length_is_valid

  ! L if usable, otherwise the 10 mm fallback. See LENGTH_FALLBACK above.
  PURE FUNCTION safe_length(L) RESULT(y)
    REAL(KIND=dp), INTENT(IN) :: L
    REAL(KIND=dp) :: y
    IF (length_is_valid(L)) THEN
      y = L
    ELSE
      y = LENGTH_FALLBACK
    END IF
  END FUNCTION safe_length

END MODULE htc_convection_constants


! Working-fluid thermophysical properties at film temperature. Input is Kelvin
! (conversion from Celsius happens at the SIF boundary). Air uses Sutherland
! viscosity + Tsilingiris polynomial fits for k, cp.
MODULE htc_fluid_properties
  USE htc_convection_constants, ONLY: dp, P_ATM, R_AIR, FLUID_AIR, &
                                      T_FILM_MIN_K, T_FILM_MAX_K, clamp
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

    REAL(KIND=dp) :: T

    ! Clamp into the validity band of the fits before any division or
    ! fractional power.
    ! WHY: every expression below is singular or complex-valued at T <= 0 K --
    ! rho and beta divide by T, mu and k_air raise T to a non-integer power.
    ! BEFORE: T_film_K was used raw, so a diverging Elmer iterate (T_surface
    ! below 0 K, or NaN) produced Inf/NaN properties that propagated into h and
    ! then into the stiffness matrix.
    ! LIMITATION: outside [100, 1500] K the returned properties are the
    ! boundary values, not an extrapolation of real air -- and real air is
    ! dissociating well before 1500 K anyway. The caller is expected to warn.
    T = clamp(T_film_K, T_FILM_MIN_K, T_FILM_MAX_K)

    rho   = P_ATM / (R_AIR * T)
    mu    = 1.458e-6_dp * T**1.5_dp / (T + 110.4_dp)
    k_air = 0.02624_dp * (T / 300.0_dp)**0.8646_dp
    cp    = 1030.5_dp - 0.19975_dp * T + 3.9734e-4_dp * T**2
    beta  = 1.0_dp / T
    Pr    = mu * cp / k_air
    nu    = mu / rho
  END SUBROUTINE compute_air_properties

END MODULE htc_fluid_properties


! Nusselt-number correlations for natural and forced convection.
! Every function applies the Nu >= 1 conduction floor before returning.
! Each function also sanitises its own arguments before the first fractional
! power and passes the result through nu_finalise.
! WHY: Ra and Re reach these functions from Elmer's live temperature field, so a
! diverging iterate can hand over a negative or non-finite value. A real
! exponent applied to a negative base yields NaN, and `MAX(Nu, 1.0)` does NOT
! remove a NaN -- with a NaN first argument gfortran may return either operand,
! so the floor silently passed NaN through under some optimisation settings.
! BEFORE: only the MAX floor existed; nu_vertical_channel additionally divided
! by El and by SQRT(El) with no guard, returning Inf/NaN for a zero gap or a
! zero channel length.
! LIMITATION: a clamped Ra/Re means the correlation is being evaluated outside
! the regime it was fitted for; the value returned is the nearest valid one, so
! it is an estimate rather than validated physics.
MODULE htc_correlations
  USE htc_convection_constants, ONLY: dp, PR_MIN, EL_MIN, NU_FLOOR, clamp, finite_or
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: nu_churchill_chu, nu_horiz_up, nu_horiz_down, &
            nu_vertical_channel, nu_parallel, nu_churchill_bernstein

CONTAINS

  ! Apply the conduction floor and drop any non-finite result. Shared by every
  ! correlation so the floor cannot be bypassed by a NaN.
  PURE FUNCTION nu_finalise(Nu_raw) RESULT(Nu)
    REAL(KIND=dp), INTENT(IN) :: Nu_raw
    REAL(KIND=dp) :: Nu
    Nu = MAX(finite_or(Nu_raw, NU_FLOOR), NU_FLOOR)
  END FUNCTION nu_finalise

  ! Churchill-Chu (1975) for a vertical isothermal plate.
  PURE FUNCTION nu_churchill_chu(Ra_in, Pr_in) RESULT(Nu)
    REAL(KIND=dp), INTENT(IN) :: Ra_in, Pr_in
    REAL(KIND=dp) :: Nu, denom, psi, Ra, Pr
    Ra = clamp(Ra_in, 0.0_dp, HUGE(1.0_dp))
    Pr = clamp(Pr_in, PR_MIN, HUGE(1.0_dp))
    IF (Ra <= 1.0e9_dp) THEN
      denom = (1.0_dp + (0.492_dp / Pr)**(9.0_dp / 16.0_dp))**(4.0_dp / 9.0_dp)
      Nu = 0.68_dp + 0.67_dp * Ra**0.25_dp / denom
    ELSE
      psi = (1.0_dp + (0.492_dp / Pr)**(9.0_dp / 16.0_dp))**(-8.0_dp / 27.0_dp)
      Nu = (0.825_dp + 0.387_dp * Ra**(1.0_dp / 6.0_dp) * psi)**2
    END IF
    Nu = nu_finalise(Nu)
  END FUNCTION nu_churchill_chu

  ! McAdams / Lloyd-Moran for hot-face-up (or cold-face-down) horizontal plate.
  PURE FUNCTION nu_horiz_up(Ra_in) RESULT(Nu)
    REAL(KIND=dp), INTENT(IN) :: Ra_in
    REAL(KIND=dp) :: Nu, Ra
    Ra = clamp(Ra_in, 0.0_dp, HUGE(1.0_dp))
    IF (Ra <= 1.0e7_dp) THEN
      Nu = 0.54_dp * Ra**0.25_dp
    ELSE
      Nu = 0.15_dp * Ra**(1.0_dp / 3.0_dp)
    END IF
    Nu = nu_finalise(Nu)
  END FUNCTION nu_horiz_up

  ! Radziemska-Lewandowski / Incropera for hot-face-down horizontal plate.
  PURE FUNCTION nu_horiz_down(Ra_in) RESULT(Nu)
    REAL(KIND=dp), INTENT(IN) :: Ra_in
    REAL(KIND=dp) :: Nu, Ra
    Ra = clamp(Ra_in, 0.0_dp, HUGE(1.0_dp))
    Nu = 0.27_dp * Ra**0.25_dp
    Nu = nu_finalise(Nu)
  END FUNCTION nu_horiz_down

  ! Bar-Cohen-Rohsenow Elenbaas for a vertical parallel-plate channel.
  ! Ra_s is the Rayleigh based on gap s; L is the streamwise channel length.
  PURE FUNCTION nu_vertical_channel(Ra_s, s, L) RESULT(Nu)
    REAL(KIND=dp), INTENT(IN) :: Ra_s, s, L
    REAL(KIND=dp) :: Nu, El
    ! El floor keeps 576/El**2 and 2.873/SQRT(El) finite. El -> EL_MIN drives
    ! Nu -> 0, which nu_finalise lifts to the conduction floor -- the physically
    ! right answer for a channel too narrow or too short to establish a plume.
    El = clamp(Ra_s * (s / L), EL_MIN, HUGE(1.0_dp))
    Nu = (576.0_dp / El**2 + 2.873_dp / SQRT(El))**(-0.5_dp)
    Nu = nu_finalise(Nu)
  END FUNCTION nu_vertical_channel

  ! Pohlhausen (laminar) / Whitaker (turbulent) flat-plate parallel flow.
  PURE FUNCTION nu_parallel(Re_in, Pr_in) RESULT(Nu)
    REAL(KIND=dp), INTENT(IN) :: Re_in, Pr_in
    REAL(KIND=dp) :: Nu, Pr13, Re, Pr
    Re = clamp(Re_in, 0.0_dp, HUGE(1.0_dp))
    Pr = clamp(Pr_in, PR_MIN, HUGE(1.0_dp))
    Pr13 = Pr**(1.0_dp / 3.0_dp)
    IF (Re < 5.0e5_dp) THEN
      Nu = 0.664_dp * Re**0.5_dp * Pr13
    ELSE
      Nu = 0.037_dp * Re**0.8_dp * Pr13
    END IF
    Nu = nu_finalise(Nu)
  END FUNCTION nu_parallel

  ! Churchill-Bernstein (1977) for bluff-body / cylinder cross-flow.
  ! Used for windward faces; leeward applies a 0.33 multiplier in the core.
  PURE FUNCTION nu_churchill_bernstein(Re_in, Pr_in) RESULT(Nu)
    REAL(KIND=dp), INTENT(IN) :: Re_in, Pr_in
    REAL(KIND=dp) :: Nu, denom1, term1, term2, Pr13, Re, Pr
    Re = clamp(Re_in, 0.0_dp, HUGE(1.0_dp))
    Pr = clamp(Pr_in, PR_MIN, HUGE(1.0_dp))
    Pr13   = Pr**(1.0_dp / 3.0_dp)
    denom1 = (1.0_dp + (0.4_dp / Pr)**(2.0_dp / 3.0_dp))**(1.0_dp / 4.0_dp)
    term1  = 0.3_dp + 0.62_dp * Re**0.5_dp * Pr13 / denom1
    term2  = (1.0_dp + (Re / 282000.0_dp)**(5.0_dp / 8.0_dp))**(4.0_dp / 5.0_dp)
    Nu = term1 * term2
    Nu = nu_finalise(Nu)
  END FUNCTION nu_churchill_bernstein

END MODULE htc_correlations


! Pure inner function -- all physics, no Elmer types -- so it stays mathematically
! total (Nu=1 conduction floor instead of aborting on degenerate Ra/Re).
!
! Total by construction: for ANY combination of arguments htc_evaluate_core
! returns a finite h in [H_MIN, H_MAX]. It never aborts and never emits NaN,
! Inf, zero or a negative value.
!
! WHY: the return value goes directly into Elmer's boundary stiffness matrix. A
! non-finite entry makes ElmerSolver abort with "Norm of solution appears to be
! NaN", and a negative h inverts the sign of the Robin term so the surface gains
! heat as it warms -- a numerical runaway on top of whatever physical one
! exists. Either way the user gets no result file at all.
!
! BEFORE: the only guard was `Nu = MAX(Nu, 1)` inside each correlation. The
! function returned +Inf for Lc_nat = 0, for VERTICAL_CHANNEL with a zero gap or
! zero channel length, and for forced convection with Lc_forced = 0; a negative
! h for Lc_nat < 0; and NaN for T_surface below 0 K. All of these are reachable
! from a legal SIF (channel_gap = Real 0.0 is emitted on every face of the
! current TRAFOLO output) or from a diverging Elmer iterate.
!
! LIMITATION: the guards make the function total, not accurate. Clamped
! temperatures, substituted lengths and floored Nusselt numbers all mean the
! correlation is being read outside its fitted regime. Callers that can report
! (the htc_evaluate shim) must warn; this module stays silent because it is
! PURE.
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
    REAL(KIND=dp) :: T_s, T_a, L_nat, L_fcd, gap, chan_L, theta, vel

    ! --- sanitise every input before it reaches an arithmetic operation ------
    T_s    = clamp(T_surface_K, T_SURF_MIN_K, T_SURF_MAX_K)
    T_a    = clamp(T_ambient_K, T_SURF_MIN_K, T_SURF_MAX_K)
    L_nat  = safe_length(Lc_nat)
    gap    = safe_length(channel_gap)
    chan_L = safe_length(channel_length)
    ! theta is folded into g_eff as COS(theta); beyond 89 deg the plate is
    ! effectively horizontal, COS goes to zero and then negative, which would
    ! feed a negative Ra into a fractional power.
    theta  = clamp(theta_deg, 0.0_dp, THETA_MAX_DEG)
    vel    = clamp(velocity_mag, 0.0_dp, HUGE(1.0_dp))

    T_film_K = 0.5_dp * (T_s + T_a)
    CALL compute_fluid_properties(fluid_id, T_film_K, rho, mu, k_fluid, cp, beta, Pr, nu_fluid)
    alpha   = k_fluid / (rho * cp)
    delta_T = ABS(T_s - T_a)
    Ra_base = STANDARD_GRAVITY * beta * delta_T / (nu_fluid * alpha)

    ! Natural-convection branch.
    Nu_nat       = NU_FLOOR
    h_div_length = L_nat
    SELECT CASE (NatOrient_in)
    CASE (HORIZONTAL_UP)
      Ra     = Ra_base * L_nat**3
      Nu_nat = nu_horiz_up(Ra)
    CASE (HORIZONTAL_DOWN)
      Ra     = Ra_base * L_nat**3
      Nu_nat = nu_horiz_down(Ra)
    CASE (VERTICAL)
      Ra     = Ra_base * L_nat**3
      Nu_nat = nu_churchill_chu(Ra, Pr)
    CASE (INCLINED)
      ! g_eff = g * cos(theta) folded into Ra; slant Lc_nat baked at SIF.
      Ra     = Ra_base * L_nat**3 * COS(theta * DEG_TO_RAD)
      Nu_nat = nu_churchill_chu(Ra, Pr)
    CASE (VERTICAL_CHANNEL)
      Ra_s         = Ra_base * gap**3
      Nu_nat       = nu_vertical_channel(Ra_s, gap, chan_L)
      h_div_length = gap   ! Elenbaas Lc is the gap s, not Lc_nat.
    CASE DEFAULT
      Nu_nat = NU_FLOOR
    END SELECT
    h_nat = Nu_nat * k_fluid / h_div_length

    ! Forced-convection branch (skip if still air or no forced orient).
    IF (FcdOrient_in == FCD_NONE .OR. vel < STILL_AIR_VELOCITY_THRESHOLD) THEN
      ! finite_or before clamp: clamp() collapses a non-finite value to its
      ! lower bound, which for h means an almost-adiabatic wall -- finite but
      ! silently wrong. With every input sanitised above, NaN should now be
      ! unreachable; if it ever is reached, H_FALLBACK is the honest answer.
      h_final = clamp(finite_or(h_nat, H_FALLBACK), H_MIN, H_MAX)
      RETURN
    END IF

    ! Lc_forced is only read on this path, so it is validated here rather than
    ! above: a still-air face legitimately carries Lc_forced = Real 0.0.
    L_fcd = safe_length(Lc_forced)
    Re = vel * L_fcd / nu_fluid
    SELECT CASE (FcdOrient_in)
    CASE (WINDWARD)
      Nu_fcd = nu_churchill_bernstein(Re, Pr)
    CASE (LEEWARD)
      Nu_fcd = MAX(0.33_dp * nu_churchill_bernstein(Re, Pr), NU_FLOOR)
    CASE (PARALLEL_FCD)
      Nu_fcd = nu_parallel(Re, Pr)
    CASE DEFAULT
      Nu_fcd = NU_FLOOR
    END SELECT
    h_forced = Nu_fcd * k_fluid / L_fcd

    ! Churchill-Usagi hybrid blend, n=3 (assisting / mixed convection).
    ! Both terms are clamped first so the cube cannot overflow to Inf.
    h_nat    = clamp(finite_or(h_nat,    H_FALLBACK), H_MIN, H_MAX)
    h_forced = clamp(finite_or(h_forced, H_FALLBACK), H_MIN, H_MAX)
    h_final  = clamp((h_nat**3 + h_forced**3)**(1.0_dp / 3.0_dp), H_MIN, H_MAX)
  END FUNCTION htc_evaluate_core

END MODULE htc_core


! Elmer shim -- external (non-module) function so it is callable by name from a
! SIF Procedure line. Reads BC parameters via DefUtils, sanitises and clamps
! them, then defers to the pure htc_evaluate_core.
!
! This function must never abort.
!
! WHY: it is called from inside Elmer's boundary assembly, once per boundary
! node per nonlinear iteration. Any Fatal here kills ElmerSolver mid-assembly,
! so no result file is written at all and the user is left with nothing --
! including for cases that would have converged perfectly well one iteration
! later. The first Picard iterate of a lossy case routinely overshoots by orders
! of magnitude, because h is evaluated at the previous iterate where dT is still
! ~0 and h is therefore at its smallest.
!
! BEFORE (DEV-1465): this shim called Fatal in six places -- no active BC, four
! missing SIF keywords, and `T_surface out of plausible range` for T above
! 1500 K. The last one was the reported failure: on
! fullbrickdeneme1_20Aug_regrouped_meshingTest the run stopped with
!   ERROR:: htc_evaluate: T_surface out of plausible range -- check SIF units
!   STOP 1
! and an empty results directory, at a point where the global temperature norm
! was only 207 -- a single hot boundary node was enough to lose the whole solve.
! Note the low-temperature half of that check was already dead code: T_surf_K
! had been raised by MAX(.., 200) on the line above, so it could never trip the
! < 100 branch.
!
! NOW: every one of those paths clamps or substitutes a documented default and
! warns once per category per MPI rank. A malformed SIF therefore degrades to a
! still-air conduction estimate and says so in the log, instead of aborting.
!
! LIMITATION: a run that logs any of these warnings is not a validated result.
! Clamping the surface temperature at 1500 K in particular means the reported
! peak temperature is a lower bound on a runaway, not a prediction -- h(1500 K)
! understates h(T_true), so the solver is being handed less cooling than the
! (already extrapolated) correlations would give. The temperature field is still
! the right thing to show the user: it makes the runaway visible.
FUNCTION htc_evaluate(Model, n, T_surface) RESULT(h_final)
  USE DefUtils
  USE htc_core,                 ONLY: htc_evaluate_core
  USE htc_convection_constants, ONLY: FLUID_AIR, T_SURF_MIN_K, T_SURF_MAX_K, &
                                      H_MIN, H_MAX, H_FALLBACK, LENGTH_FALLBACK, &
                                      THETA_MAX_DEG, CELSIUS_TO_KELVIN, FCD_NONE, &
                                      VERTICAL_CHANNEL, STILL_AIR_VELOCITY_THRESHOLD, &
                                      length_is_valid, clamp
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
  REAL(KIND=dp) :: channel_gap, channel_length, T_surf_K, T_amb_K, T_raw, h_raw
  LOGICAL       :: Found

  ! Default ambient used only when the SIF supplies none. 20 C is the value
  ! TRAFOLO emits on every convective face today.
  REAL(KIND=dp), PARAMETER :: T_AMBIENT_DEFAULT_C = 20.0_dp

  ! Warn-once slots. SAVEd locals persist for the life of the process, so each
  ! MPI rank reports each category at most once instead of once per boundary
  ! node per iteration (which would be tens of thousands of lines).
  INTEGER, PARAMETER :: W_NO_BC   =  1, W_NAT     =  2, W_FCD     =  3
  INTEGER, PARAMETER :: W_LCNAT   =  4, W_TAMB    =  5, W_LCBAD   =  6
  INTEGER, PARAMETER :: W_LCFCD   =  7, W_CHAN    =  8, W_TNAN    =  9
  INTEGER, PARAMETER :: W_THI     = 10, W_TLO     = 11, W_HCLAMP  = 12
  INTEGER, PARAMETER :: W_THETA   = 13, W_BANNER  = 14
  INTEGER, PARAMETER :: N_WARN    = 14
  LOGICAL, SAVE      :: warned(N_WARN) = .FALSE.

  ! 512 chars, not MAX_NAME_LEN: Elmer's MAX_NAME_LEN is shorter than the
  ! diagnostic lines below, and an internal WRITE that overruns its record
  ! aborts the run with "Fortran runtime error: End of record" -- the exact
  ! class of failure this file exists to prevent.
  CHARACTER(LEN=512) :: buf

  CALL say_once(W_BANNER, 'htc_udf guard rails active: T_surface clamped to ' // &
                          '[100,1500] K, h clamped to [1e-3,1e5] W/m2K, no Fatal path.')

  BC => GetBC()
  IF (.NOT. ASSOCIATED(BC)) THEN
    ! Nothing can be read, so there is no defensible h. Return the still-air
    ! order of magnitude rather than aborting; the warning is the audit trail.
    CALL say_once(W_NO_BC, 'No active boundary condition; using fallback h.')
    h_final = H_FALLBACK
    RETURN
  END IF
  Constants => GetConstants()

  NatOrient_in = GetInteger(BC, 'NatOrient', Found)
  IF (.NOT. Found) THEN
    ! 0 is not a valid orientation code, so htc_evaluate_core takes its
    ! CASE DEFAULT branch and returns the pure-conduction limit Nu = 1.
    CALL say_once(W_NAT, 'Missing SIF parameter NatOrient; using conduction limit (Nu=1).')
    NatOrient_in = 0
  END IF

  FcdOrient_in = GetInteger(BC, 'FcdOrient', Found)
  IF (.NOT. Found) THEN
    CALL say_once(W_FCD, 'Missing SIF parameter FcdOrient; assuming still air.')
    FcdOrient_in = FCD_NONE
  END IF

  Lc_nat = GetConstReal(BC, 'Lc_nat', Found)
  IF (.NOT. Found) THEN
    CALL say_once(W_LCNAT, 'Missing SIF parameter Lc_nat; using 10 mm fallback.')
    Lc_nat = LENGTH_FALLBACK
  ELSE IF (.NOT. length_is_valid(Lc_nat)) THEN
    WRITE(buf,'(A,ES12.4,A)') 'Lc_nat = ', Lc_nat, &
      ' is not a usable length; using 10 mm fallback.'
    CALL say_once(W_LCBAD, buf)
  END IF

  ! T_ambient, velocity_mag and Fluid Id are uniform across the case in this
  ! workflow. Read from the BC first (per-face override); fall back to the SIF
  ! Constants block. Temperatures are Celsius (matches Elmer's Temperature).
  T_ambient_C = GetConstReal(BC, 'T_ambient', Found)
  IF (.NOT. Found .AND. ASSOCIATED(Constants)) &
    T_ambient_C = GetConstReal(Constants, 'T_ambient', Found)
  IF (.NOT. Found) &
    T_ambient_C = GetConstReal(BC, 'External Temperature', Found)
  IF (.NOT. Found) THEN
    ! External Temperature drives the Robin term Elmer assembles; T_ambient only
    ! sets the film temperature for the fluid properties, so a wrong guess here
    ! is a second-order error, not a wrong boundary condition.
    CALL say_once(W_TAMB, 'Missing T_ambient (BC, Constants and External Temperature); using 20 C.')
    T_ambient_C = T_AMBIENT_DEFAULT_C
  END IF

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

  ! Report the degenerate geometry only on the branch that actually uses it --
  ! a still-air face legitimately carries Lc_forced = 0 and channel_gap = 0,
  ! which is what TRAFOLO emits on every face of the current output.
  IF (FcdOrient_in /= FCD_NONE .AND. velocity_mag >= STILL_AIR_VELOCITY_THRESHOLD) THEN
    IF (.NOT. length_is_valid(Lc_forced)) THEN
      WRITE(buf,'(A,ES12.4,A)') 'Forced convection requested but Lc_forced = ', &
        Lc_forced, '; using 10 mm fallback.'
      CALL say_once(W_LCFCD, buf)
    END IF
  END IF
  IF (NatOrient_in == VERTICAL_CHANNEL) THEN
    IF (.NOT. length_is_valid(channel_gap) .OR. .NOT. length_is_valid(channel_length)) THEN
      WRITE(buf,'(A,ES12.4,A,ES12.4,A)') 'Vertical-channel face with channel_gap = ', &
        channel_gap, ' and channel_length = ', channel_length, &
        '; using 10 mm fallback(s).'
      CALL say_once(W_CHAN, buf)
    END IF
  END IF
  IF (theta_deg < 0.0_dp .OR. theta_deg > THETA_MAX_DEG) THEN
    WRITE(buf,'(A,F9.3,A)') 'theta_deg = ', theta_deg, &
      ' outside [0,89]; clamped (g_eff = g*cos(theta) must stay positive).'
    CALL say_once(W_THETA, buf)
  END IF

  ! SIF Temperature variable and T_ambient are both Celsius; convert to Kelvin
  ! before the core (fluid-property correlations require absolute temperature).
  T_amb_K = clamp(T_ambient_C + CELSIUS_TO_KELVIN, T_SURF_MIN_K, T_SURF_MAX_K)

  T_raw = T_surface + CELSIUS_TO_KELVIN
  IF (T_raw /= T_raw) THEN
    ! NaN in the temperature field means the previous linear solve diverged.
    ! Fall back to the ambient so the assembly stays finite and the next iterate
    ! has a chance to recover; aborting here guarantees no result.
    CALL say_once(W_TNAN, 'Temperature field is NaN (previous solve diverged); ' // &
                          'evaluating h at ambient for this node.')
    T_surf_K = T_amb_K
  ELSE
    IF (T_raw > T_SURF_MAX_K) THEN
      WRITE(buf,'(A,F12.1,A)') 'Surface temperature reached ', T_raw - CELSIUS_TO_KELVIN, &
        ' C, above the 1226.85 C validity limit of the air-property fits.' // &
        ' h clamped at its 1500 K value; treat the temperature field as a' // &
        ' lower bound on a thermal runaway, not a prediction.'
      CALL say_once(W_THI, buf)
    ELSE IF (T_raw < T_SURF_MIN_K) THEN
      WRITE(buf,'(A,F12.1,A)') 'Surface temperature reached ', T_raw - CELSIUS_TO_KELVIN, &
        ' C, below the -173.15 C validity limit; h clamped at its 100 K value.'
      CALL say_once(W_TLO, buf)
    END IF
    T_surf_K = clamp(T_raw, T_SURF_MIN_K, T_SURF_MAX_K)
  END IF

  h_raw = htc_evaluate_core(fluid_id, T_surf_K, T_amb_K, NatOrient_in, FcdOrient_in, &
                            Lc_nat, Lc_forced, theta_deg, velocity_mag, &
                            channel_gap, channel_length)

  ! Belt and braces. htc_evaluate_core is total and already clamps, so this only
  ! fires if a future edit breaks that guarantee -- but the cost of it firing
  ! silently is an aborted solve, so the check stays.
  IF (h_raw /= h_raw .OR. h_raw < H_MIN .OR. h_raw > H_MAX) THEN
    WRITE(buf,'(A,ES14.5,A)') 'htc_evaluate_core returned h = ', h_raw, &
      ' outside [1e-3,1e5] W/m2K; clamped. This is a bug, please report.'
    CALL say_once(W_HCLAMP, buf)
    IF (h_raw /= h_raw) THEN
      h_final = H_FALLBACK
    ELSE
      h_final = clamp(h_raw, H_MIN, H_MAX)
    END IF
  ELSE
    h_final = h_raw
  END IF

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

  ! Emit msg through Elmer's Warn (or Info for the banner) the first time slot is
  ! used, then stay quiet. Host association gives access to `warned`.
  SUBROUTINE say_once(slot, msg)
    INTEGER,      INTENT(IN) :: slot
    CHARACTER(*), INTENT(IN) :: msg
    IF (warned(slot)) RETURN
    warned(slot) = .TRUE.
    IF (slot == W_BANNER) THEN
      CALL Info('htc_evaluate', msg, Level=5)
    ELSE
      CALL Warn('htc_evaluate', msg)
    END IF
  END SUBROUTINE say_once

END FUNCTION htc_evaluate
