! Check the full linear model's energy/coenergy and its reactive terminal power.
SUBROUTINE CheckEnergy(Model, Solver, dt, Transient)
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient, Found
  INTEGER, SAVE :: step = 0
  REAL(KIND=dp), SAVE :: ew = 0._dp, vs = 0._dp, vc = 0._dp, is = 0._dp, ic = 0._dp
  REAL(KIND=dp) :: w, co, v, current, phase, omega, omega_bdf, circuit_w, err

  step = step + 1
  w = GetConstReal(Model % Simulation, 'res: Magnetic Field Energy', Found)
  IF (.NOT. Found) CALL Fatal('CheckEnergy', 'Magnetic energy missing')
  co = GetConstReal(Model % Simulation, 'res: Magnetic Coenergy', Found)
  IF (.NOT. Found) CALL Fatal('CheckEnergy', 'Magnetic coenergy missing')
  ! A guard, not the check: CalcFields sets the litz block's coenergy to its
  ! energy, so this only fails on a broken total or a nonlinear material.
  ! The comparisons are written so that a NaN fails too.
  IF (.NOT. (ABS(w-co) <= 1.e-12_dp * MAX(ABS(w), 1.e-20_dp))) &
      CALL Fatal('CheckEnergy', 'Linear full-state energy and coenergy differ')

  ! Discard the first 50 Hz cycle, then Fourier analyze the next 40 samples.
  omega = 100._dp * ACOS(-1._dp)
  IF (step > 40) THEN
    v = GetConstReal(Model % Simulation, 'res: v_component(1)', Found)
    IF (.NOT. Found) CALL Fatal('CheckEnergy', 'Circuit voltage missing')
    current = GetConstReal(Model % Simulation, 'res: i_component(1)', Found)
    IF (.NOT. Found) CALL Fatal('CheckEnergy', 'Circuit current missing')
    phase = omega * GetTime()
    ew = ew + w / 40._dp
    vs = vs + v * SIN(phase) / 20._dp
    vc = vc + v * COS(phase) / 20._dp
    is = is + current * SIN(phase) / 20._dp
    ic = ic + current * COS(phase) / 20._dp
  END IF

  ! For sinusoidal steady state <W> = L I_peak^2/4 = Q/(2 omega).
  ! The mesh covers half the circuit (Circuit Equation Voltage Factor = 0.5).
  ! Use the imaginary part of the BDF2 differentiation symbol to remove its
  ! leading frequency error; the remaining real part is numerical dissipation.
  omega_bdf = (4._dp*SIN(omega*dt) - SIN(2._dp*omega*dt)) / (2._dp*dt)
  circuit_w = 0.5_dp * (vc*is - vs*ic) / (4._dp*omega_bdf)
  CALL ListAddConstReal(Model % Simulation, 'res: Mean Stored Energy', ew)
  CALL ListAddConstReal(Model % Simulation, 'res: Circuit Stored Energy', circuit_w)
  IF (step == 80) THEN
    err = ABS(ew-circuit_w) / ABS(circuit_w)
    WRITE(Message,'(A,4ES22.12)') 'Final energy, mean energy, circuit energy, relative error: ', &
        w, ew, circuit_w, err
    CALL Info('CheckEnergy', Message, Level=1)
    IF (.NOT. (err <= 1.e-5_dp)) CALL Fatal('CheckEnergy', 'Stored energy does not match circuit reactive power')
  END IF
END SUBROUTINE CheckEnergy
