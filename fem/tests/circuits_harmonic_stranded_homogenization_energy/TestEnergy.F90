! Check the harmonic energy of a homogenized coil against the circuit: with peak
! phasors W = VF Q/omega, Q = Im(V conj(I))/2 of the coil and VF the share of
! the coil the mesh holds. Energy and coenergy agree for linear materials. Also
! built by circuits_harmonic_foilsheet_energy.
SUBROUTINE CheckHarmonicEnergy(Model, Solver, dt, Transient)
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient, Found
  REAL(KIND=dp) :: w, co, vre, vim, ire, iim, omega, vf, circuit_w, err

  w = GetConstReal(Model % Simulation, 'res: Magnetic Field Energy', Found)
  IF (.NOT. Found) CALL Fatal('CheckHarmonicEnergy', 'Magnetic energy missing')
  co = GetConstReal(Model % Simulation, 'res: Magnetic Coenergy', Found)
  IF (.NOT. Found) CALL Fatal('CheckHarmonicEnergy', 'Magnetic coenergy missing')
  vre = GetConstReal(Model % Simulation, 'res: v_component(1) re', Found)
  IF (Found) vim = GetConstReal(Model % Simulation, 'res: v_component(1) im', Found)
  IF (Found) ire = GetConstReal(Model % Simulation, 'res: i_component(1) re', Found)
  IF (Found) iim = GetConstReal(Model % Simulation, 'res: i_component(1) im', Found)
  IF (.NOT. Found) CALL Fatal('CheckHarmonicEnergy', 'Circuit voltage or current missing')
  omega = GetAngularFrequency(Found=Found)
  IF (.NOT. Found) CALL Fatal('CheckHarmonicEnergy', 'Angular Frequency or Frequency missing')
  vf = GetConstReal(Model % Components(1) % Values, 'Circuit Equation Voltage Factor', Found)
  IF (.NOT. Found) vf = 1._dp

  circuit_w = vf * 0.5_dp * (vim*ire - vre*iim) / omega
  err = ABS(w - circuit_w) / ABS(circuit_w)
  WRITE(Message,'(A,3ES22.12)') 'Energy, circuit energy, relative error: ', w, circuit_w, err
  CALL Info('CheckHarmonicEnergy', Message, Level=1)
  ! Written so that a NaN fails too.
  IF (.NOT. (ABS(w-co) <= 1.e-12_dp * ABS(w))) &
      CALL Fatal('CheckHarmonicEnergy', 'Linear energy and coenergy differ')
  IF (.NOT. (err <= 1.e-5_dp)) &
      CALL Fatal('CheckHarmonicEnergy', 'Stored energy does not match circuit reactive power')
END SUBROUTINE CheckHarmonicEnergy
