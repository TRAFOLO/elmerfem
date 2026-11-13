! Check the harmonic field power against the circuit: with peak phasors
! P = VF Re(V conj(I))/2 of the coil, VF the share of the coil the mesh holds.
! The magnetic loss of the core has to be a real part of it for the check to mean anything.
SUBROUTINE CheckHarmonicPower(Model, Solver, dt, Transient)
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient, Found
  REAL(KIND=dp) :: p, joule, vre, vim, ire, iim, vf, circuit_p, err

  p = GetConstReal(Model % Simulation, 'res: Eddy current power', Found)
  IF (.NOT. Found) CALL Fatal('CheckHarmonicPower', 'Eddy current power missing')
  joule = GetConstReal(Model % Simulation, 'res: joule loss', Found)
  IF (.NOT. Found) CALL Fatal('CheckHarmonicPower', 'Joule loss missing')
  vre = GetConstReal(Model % Simulation, 'res: v_component(1) re', Found)
  IF (Found) vim = GetConstReal(Model % Simulation, 'res: v_component(1) im', Found)
  IF (Found) ire = GetConstReal(Model % Simulation, 'res: i_component(1) re', Found)
  IF (Found) iim = GetConstReal(Model % Simulation, 'res: i_component(1) im', Found)
  IF (.NOT. Found) CALL Fatal('CheckHarmonicPower', 'Circuit voltage or current missing')
  vf = GetConstReal(Model % Components(1) % Values, 'Circuit Equation Voltage Factor', Found)
  IF (.NOT. Found) vf = 1._dp

  circuit_p = vf * 0.5_dp * (vre*ire + vim*iim)
  err = ABS(p - circuit_p) / ABS(circuit_p)
  WRITE(Message,'(A,4ES22.12)') 'Power, Joule part, circuit power, relative error: ', p, joule, circuit_p, err
  CALL Info('CheckHarmonicPower', Message, Level=1)
  ! Written so that a NaN fails too.
  IF (.NOT. (p - joule >= 0.2_dp * p)) &
      CALL Fatal('CheckHarmonicPower', 'The magnetic loss is not a fifth of the power')
  IF (.NOT. (err <= 1.e-3_dp)) &
      CALL Fatal('CheckHarmonicPower', 'Field power does not match circuit active power')
END SUBROUTINE CheckHarmonicPower
