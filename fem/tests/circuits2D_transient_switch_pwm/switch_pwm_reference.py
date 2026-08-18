"""Independent lumped reference for circuits2D_transient_switch_pwm.

Run with any Python 3 (no third-party packages):  python switch_pwm_reference.py

Circuit (identical to the sif):
    V0 (DC)  --  L, R_coil  --  switch      (series loop)
    s(t) = t*f - phase/360 ;  closed iff  MODULO(s,1) < duty
    R_switch = Ron when closed, Roff when open        (time driven only)

Because the switch resistance does not depend on the solution, the timestep
relation is LINEAR and closed form -- there is no fixed point and no iteration
anywhere in this file, which mirrors the fact that the Elmer run uses
"Steady State Max Iterations = 1":

    L*(I_n - I_{n-1})/dt + (Rcoil + R_switch(t_n))*I_n = V0
    =>  I_n = (V0 + a*I_{n-1}) / (a + Rcoil + R_switch(t_n)),   a = L/dt

R is evaluated at t_n, the END of the step, which is the BDF1 implicit point and
exactly what the Fortran does (it reads GetTime(), the current solution time).

L is not a free parameter. It was calibrated from an Elmer run of this exact
mesh and coil with the switch replaced by a fixed 1 ohm resistor and the
sinusoidal drive of the diode tests, 250 steps.  Solving the BDF1 recurrence
for L step by step gives min 2.395782954 H / max 2.395782960 H /
median 2.3957829546 H, and replaying the recurrence with that constant L
reproduces the Elmer currents to 8.1e-13 A (amplitude 3.17 A).  This is an
independent re-derivation of the inductance the diode tests calibrated
(2.3957829553 H); the two agree to 6.8e-10 H.

EDGE PLACEMENT.  With f = 0.2 Hz and dt = T/100 the grid points sit at
s = n/100.  phase = 1.8 deg = 0.005 of a period puts both edges (closing at
s = 0.005, opening at s = 0.405) exactly halfway between grid points, so the
closed/open decision has half a step of margin and cannot be flipped by
round-off in either code.  Schedule: n mod 100 = 1..40 closed, 41..100 open.

RECORDED OUTPUT (2026-08-17, CPython 3.11.9):

    params: V0=10 f=0.2 duty=0.4 phase=1.8deg Ron=0.1 Roff=100000
            L=2.3957829546 Rcoil=1 nper=100 dt=0.05
    schedule of one period: closed on steps 1..40, open on steps 41..100

    --- third switching period ---     (step 241 is the asserted one)
    n=239  t= 11.950  closed I=5.3397567578e+00 A   R_switch=0.1
    n=240  t= 12.000  closed I=5.4239393974e+00 A   R_switch=0.1
    n=241  t= 12.050  open   I=2.6975967638e-03 A   R_switch=100000
    n=242  t= 12.100  open   I=1.0124304757e-04 A   R_switch=100000
    n=300  t= 15.000  open   I=9.9999000010e-05 A   R_switch=100000
    n=301  t= 15.050  closed I=2.0411418929e-01 A   R_switch=0.1

    open-state Ohm's-law leakage V0/(Rcoil+Roff) = 9.999900001e-05 A
    step 241 is 96.3% carried by the closed-phase current:
       a*I_240 = 2.598916e+02  vs  V0 = 1.000000e+01
    collapse factor across the opening edge: I_240/I_241 = 2010.7

    --- periodicity (same phase in successive periods) ---
    n= 41  I=2.6975774456e-03
    n=141  I=2.6975967638e-03
    n=241  I=2.6975967638e-03      (periodic to 10 digits by the 3rd period)

Elmer reaches 2.69759676E-03 A at step 241, a relative deviation of 4.5e-12.
"""
import math

V0    = 10.0
fsw   = 0.2
duty  = 0.4
phase = 1.8                   # degrees
Ron   = 0.1
Roff  = 1.0e5
L     = 2.3957829546          # calibrated, see docstring
Rcoil = 1.0
nper  = 100
dt    = (1.0 / fsw) / nper
a     = L / dt


def closed(t):
    """Gate state at time t, exactly as the Fortran evaluates it."""
    return math.fmod(t * fsw - phase / 360.0, 1.0) % 1.0 < duty


def Rsw(t):
    return Ron if closed(t) else Roff


def run_bdf1(nsteps):
    """Replay the closed-form BDF1 recurrence.  No iteration of any kind."""
    I, out = 0.0, []
    for n in range(1, nsteps + 1):
        t = n * dt
        I = (V0 + a * I) / (a + Rcoil + Rsw(t))
        out.append((n, t, closed(t), I))
    return out


if __name__ == '__main__':
    print('params: V0=%g f=%g duty=%g phase=%gdeg Ron=%g Roff=%g L=%.10f '
          'Rcoil=%g nper=%d dt=%g' % (V0, fsw, duty, phase, Ron, Roff, L,
                                      Rcoil, nper, dt))
    sol = run_bdf1(301)
    on = [n for n, t, c, I in sol if c and n <= 100]
    print('schedule of one period: closed on steps %d..%d, open on steps %d..%d'
          % (on[0], on[-1], on[-1] + 1, 100))

    print('\n--- third switching period ---')
    for n, t, c, I in sol:
        if n in (239, 240, 241, 242, 300, 301):
            print('n=%3d  t=%7.3f  %-6s I=%.10e A   R_switch=%g'
                  % (n, t, 'closed' if c else 'open', I, Rsw(t)))

    I240 = [I for n, t, c, I in sol if n == 240][0]
    I241 = [I for n, t, c, I in sol if n == 241][0]
    print('\nopen-state Ohm\'s-law leakage V0/(Rcoil+Roff) = %.9e A'
          % (V0 / (Rcoil + Roff)))
    print('step 241 is %.1f%% carried by the closed-phase current:'
          % (100.0 * a * I240 / (V0 + a * I240)))
    print('   a*I_240 = %e  vs  V0 = %e' % (a * I240, V0))
    print('collapse factor across the opening edge: I_240/I_241 = %.1f'
          % (I240 / I241))

    print('\n--- periodicity (same phase in successive periods) ---')
    for n, t, c, I in sol:
        if n in (41, 141, 241):
            print('n=%3d  I=%.10e' % (n, I))
