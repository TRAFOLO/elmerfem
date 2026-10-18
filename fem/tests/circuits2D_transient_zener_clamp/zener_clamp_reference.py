"""Independent lumped reference for circuits2D_transient_zener_clamp.

Run with any Python 3 (no third-party packages):  python zener_clamp_reference.py

Circuit (identical to the sif):
    V(t) = V0*sin(2*pi*f*t)  --  L, R_coil  --  zener      (series loop)
    R(V) = Ron + 0.25*(Roff-Ron)*(1 - tanh((V-Vf)/Vs))*(1 + tanh((V+Vz)/Vs))
    operating point:  V_Z = R(V_Z) * I                     (implicit)

This is the SAME method as ../circuits2D_transient_diode_rectifier/
diode_reference.py -- only the resistance law differs (two knees instead of
one).  Nothing here is shared with Elmer: each timestep is solved by scanning
the scalar residual g(V) = R(V)*I(V) - V for sign changes and bisecting, never
by Picard.  The scan doubles as a diagnostic: it reports how many fixed points
the timestep has, which matters more here than for the plain diode because a
clamp is in principle capable of having several.

WHY Vs = 0.1 IS NOW AFFORDABLE (it was 1.0 while the loop was Picard)
---------------------------------------------------------------------
The knee width Vs is what sets the clamp accuracy: the flank of the smoothed
law sits at

    V_clamp = -Vz - 0.5*Vs*ln(Roff/R)

so the overshoot beyond the nominal breakdown voltage is proportional to Vs
(here 0.5*0.1*ln(1e5/9.86) = 0.461 V, and the exact BDF1 solve below lands on
-5.4617 V for Vz = 5).  Vs could not be made small while the coupled loop was
PICARD -- R frozen at the (relaxed) previous iterate, linear solve, repeat --
because on this series loop that is the scalar map

    V <- G(V) = R(V)*S/(D + R(V)),   S = Vsrc(t) + a*I_{n-1},  D = a + Rcoil,
    G'(V) = S*D*R'(V)/(D + R(V))**2 ,                          a = L/dt

which converges under relaxation w only while w < 2/(1 + |G'|).  In reverse
conduction the branch is held at |V| ~ Vz by a resistance sitting on the
exponential flank, where R' = 2R/Vs, so

    |G'|  =  2 * |V_clamp| * D / ( Vs * (D + R) )   ~   2*Vz/Vs   for R << D,

i.e. the stiffness is inversely proportional to exactly the parameter that
controls the accuracy.  At Vs = 0.1 this run has max|G'| = 101.7 (computed from
the converged trajectory below), needing w < 0.020, and the forward lobe then
crawls at |1 - w*(1+|G'|)| ~ 0.93 per iteration.  That is why the Picard
version of this test had to open Vs up to 1.0 and settled at a clamp voltage of
-8.85 V for a 5 V part -- an 77% overshoot.

Slice 4 stamps the analytic dR/dV instead (NEWTON, the default; see
'Diode Newton Linearization' in CircuitsAndDynamics.F90), so the flank is no
longer a stability limit and Vs goes back to 0.1.  MEASURED in Elmer over the
full 75-step run, Vs = 0.1, 'Nonlinear Circuit Residual Tolerance = 1e-8',
'Steady State Min Iterations = 3', 'Steady State Max Iterations = 100'
(timesteps that exhaust the limit / min / median / max coupled iterations,
wall clock):

    w = 1.0  ->  73/75 exhausted   4/100/100      80.3 s   diverges
    w = 0.7  ->   0/75 exhausted   4/  7/ 29       2.6 s   <== used here
    w = 0.5  ->   1/75 exhausted   4/ 11/100       5.2 s
    w = 0.3  ->   0/75 exhausted   5/ 19/ 92      18.3 s

w = 0.7 is both the largest value of the sweep that is stable and the cheapest.
Newton still needs SOME damping here (w = 1.0 chatters) because the expansion
point is the relaxed iterate and the two knees are only 5.7 V apart at a knee
width of 0.1 V; the freewheel sibling, whose diode has a single knee, runs at
w = 1.0 with no relaxation at all.  For the record, the Picard version of this
test needed w = 0.10 AND Vs = 1.0 and took 14.9 s with median 8 iterations.

WHAT THIS BUYS, IN ONE LINE: the clamp moves from -8.853 V to -5.462 V against
a nominal -Vz = -5 V, i.e. the overshoot drops from 77% to 9.2%, and the
residual overshoot is now a pure knee-width effect (4.6*Vs) rather than a
stability compromise.

L is not a free parameter.  It was calibrated from an Elmer run of this exact
mesh and coil with the nonlinear element replaced by a fixed 1 ohm resistor and
the same sinusoidal drive, 250 steps.  Solving the BDF1 recurrence
L*(I_n - I_{n-1})/dt = V(t_n) - R_tot*I_n for L step by step gives
min 2.395782954 H / max 2.395782960 H / median 2.3957829546 H; a golden-section
fit that minimises the worst replay error over the whole run lands on the same
2.3957829546 H and replays the Elmer currents to 5.95e-13 A on a 3.172770 A
amplitude.  (The diode tests calibrated 2.3957829553 H by the same recipe; the
two agree to 7e-10 H, and replaying with the diode value costs only 7.1e-10 A.)

RECORDED OUTPUT (2026-08-17, CPython 3.11.9):

    params: Vf=0.7 Vz=5 Ron=0.1 Roff=100000 Vs=0.1 L=2.3957829546 Rcoil=1 V0=10
            f=0.2 nper=100 dt=0.05
    R(V) at the characteristic points:
      R(+9)=1.000000e-01  R(+5)=1.000000e-01  R(+0.7)=5.000005e+04
      R(-2.15)=1.000000e+05  R(-5)=5.000005e+04  R(-9)=1.000000e-01
    max number of fixed points found in any step: 1

    --- deep reverse, period 1 (t = 0.75 T, source at its negative peak) ---
    step  75  t=3.7500  Vsrc=-10.000000  I=-5.537790488e-01 A
              Vz=-5.461713 V  Rz=9.862622 ohm
       neighbour step 74 t=3.7000 I=-4.706223577e-01
       neighbour step 76 t=3.8000 I=-6.346926862e-01
       plain single-knee diode at the same instant: I=-9.999890537e-05 A
       -> the clamp carries 5538 times more reverse current
       source follows to -10.000 V but the branch stops at -5.462 V
       the same instant one period later (step 175) gives -5.578004664e-01 A,
       i.e. period-to-period drift of 0.73% -- the run is asserted
       step by step against this reference, not against a periodic state

    --- forward peak, period 1 ---
    step  41  t=2.0500  Vsrc=+5.358268  I=3.579523441 A  Vz=+1.340298 V
              Rz=0.374435 ohm

    --- timestep refinement (continuous limit) ---
    nper=100  I(3.75)=-0.553779049  I(2.05)=3.579523441
    nper=200  I(3.75)=-0.548266104  I(2.05)=3.591141930
    nper=400  I(3.75)=-0.548153499  I(2.05)=3.596960276
    Richardson (2*fine-coarse) at t=3.75: -0.548040895  (1.0% from nper=100)
    Richardson (2*fine-coarse) at t=2.05: 3.602778621  (0.6% from nper=100)

Note how much crisper the law is at Vs = 0.1 than the Vs = 1.0 version this
replaces: R(+5) and R(-9) are now exactly Ron and R(-2.15) is exactly Roff,
whereas at Vs = 1.0 they were 18.5, 33.6 and 9.93e4 -- the knees no longer
bleed into the plateaus.

The gap to the continuous limit is BDF1 time-discretisation error, present in
Elmer and in this reference alike; it cancels in the comparison because both use
the same BDF1 grid.  The asserted instant is the deep reverse point, where the
two-knee law does something the single-knee law cannot: the plain diode leaks
-1.0e-04 A and lets its branch follow the source to -10 V, while the zener
conducts -0.554 A at 9.86 ohm and holds the branch 4.54 V off the source.

WHAT ELMER PRODUCES (build of 2026-08-17, Vs = 0.1, w = 0.7, Newton on):
    I       = -5.537790487613e-01 A   (this reference: -5.537790488e-01)
    v_zener = -5.461713292394e+00 V   (this reference: -5.461713292)
i.e. 11 significant digits on both, against a scheme that shares no code with
Elmer and never uses Picard.
"""
import math

# ---------------------------------------------------------------- parameters
V0    = 10.0
f     = 0.2
L     = 2.3957829546          # calibrated, see docstring
Rcoil = 1.0
Vf    = 0.7
Vz    = 5.0
Ron   = 0.1
Roff  = 1.0e5
Vs    = 0.1                   # see "WHY Vs = 0.1 IS NOW AFFORDABLE" above


def Rz(v):
    """Two-knee (zener) law, exactly the expression the Fortran evaluates."""
    return Ron + 0.25 * (Roff - Ron) * (1.0 - math.tanh((v - Vf) / Vs)) \
                                     * (1.0 + math.tanh((v + Vz) / Vs))


def Rd(v):
    """Single-knee law of the plain diode, for the side-by-side comparison.
    Vf and Vs are the ones the four existing diode tests use (0 and 0.1), so
    this reproduces circuits2D_transient_diode_blocking exactly."""
    return Ron + 0.5 * (Roff - Ron) * (1.0 - math.tanh(v / 0.1))


def Vsrc(t):
    return V0 * math.sin(2.0 * math.pi * f * t)


def _grid():
    """Coarse scan over the whole plausible range plus a fine scan across each
    of the two transitions (width ~Vs around Vf and around -Vz)."""
    g = [(-5.0 * V0) + (10.0 * V0) * k / 8000.0 for k in range(8001)]
    for c in (0.0, Vf, -Vz):
        g += [c - 4.0 + 8.0 * k / 40000.0 for k in range(40001)]
    return sorted(set(g))

GRID = _grid()


def _bisect(g, A, B, gA):
    for _ in range(200):
        M = 0.5 * (A + B)
        gm = g(M)
        if gA * gm <= 0.0:
            B = M
        else:
            A, gA = M, gm
        if B - A < 1e-14 * max(1.0, abs(M)):
            break
    return 0.5 * (A + B)


def step_bdf1(Vn, Iprev, a, R=Rz, scan=True):
    """One implicit BDF1 step.  Returns (I_n, V_branch, n_fixed_points).

    L*(I_n - I_{n-1})/dt + Rcoil*I_n + V = V(t_n),   V = R(V)*I_n
    =>  I_n(V) = (V(t_n) + a*I_{n-1} - V) / (a + Rcoil),   a = L/dt
    =>  g(V)   = R(V)*I_n(V) - V = 0
    """
    def In(v):
        return (Vn + a * Iprev - v) / (a + Rcoil)

    def g(v):
        return R(v) * In(v) - v

    if not scan:                       # single-root fast path
        lo, hi = -5.0 * V0, 5.0 * V0
        v = _bisect(g, lo, hi, g(lo))
        return In(v), v, 1

    roots = []
    vprev, gprev = GRID[0], g(GRID[0])
    for v in GRID[1:]:
        gv = g(v)
        if gprev == 0.0:
            roots.append(vprev)
        elif gprev * gv < 0.0:
            roots.append(_bisect(g, vprev, v, gprev))
        vprev, gprev = v, gv
    if not roots:
        raise RuntimeError('no root for Vn=%g Iprev=%g' % (Vn, Iprev))
    return In(roots[0]), roots[0], len(roots)


def run_bdf1(nsteps, nper=100, R=Rz, scan=True):
    dt = (1.0 / f) / nper
    a = L / dt
    I, out, maxroots = 0.0, [], 1
    for n in range(1, nsteps + 1):
        t = n * dt
        I, v, nr = step_bdf1(Vsrc(t), I, a, R=R, scan=scan)
        maxroots = max(maxroots, nr)
        out.append((n, t, Vsrc(t), I, v, R(v)))
    return out, maxroots


def at_time(tend, nper, R=Rz):
    n = int(round(tend / ((1.0 / f) / nper)))
    return run_bdf1(n, nper, R=R, scan=False)[0][-1][3]


if __name__ == '__main__':
    T, nper = 1.0 / f, 100
    dt = T / nper
    print('params: Vf=%g Vz=%g Ron=%g Roff=%g Vs=%g L=%.10f Rcoil=%g V0=%g '
          'f=%g nper=%d dt=%g' % (Vf, Vz, Ron, Roff, Vs, L, Rcoil, V0, f,
                                  nper, dt))
    print('R(V) at the characteristic points:')
    print('  R(+9)=%e  R(+5)=%e  R(+0.7)=%e' % (Rz(9.0), Rz(5.0), Rz(Vf)))
    print('  R(-2.15)=%e  R(-5)=%e  R(-9)=%e' % (Rz(-2.15), Rz(-Vz), Rz(-9.0)))

    sol, maxroots = run_bdf1(76, nper)
    print('max number of fixed points found in any step: %d' % maxroots)

    rv = sol[74]                       # step 75, t = 3.75 s = 0.75 T
    print('\n--- deep reverse, period 1 (t = 0.75 T, source at its negative '
          'peak) ---')
    print('step %3d  t=%.4f  Vsrc=%+.6f  I=%.9e A  Vz=%+.6f V  Rz=%.6f ohm'
          % (rv[0], rv[1], rv[2], rv[3], rv[4], rv[5]))
    for s in sol:
        if s[0] in (74, 76):
            print('   neighbour step %d t=%.4f I=%.9e' % (s[0], s[1], s[3]))
    dsol, _ = run_bdf1(75, nper, R=Rd, scan=False)
    print('   plain single-knee diode at the same instant: I=%.9e A'
          % dsol[-1][3])
    print('   -> the clamp carries %.0f times more reverse current'
          % abs(rv[3] / dsol[-1][3]))
    print('   source follows to %.3f V but the branch stops at %.3f V'
          % (rv[2], rv[4]))
    later, _ = run_bdf1(175, nper, scan=False)
    print('   the same instant one period later (step 175) gives %.9e A,'
          % later[-1][3])
    print('   i.e. period-to-period drift of %.2f%% -- the run is asserted'
          % (100.0 * abs(later[-1][3] - rv[3]) / abs(rv[3])))
    print('   step by step against this reference, not against a periodic state')

    pk = max((s for s in sol if s[1] <= T), key=lambda s: s[3])
    print('\n--- forward peak, period 1 ---')
    print('step %3d  t=%.4f  Vsrc=%+.6f  I=%.9f A  Vz=%+.6f V  Rz=%.6f ohm'
          % (pk[0], pk[1], pk[2], pk[3], pk[4], pk[5]))

    print('\n--- timestep refinement (continuous limit) ---')
    vals = {}
    for np_ in (100, 200, 400):
        vals[np_] = (at_time(0.75 * T, np_), at_time(pk[1], np_))
        print('nper=%-4d I(%.2f)=%.9f  I(%.2f)=%.9f'
              % (np_, 0.75 * T, vals[np_][0], pk[1], vals[np_][1]))
    for k, lbl in ((0, 0.75 * T), (1, pk[1])):
        rich = 2 * vals[400][k] - vals[200][k]
        print('Richardson (2*fine-coarse) at t=%.2f: %.9f  (%.1f%% from '
              'nper=100)' % (lbl, rich,
                             100.0 * abs(rich - vals[100][k]) / abs(vals[100][k])))
