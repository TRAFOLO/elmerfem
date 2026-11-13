"""Independent lumped reference for circuits2D_transient_diode_{rectifier,blocking}.

Run with any Python 3 (no third-party packages):  python diode_reference.py

Circuit (identical to the sif):
    V(t) = V0*sin(2*pi*f*t)  --  L, R_coil  --  diode      (series loop)
    R(V) = Ron + 0.5*(Roff-Ron)*(1 - tanh(V/Vsm)),  V = diode voltage
    operating point:  V_D = R(V_D) * I                     (implicit)

L is not a free parameter: it was calibrated from an Elmer run of this exact
case with the diode replaced by a fixed 1 ohm resistor.  Solving the BDF1
recurrence  L*(I_n - I_{n-1})/dt = V(t_n) - R_tot*I_n  for L step by step over
250 steps gave min 2.395726 H / max 2.396014 H / median 2.3957830 H, and
replaying the recurrence with that constant L reproduces the Elmer currents to
1.3e-7 A (amplitude 2.73 A).  The continuous-time amplitude
V0/sqrt(R^2+(wL)^2) = 2.7667 A is 1.44% above the Elmer/BDF1 amplitude of
2.7268 A -- that gap is the BDF1 damping at 100 steps/period, which is why the
asserted reference below is evaluated on the SAME BDF1 grid Elmer uses.

Nothing here is shared with Elmer: each timestep is solved by scanning the
scalar residual g(V) = R(V)*I(V) - V for sign changes and bisecting, never by
Picard.  The scan doubles as a diagnostic -- it reports how many fixed points
the timestep has.

RECORDED OUTPUT (2026-08-16, CPython 3.11.9):

    params: Ron=0.1 Roff=100000 Vsm=0.1 L=2.3957829553 Rcoil=1 V0=10 f=0.2
            nper=100 dt=0.05
    max number of fixed points found in any step: 1

    --- forward peak, period 2 ---
    step 142  t=7.1000  Vsrc=+4.817537  I=3.946375414 A  Vd=+0.703095 V
              Rd=0.178162 ohm
       neighbour step 141 t=7.0500 I=3.942867863
       neighbour step 143 t=7.1500 I=3.938372286

    --- deep reverse (t = 1.75 T) ---
    step 175  t=8.7500  Vsrc=-10.000000  I=-9.999890537e-05 A
              Vd=-9.999891 V  Rd=100000 ohm
       pure-resistive leakage Vsrc/(Roff+Rcoil) = -9.999900001e-05 A
       |wL*I| = 3.011e-04 V vs |Vsrc| = 10.000 V  -> the inductive term is
       3e-5 of the source, so the leakage is the designed Ohm's-law value

    --- timestep refinement (continuous limit) ---
    nper=100  I(7.10)=3.946375414  I(8.75)=-9.999890537e-05
    nper=200  I(7.10)=3.960334320  I(8.75)=-9.999895263e-05
    nper=400  I(7.10)=3.967331250  I(8.75)=-9.999897627e-05
    Richardson (2*fine-coarse) at t=7.10: 3.974328181  (0.7% above nper=100)

The 0.7% gap at the forward peak is BDF1 time-discretisation error, present in
Elmer and in this reference alike; it cancels in the comparison because both
use the same grid.  The reverse leakage is insensitive to it (7 digits stable).
A plain explicit RK4 integration of the continuous ODE was
tried first and is NOT usable here: in the blocking state the system time
constant is L/(Roff+Rcoil) = 2.4e-5 s, so RK4 is unstable at any step size
that is affordable over 8.75 s, and it diverged.  The refinement study above
replaces it.
"""
import math, sys

# ---------------------------------------------------------------- parameters
V0    = 10.0
f     = 0.2
L     = 2.3957829553          # calibrated, see docstring
Rcoil = 1.0
Ron   = 0.1
Roff  = 1.0e5
Vsm   = 0.1


def Rd(v):
    return Ron + 0.5 * (Roff - Ron) * (1.0 - math.tanh(v / Vsm))


def Vsrc(t):
    return V0 * math.sin(2.0 * math.pi * f * t)


def _grid():
    """Coarse scan over the whole plausible range plus a fine scan across the
    diode transition (width ~Vsm around 0)."""
    g = [(-5.0 * V0) + (10.0 * V0) * k / 8000.0 for k in range(8001)]
    w = max(20.0 * Vsm, 2.0)
    g += [-w + 2.0 * w * k / 40000.0 for k in range(40001)]
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


def step_bdf1(Vn, Iprev, a, scan=True):
    """One implicit BDF1 step.  Returns (I_n, V_D, n_fixed_points).

    L*(I_n - I_{n-1})/dt + Rcoil*I_n + V_D = V(t_n),  V_D = R(V_D)*I_n
    =>  I_n(V) = (V(t_n) + a*I_{n-1} - V) / (a + Rcoil),  a = L/dt
    =>  g(V)   = R(V)*I_n(V) - V = 0
    """
    def In(v):
        return (Vn + a * Iprev - v) / (a + Rcoil)

    def g(v):
        return Rd(v) * In(v) - v

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


def run_bdf1(nsteps, nper=100, scan=True):
    dt = (1.0 / f) / nper
    a = L / dt
    I, out, maxroots = 0.0, [], 1
    for n in range(1, nsteps + 1):
        t = n * dt
        I, v, nr = step_bdf1(Vsrc(t), I, a, scan=scan)
        maxroots = max(maxroots, nr)
        out.append((n, t, Vsrc(t), I, v, Rd(v)))
    return out, maxroots


def at_time(tend, nper):
    n = int(round(tend / ((1.0 / f) / nper)))
    return run_bdf1(n, nper, scan=False)[0][-1][3]


if __name__ == '__main__':
    T, nper = 1.0 / f, 100
    dt = T / nper
    print('params: Ron=%g Roff=%g Vsm=%g L=%.10f Rcoil=%g V0=%g f=%g nper=%d '
          'dt=%g' % (Ron, Roff, Vsm, L, Rcoil, V0, f, nper, dt))
    sol, maxroots = run_bdf1(175, nper)
    print('max number of fixed points found in any step: %d' % maxroots)

    pk = max((s for s in sol if s[1] > T + 1e-12), key=lambda s: s[3])
    print('\n--- forward peak, period 2 ---')
    print('step %d  t=%.4f  Vsrc=%+.6f  I=%.9f A  Vd=%+.6f V  Rd=%.6f ohm'
          % (pk[0], pk[1], pk[2], pk[3], pk[4], pk[5]))
    for s in sol:
        if s[0] in (pk[0] - 1, pk[0] + 1):
            print('   neighbour step %d t=%.4f I=%.9f' % (s[0], s[1], s[3]))

    rv = sol[int(round(1.75 * T / dt)) - 1]
    print('\n--- deep reverse (t = 1.75 T) ---')
    print('step %d  t=%.4f  Vsrc=%+.6f  I=%.9e A  Vd=%+.6f V  Rd=%.6g ohm'
          % (rv[0], rv[1], rv[2], rv[3], rv[4], rv[5]))
    print('   pure-resistive leakage Vsrc/(Roff+Rcoil) = %.9e A'
          % (rv[2] / (Roff + Rcoil)))
    print('   |wL*I| = %.3e V vs |Vsrc| = %.3f V'
          % (abs(2 * math.pi * f * L * rv[3]), abs(rv[2])))

    print('\n--- timestep refinement (continuous limit) ---')
    vals = {}
    for np_ in (100, 200, 400):
        vals[np_] = (at_time(pk[1], np_), at_time(1.75 * T, np_))
        print('nper=%-4d I(%.2f)=%.9f  I(%.2f)=%.9e'
              % (np_, pk[1], vals[np_][0], 1.75 * T, vals[np_][1]))
    rich = 2 * vals[400][0] - vals[200][0]
    print('Richardson (2*fine-coarse) at t=%.2f: %.9f  (%.1f%% above nper=100)'
          % (pk[1], rich, 100.0 * (rich - vals[100][0]) / vals[100][0]))
