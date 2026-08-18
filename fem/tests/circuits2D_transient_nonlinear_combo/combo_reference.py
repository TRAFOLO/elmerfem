"""Independent lumped reference for circuits2D_transient_nonlinear_combo.

Run with any Python 3 (no third-party packages):  python combo_reference.py

Circuit (identical to the sif) -- the ZENER-CLAMPED BUCK (demagnetization)
cell, the buck cell of ../circuits2D_transient_switch_diode_freewheel with the
single freewheel diode replaced by a diode in series with a zener:

           +--- switch ---+-------------+
           |    S       A |             |
  V0 (DC)--+              |           L, Rcoil
           |              |             |
           |        zener |             |
           +----- M ------+             |
           |  diode                     |
           +----------------------------+
                   G

    switch:  R_s(t) = Ron_s if MODULO(t*fsw - phase/360, 1) < duty else Roff_s
    diode :  R_d(V) = Ron + 0.5*(Roff-Ron)*(1 - tanh((V-Vf)/Vs))
    zener :  R_z(V) = Ron + 0.25*(Roff-Ron)*(1-tanh((V-Vf)/Vs))*(1+tanh((V+Vz)/Vs))
    coil  :  L*(i1_n - i1_{n-1})/dt + Rcoil*i1_n = v1

The zener points the OPPOSITE way round in the branch, so the clamp current
runs FORWARD through the diode and BACKWARD through the zener, i.e. the zener
is in reverse BREAKDOWN and the branch holds node A at
-(Vz + Vf + I*(Ron_d + Ron_z)) rather than at one forward drop.

THE STEP, AND WHY IT IS STILL ONE SCALAR EQUATION.  The clamp branch holds two
solution-dependent elements that share one current but have separate voltages,
so a naive formulation has two unknowns.  It reduces to one because each
element's characteristic

    i(V) = V / R(V),    di/dV = (R - V*R')/R**2 > 0

is STRICTLY MONOTONE for both laws (R > 0 everywhere, and V*R' <= 0 on every
branch of either law), hence invertible.  Taking the diode's own voltage v3 as
the unknown:

    i_clamp = v3 / R_d(v3)                       (diode characteristic)
    v4      = R_z-inverse(-i_clamp)              (zener carries -i_clamp)
    v1      = v4 - v3                            (KVL of the clamp loop)
    i_coil  = (a*i_coil_prev + v1) / (a + Rcoil),        a = L/dt
    i_sw    = (V0 - v1) / R_s(t)                 (KVL of the source loop)
    g(v3)   = i_coil - i_sw - i_clamp            (KCL at node A)

g is strictly decreasing in v3, so the step has exactly one solution.  Nothing
here is shared with Elmer: the inner inverse and the outer root are both found
by bisection, never by Picard, and the outer residual is additionally SCANNED
over a grid every timestep to count sign changes and to re-check that g really
is monotone -- the reported "max number of fixed points" is that count.

R_s is evaluated at t_n, the END of the step, which is the BDF1 implicit point
and exactly what the Fortran does (it reads GetTime()).  With fsw = 0.2 Hz,
dt = Tsw/50 and phase = 3.6 deg = 0.01 of a period both switching edges sit
exactly halfway between grid points, so the closed/open decision has half a
step of margin in either code.  Schedule: steps n mod 50 = 1..40 closed,
41..50 open.

DUTY IS 0.8, not the 0.4 of the freewheel sibling: a 5.7 V clamp empties the
coil about five times faster than a 0.7 V freewheel, so at duty 0.4 the coil
current would hit zero 15 steps into every OFF phase and the cell would run in
discontinuous conduction, where "the current at the end of an OFF phase"
asserts nothing.  At duty 0.8 the cell stays in continuous conduction.

L is not a free parameter.  It was re-derived for this test from an Elmer run
of this exact mesh and coil with a fixed 1 ohm resistor in place of the
nonlinear elements and a 10*sin(2*pi*0.2*t) drive, 250 steps at dt = 0.05, so
the loop is linear and Rtot = 2.0 exactly.  Solving the BDF1 recurrence
L*(I_n - I_{n-1})/dt = V(t_n) - Rtot*I_n for L step by step gives
min 2.3957829537 H / max 2.3957829601 H / median 2.3957829546 H, and a
golden-section fit minimising the worst replay error over the whole run lands
on the same 2.3957829546 H with a replay error of 1.7e-11 A on a 3.172770 A
amplitude.  That is the value the zener and freewheel siblings calibrated; the
diode tests' 2.3957829553 H replays this run to 7.1e-10 A.

NO TIMESTEP-REFINEMENT STUDY is reported here, unlike the zener sibling: the
switch phase is chosen so that the edges fall midway between GRID POINTS, so
halving dt puts them exactly ON grid points and changes the schedule rather
than refining it.  Refining would need the phase re-chosen with it, which is a
different run, not the same run on a finer grid.

RECORDED OUTPUT (2026-08-18, CPython 3.11.9):

    source V0=10 | switch fsw=0.2 duty=0.8 phase=3.6deg Ron=0.1 Roff=100000
    diode  Vf=0.7 Ron=0.1 Roff=100000 | zener Vz=5 Vf=0 Ron=0.1 Roff=100000 | Vs=0.1
    coil   L=2.3957829546 Rcoil=1 | grid nper=50 dt=0.1 nsteps=150
    schedule of one period: closed on steps 1..40, open on steps 41..50
    max number of fixed points found in any step: 1  (g(v3) strictly decreasing at every step)

    --- the four ASSERTED quantities ---
    step 140 t= 14.00 closed END OF AN ON PHASE
        i_coil =8.0812231477e+00  i_switch=+8.081313e+00  i_clamp=-9.0004208071e-05
        v_coil =+9.1918686848  v_diode=-9.000421 (R=100000)  v_zener=+0.191448 (R=2127.1)
    step 150 t= 15.00 open   END OF AN OFF PHASE
        i_coil =3.0415318831e+00  i_switch=+1.687839e-04  i_clamp=+3.0413630992e+00
        v_coil =-6.8783891882  v_diode=+1.329985 (R=0.437299)  v_zener=-5.548404 (R=1.82431)

    --- clamp check at step 150 (the OFF-phase node voltage) ---
      ideal clamp -(Vz + Vf + I*(Ron_d+Ron_z)) = -6.308306 V
      solved v_coil = v(A) - v(G)              = -6.878389 V
      knee-width overshoot                     = 0.570083 V = 5.7*Vs

    --- blocking check at step 140 (the ON-phase clamp leakage) ---
      node A at +9.191869 V, clamp branch leaks -9.000421e-05 A
      |v_A|/Roff       = 9.191869e-05 A   ratio 0.9792
      |v_A|/(2*Roff)   = 4.595934e-05 A   ratio 1.9583
      the SERIES DIODE blocks alone: R_diode=100000 ohm at -9.0004 V, while
      the zener sits on its FORWARD knee: R_zener=2127.1 ohm at +0.1914 V

    --- trajectory across the opening edge of period 3 ---
      n=139 t= 13.90 closed i_coil=+8.0348647940e+00 i_sw=+8.0350e+00 i_clamp=-9.0050e-05 v_coil= +9.19650
      n=140 t= 14.00 closed i_coil=+8.0812231477e+00 i_sw=+8.0813e+00 i_clamp=-9.0004e-05 v_coil= +9.19187
      n=141 t= 14.10 open   i_coil=+7.4771527129e+00 i_sw=+1.6995e-04 i_clamp=+7.4770e+00 v_coil= -6.99506
      n=142 t= 14.20 open   i_coil=+6.8977799916e+00 i_sw=+1.6983e-04 i_clamp=+6.8976e+00 v_coil= -6.98273
      n=145 t= 14.50 open   i_coil=+5.2981315601e+00 i_sw=+1.6946e-04 i_clamp=+5.2980e+00 v_coil= -6.94591
      n=149 t= 14.90 open   i_coil=+3.4555894692e+00 i_sw=+1.6893e-04 i_clamp=+3.4554e+00 v_coil= -6.89314
      n=150 t= 15.00 open   i_coil=+3.0415318831e+00 i_sw=+1.6878e-04 i_clamp=+3.0414e+00 v_coil= -6.87839

    --- periodicity (same phase in successive periods) ---
      n= 40 t=  4.00 closed i_coil=7.5816280142e+00
      n= 90 t=  9.00 closed i_coil=8.0320915159e+00
      n=140 t= 14.00 closed i_coil=8.0812231477e+00
      n= 50 t=  5.00 open   i_coil=2.7133057680e+00
      n=100 t= 10.00 open   i_coil=3.0092434842e+00
      n=150 t= 15.00 open   i_coil=3.0415318831e+00

    --- EVERY ELEMENT IS LOAD BEARING: the same run with one of them changed ---
      Each variant is its own trajectory from t = 0, so the comparison is the
      FRACTION of the coil current the OFF window removes, not the raw current.
      step   this test         breakdown removed    single-diode freewheel
             i_coil  v_coil    i_coil    v_coil     i_coil    v_coil
      140  8.081223   +9.1919  7.581645   +9.2418  8.424645   +9.1575
      141  7.477153   -6.9951  4.819290  -61.3607  8.030788   -1.4052
      145  5.298132   -6.9459  0.000100   +0.0001  6.608838   -1.3861
      150  3.041532   -6.8784  0.000100   +0.0001  5.132850   -1.3654
      this test          : OFF window 140 -> 150 removes  62.4% of the coil current
      breakdown removed  : OFF window 140 -> 150 removes 100.0% of the coil current
      single-diode freewh: OFF window 140 -> 150 removes  39.1% of the coil current
      BREAKDOWN REMOVED: the two elements are then plain diodes in ANTI-SERIES, so the
      branch cannot conduct in EITHER direction. The coil is left with nothing but
      the two 1e5 ohm off-state paths, node A is thrown to -61.4 V and the current
      is gone by step 145 -- exactly the freewheel-less behaviour of the sibling
      circuits2D_transient_switch_pwm. The zener's SECOND KNEE is the only thing
      that opens this path, and the series diode is the only thing that closes it
      again during the ON phase.
      SINGLE-DIODE FREEWHEEL (the sibling topology, same 0.7 V part): the coil sees
      only -1.3654 V instead of -6.8784 V during the OFF phase, so it demagnetizes by
      39.1% instead of 62.4% in the same 1 s window. THAT is the accelerated
      demagnetization the zener clamp exists for.
      SERIES DIODE DELETED (zener straight across the coil): its Vf = 0 forward knee
      then SHORTS the coil during the ON phase -- node A collapses from +9.1919 V
      to +4.7973 V (the Ron_s/Ron_z divider) while the switch carries 52.0 A of
      shoot-through into a clamp branch drawing 48.0 A. The coil, no longer driven,
      reaches only 4.054777 A against 8.081223 A at step 140. The series diode is what
      stops the clamp branch from eating the ON phase; the OFF phase clamps at
      -5.4826 V either way.

WHAT ELMER PRODUCES (build of 2026-08-18, w = 0.3, Newton on, residual gate 1e-8),
against this reference, which shares no code with it:
    i_coil  step 140 = +8.0812231476610e+00 A   rel. deviation 4.8e-12
    i_clamp step 140 = -9.0004208068550e-05 A   rel. deviation 2.7e-11
    i_coil  step 150 = +3.0415318829040e+00 A   rel. deviation 6.4e-11
    v_coil  step 150 = -6.8783891887220e+00 V   rel. deviation 7.6e-11
and every row of the circuit is satisfied to machine precision in the solved
dofs (KCL at node M 0.0e+00, KCL at node A 1.1e-13, KVL of the clamp loop
0.0e+00, KVL of the source loop 2.0e-12).
"""
import math

# ---------------------------------------------------------------- parameters
V0     = 10.0                 # DC source [V]
fsw    = 0.2                  # switching frequency [Hz]
duty   = 0.8
phase  = 3.6                  # [deg] = 0.01 of a period
Ron_s  = 0.1                  # switch on  resistance
Roff_s = 1.0e5                # switch off resistance

Vf_d   = 0.7                  # series diode: forward knee [V]
Ron_d  = 0.1
Roff_d = 1.0e5

Vz     = 5.0                  # zener: breakdown voltage [V]
Vf_z   = 0.0                  # zener: forward knee [V]
Ron_z  = 0.1
Roff_z = 1.0e5

Vs     = 0.1                  # smoothing voltage, both elements, both knees
L      = 2.3957829546         # calibrated, see docstring
Rcoil  = 1.0
nper   = 50
dt     = (1.0 / fsw) / nper
a      = L / dt
NSTEPS = 150


def closed(t):
    """Gate state at time t, exactly as the Fortran evaluates it."""
    return math.fmod(t * fsw - phase / 360.0, 1.0) % 1.0 < duty


def Rs(t):
    return Ron_s if closed(t) else Roff_s


def Rd(v):
    """Series diode: the single-knee law of the plain Diode component."""
    return Ron_d + 0.5 * (Roff_d - Ron_d) * (1.0 - math.tanh((v - Vf_d) / Vs))


def Rz(v, breakdown=True):
    """Zener: the two-knee law.  breakdown=False is the same element with
    'Diode Breakdown Voltage' left out, i.e. the plain single-knee law with
    Vf = Vf_z; the clamp branch is then two diodes in ANTI-SERIES and cannot
    conduct at all, which is one of the variants computed at the bottom."""
    if breakdown:
        return Ron_z + 0.25 * (Roff_z - Ron_z) \
            * (1.0 - math.tanh((v - Vf_z) / Vs)) \
            * (1.0 + math.tanh((v + Vz) / Vs))
    return Ron_z + 0.5 * (Roff_z - Ron_z) * (1.0 - math.tanh((v - Vf_z) / Vs))


# ------------------------------------------------------------------ helpers
# i(v) = v/R(v) is STRICTLY INCREASING for both laws: di/dv = (R - v*R')/R**2,
# and R - v*R' > 0 everywhere (R > 0, and v*R' <= 0 on every branch of either
# law).  The element characteristic is therefore invertible, which is what lets
# the clamp branch -- two elements sharing one current -- be reduced to a
# single scalar equation below.  Nothing is taken on trust: the consequence of
# it -- that the outer residual g(v3) is strictly decreasing, hence has exactly
# one root -- is re-checked numerically on the scan grid at every timestep, and
# step_bdf1 raises if it ever fails.
BRACKET = 60.0                # |v| bracket for the element inverses [V]


def inv_z(itarget, breakdown=True):
    """Zener branch voltage that carries the current itarget (bisection on the
    monotone characteristic i = v/Rz(v))."""
    lo, hi = -BRACKET, BRACKET
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if m / Rz(m, breakdown) < itarget:
            lo = m
        else:
            hi = m
        if hi - lo < 1e-14 * max(1.0, abs(m)):
            break
    return 0.5 * (lo + hi)


def zener_v(ic, mode):
    """Branch voltage of the SECOND clamp element, given that the clamp branch
    carries ic in the G -> M -> A direction (so the element, which points
    A -> M, carries -ic).

      'zener'     the test itself: the two-knee law, reverse breakdown
      'nobreak'   the same element with 'Diode Breakdown Voltage' removed
      'freewheel' the element deleted altogether, i.e. the single-diode clamp
                  of the sibling circuits2D_transient_switch_diode_freewheel
    """
    if mode == 'freewheel':
        return 0.0
    return inv_z(-ic, mode == 'zener')


def step_bdf1(t, i1prev, mode='zener', scan=False):
    """One implicit BDF1 step of the zener-clamped buck cell.

    Unknown: v3, the SERIES DIODE's own branch voltage.  Everything else
    follows explicitly (only the zener inverse needs a bisection):

        i_clamp = v3 / Rd(v3)                    (diode characteristic)
        v4      = inv_z(-i_clamp)                (zener carries -i_clamp)
        v1      = v4 - v3                        (KVL of the clamp loop)
        i_coil  = (a*i_coil_prev + v1)/(a + Rcoil)      a = L/dt
        i_sw    = (V0 - v1)/Rs(t)                (KVL of the source loop)
        g(v3)   = i_coil - i_sw - i_clamp        (KCL at node A)

    g is strictly decreasing in v3, so the step has exactly one solution.
    Returns (i_coil, i_switch, i_clamp, v1, v3, v4, nroots).
    """
    R_s = Rs(t)

    def parts(v3):
        ic = v3 / Rd(v3)
        v4 = zener_v(ic, mode)
        v1 = v4 - v3
        i1 = (a * i1prev + v1) / (a + Rcoil)
        i2 = (V0 - v1) / R_s
        return i1, i2, ic, v1, v4

    def g(v3):
        i1, i2, ic, v1, v4 = parts(v3)
        return i1 - i2 - ic

    nroots = 1
    if scan:
        nroots, prev_g, decreasing = 0, None, True
        for v in GRID:
            gv = g(v)
            if prev_g is not None:
                if prev_g * gv < 0.0 or gv == 0.0:
                    nroots += 1
                if gv > prev_g:
                    decreasing = False
            prev_g = gv
        if not decreasing:
            raise RuntimeError('g(v3) is not monotone at t=%g' % t)

    lo, hi = -BRACKET, BRACKET
    glo = g(lo)
    for _ in range(300):
        m = 0.5 * (lo + hi)
        gm = g(m)
        if glo * gm <= 0.0:
            hi = m
        else:
            lo, glo = m, gm
        if hi - lo < 1e-14 * max(1.0, abs(m)):
            break
    v3 = 0.5 * (lo + hi)
    i1, i2, ic, v1, v4 = parts(v3)
    return i1, i2, ic, v1, v3, v4, nroots


def _grid():
    """Coarse scan over the whole plausible range of the diode voltage plus a
    fine scan across its knee (width ~Vs around Vf_d) and across zero, where
    the blocking-phase operating point sits."""
    g = [-BRACKET + 2.0 * BRACKET * k / 600.0 for k in range(601)]
    for c in (Vf_d, 0.0):
        g += [c - 2.0 + 4.0 * k / 800.0 for k in range(801)]
    return sorted(set(g))

GRID = _grid()


def step_nodiode(t, i1prev):
    """One BDF1 step of the SAME cell with the SERIES DIODE deleted, i.e. the
    zener sitting straight across the coil.  The unknown is then the zener's
    own voltage v4 (= v_coil, since node M collapses onto the return rail):

        i_clamp = -v4 / Rz(v4)                   (clamp current G -> A)
        i_coil  = (a*i1prev + v4)/(a + Rcoil)
        i_sw    = (V0 - v4)/Rs(t)
        g(v4)   = i_coil - i_sw - i_clamp        (KCL at node A)

    g is strictly increasing in v4, so again exactly one solution.
    """
    R_s = Rs(t)

    def g(v4):
        return (a * i1prev + v4) / (a + Rcoil) - (V0 - v4) / R_s + v4 / Rz(v4)

    lo, hi = -BRACKET, BRACKET
    for _ in range(300):
        m = 0.5 * (lo + hi)
        if g(m) < 0.0:
            lo = m
        else:
            hi = m
        if hi - lo < 1e-14 * max(1.0, abs(m)):
            break
    v4 = 0.5 * (lo + hi)
    return ((a * i1prev + v4) / (a + Rcoil), (V0 - v4) / R_s,
            -v4 / Rz(v4), v4)


def run_nodiode(nsteps=NSTEPS):
    i1, out = 0.0, []
    for n in range(1, nsteps + 1):
        t = n * dt
        i1, i2, ic, v4 = step_nodiode(t, i1)
        out.append((n, t, closed(t), i1, i2, ic, v4, 0.0, v4))
    return out


def run_bdf1(nsteps=NSTEPS, mode='zener', scan=False):
    i1, out, maxroots = 0.0, [], 0
    for n in range(1, nsteps + 1):
        t = n * dt
        i1, i2, ic, v1, v3, v4, nr = step_bdf1(t, i1, mode, scan)
        maxroots = max(maxroots, nr)
        out.append((n, t, closed(t), i1, i2, ic, v1, v3, v4))
    return out, maxroots


if __name__ == '__main__':
    print('source V0=%g | switch fsw=%g duty=%g phase=%gdeg Ron=%g Roff=%g'
          % (V0, fsw, duty, phase, Ron_s, Roff_s))
    print('diode  Vf=%g Ron=%g Roff=%g | zener Vz=%g Vf=%g Ron=%g Roff=%g | '
          'Vs=%g' % (Vf_d, Ron_d, Roff_d, Vz, Vf_z, Ron_z, Roff_z, Vs))
    print('coil   L=%.10f Rcoil=%g | grid nper=%d dt=%g nsteps=%d'
          % (L, Rcoil, nper, dt, NSTEPS))

    sol, maxroots = run_bdf1(scan=True)
    on = [n for n, t, c, i1, i2, ic, v1, v3, v4 in sol if c and n <= nper]
    print('schedule of one period: closed on steps %d..%d, open on steps %d..%d'
          % (on[0], on[-1], on[-1] + 1, nper))
    print('max number of fixed points found in any step: %d  '
          '(g(v3) strictly decreasing at every step)' % maxroots)

    print('\n--- the four ASSERTED quantities ---')
    for want, tag in ((140, 'END OF AN ON PHASE'), (150, 'END OF AN OFF PHASE')):
        n, t, c, i1, i2, ic, v1, v3, v4 = sol[want - 1]
        print('step %3d t=%6.2f %-6s %s' % (n, t, 'closed' if c else 'open', tag))
        print('    i_coil =%.10e  i_switch=%+.6e  i_clamp=%+.10e'
              % (i1, i2, ic))
        print('    v_coil =%+.10f  v_diode=%+.6f (R=%.6g)  '
              'v_zener=%+.6f (R=%.6g)' % (v1, v3, Rd(v3), v4, Rz(v4)))

    n, t, c, i1, i2, ic, v1, v3, v4 = sol[149]
    ideal = -(Vz + Vf_d + i1 * (Ron_d + Ron_z))
    print('\n--- clamp check at step 150 (the OFF-phase node voltage) ---')
    print('  ideal clamp -(Vz + Vf + I*(Ron_d+Ron_z)) = %+.6f V' % ideal)
    print('  solved v_coil = v(A) - v(G)              = %+.6f V' % v1)
    print('  knee-width overshoot                     = %.6f V = %.1f*Vs'
          % (abs(v1 - ideal), abs(v1 - ideal) / Vs))

    n, t, c, i1, i2, ic, v1, v3, v4 = sol[139]
    print('\n--- blocking check at step 140 (the ON-phase clamp leakage) ---')
    print('  node A at %+.6f V, clamp branch leaks %+.6e A' % (v1, ic))
    print('  |v_A|/Roff       = %.6e A   ratio %.4f'
          % (abs(v1) / Roff_d, abs(ic) / (abs(v1) / Roff_d)))
    print('  |v_A|/(2*Roff)   = %.6e A   ratio %.4f'
          % (abs(v1) / (2 * Roff_d), abs(ic) / (abs(v1) / (2 * Roff_d))))
    print('  the SERIES DIODE blocks alone: R_diode=%.6g ohm at %+.4f V, while'
          % (Rd(v3), v3))
    print('  the zener sits on its FORWARD knee: R_zener=%.6g ohm at %+.4f V'
          % (Rz(v4), v4))

    print('\n--- trajectory across the opening edge of period 3 ---')
    for n, t, c, i1, i2, ic, v1, v3, v4 in sol:
        if n in (139, 140, 141, 142, 145, 149, 150):
            print('  n=%3d t=%6.2f %-6s i_coil=%+.10e i_sw=%+.4e '
                  'i_clamp=%+.4e v_coil=%+9.5f'
                  % (n, t, 'closed' if c else 'open', i1, i2, ic, v1))

    print('\n--- periodicity (same phase in successive periods) ---')
    for want in (40, 90, 140, 50, 100, 150):
        n, t, c, i1 = sol[want - 1][0], sol[want - 1][1], sol[want - 1][2], \
            sol[want - 1][3]
        print('  n=%3d t=%6.2f %-6s i_coil=%.10e'
              % (n, t, 'closed' if c else 'open', i1))

    print('\n--- EVERY ELEMENT IS LOAD BEARING: the same run with one of them '
          'changed ---')
    nb, _ = run_bdf1(mode='nobreak')
    fw, _ = run_bdf1(mode='freewheel')
    print('  Each variant is its own trajectory from t = 0, so the comparison '
          'is the')
    print('  FRACTION of the coil current the OFF window removes, not the raw '
          'current.')
    print('  step   this test         breakdown removed    single-diode '
          'freewheel')
    print('         i_coil  v_coil    i_coil    v_coil     i_coil    v_coil')
    for want in (140, 141, 145, 150):
        _, _, _, i1c, _, _, v1c, _, _ = sol[want - 1]
        _, _, _, i1n, _, _, v1n, _, _ = nb[want - 1]
        _, _, _, i1f, _, _, v1f, _, _ = fw[want - 1]
        print('  %3d  %8.6f %+9.4f  %8.6f %+9.4f  %8.6f %+9.4f'
              % (want, i1c, v1c, i1n, v1n, i1f, v1f))
    for tag, s in (('this test          ', sol), ('breakdown removed  ', nb),
                   ('single-diode freewh', fw)):
        print('  %s: OFF window 140 -> 150 removes %5.1f%% of the coil current'
              % (tag, 100.0 * (s[139][3] - s[149][3]) / s[139][3]))
    print('  BREAKDOWN REMOVED: the two elements are then plain diodes in ANTI-'
          'SERIES, so the')
    print('  branch cannot conduct in EITHER direction. The coil is left with '
          'nothing but')
    print('  the two 1e5 ohm off-state paths, node A is thrown to %.1f V and '
          'the current' % nb[140][6])
    print('  is gone by step 145 -- exactly the freewheel-less behaviour of the'
          ' sibling')
    print('  circuits2D_transient_switch_pwm. The zener\'s SECOND KNEE is the '
          'only thing')
    print('  that opens this path, and the series diode is the only thing that '
          'closes it')
    print('  again during the ON phase.')
    print('  SINGLE-DIODE FREEWHEEL (the sibling topology, same 0.7 V part): '
          'the coil sees')
    print('  only %+.4f V instead of %+.4f V during the OFF phase, so it '
          'demagnetizes by' % (fw[149][6], sol[149][6]))
    print('  %.1f%% instead of %.1f%% in the same 1 s window. THAT is the '
          'accelerated' % (100.0 * (fw[139][3] - fw[149][3]) / fw[139][3],
                           100.0 * (sol[139][3] - sol[149][3]) / sol[139][3]))
    print('  demagnetization the zener clamp exists for.')

    nd = run_nodiode()
    print('  SERIES DIODE DELETED (zener straight across the coil): its Vf = 0 '
          'forward knee')
    print('  then SHORTS the coil during the ON phase -- node A collapses from '
          '%+.4f V' % sol[139][6])
    print('  to %+.4f V (the Ron_s/Ron_z divider) while the switch carries '
          '%.1f A of' % (nd[139][6], nd[139][4]))
    print('  shoot-through into a clamp branch drawing %.1f A. The coil, no '
          'longer driven,' % abs(nd[139][5]))
    print('  reaches only %.6f A against %.6f A at step 140. The series diode '
          'is what' % (nd[139][3], sol[139][3]))
    print('  stops the clamp branch from eating the ON phase; the OFF phase '
          'clamps at')
    print('  %+.4f V either way.' % nd[149][6])
