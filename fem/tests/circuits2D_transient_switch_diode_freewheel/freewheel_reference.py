"""Independent lumped reference for circuits2D_transient_switch_diode_freewheel.

Run with any Python 3 (no third-party packages):  python freewheel_reference.py

Circuit (identical to the sif) -- the buck / chopper cell:

           +--- switch ---+---------+
           |    S       A |         |
  V0 (DC)--+              |        L, Rcoil
           |              |         |
           +--------------+         |
           |        diode |         |
           +--------------+---------+
                   G

    switch:  R_s(t) = Ron_s if MODULO(t*fsw - phase/360, 1) < duty else Roff_s
    diode :  R_d(V) = Ron + 0.5*(Roff-Ron)*(1 - tanh((V - Vf)/Vs)),  V = -v_coil
    coil  :  L*(i1_n - i1_{n-1})/dt + Rcoil*i1_n = v_coil = -V

One BDF1 step in the single unknown V = v_diode (all three branch currents
follow from it):

    i_coil   = (a*i1_prev - V) / (a + Rcoil),      a = L/dt
    i_diode  = V / R_d(V)
    i_switch = (V0 + V) / R_s(t)                   (v_switch = V0 - v_coil = V0 + V)
    KCL at node A:  i_switch + i_diode - i_coil = 0

so the residual is  g(V) = (V0 + V)/R_s(t) - (i_coil - i_diode), solved by
scanning for sign changes and bisecting -- never by Picard, and with no code in
common with Elmer.  The scan doubles as a diagnostic: it reports how many fixed
points each timestep has.

R_s is evaluated at t_n, the END of the step, which is the BDF1 implicit point
and exactly what the Fortran does (it reads GetTime(), the current solution
time).  With fsw = 0.2 Hz, dt = Tsw/50 and phase = 3.6 deg = 0.01 of a period
both switching edges sit exactly halfway between grid points, so the closed/open
decision has half a step of margin in either code.  Schedule: steps
n mod 50 = 1..20 closed, 21..50 open.

L is not a free parameter.  It was calibrated from an Elmer run of this exact
mesh and coil with a fixed 1 ohm resistor in place of the nonlinear element,
250 steps: per-step solution of the BDF1 recurrence gives min 2.395782954 H /
max 2.395782960 H / median 2.3957829546 H, and a golden-section fit minimising
the worst replay error over the whole run lands on the same 2.3957829546 H with
a replay error of 5.95e-13 A on a 3.172770 A amplitude.

RECORDED OUTPUT (2026-08-17, CPython 3.11.9):

    params: V0=10 fsw=0.2 duty=0.4 phase=3.6deg Ron_s=0.1 Roff_s=100000 | Vf=0 Ron=0.1 Roff=100000 Vs=0.1 | L=2.3957829546 Rcoil=1 nper=50 dt=0.1
    schedule of one period: closed on steps 1..20, open on steps 21..50
    max number of fixed points found in any step: 1

    --- the two ASSERTED instants ---
    step 120  t= 12.000  closed i_coil=5.8933886201e+00  i_switch=+5.893483e+00  i_diode=-9.410652e-05  v_diode=-9.410652  R_diode=100000
    step 150  t= 15.000  open   i_coil=1.2581799855e+00  i_switch=+1.062215e-04  i_diode=+1.258074e+00  v_diode=+0.622150  R_diode=0.494526

    --- trajectory around the switching edges of period 3 ---
      n=119 t= 11.900 closed i_coil=+5.7465780295e+00 i_sw=+5.7467e+00 i_d=-9.4253e-05 v_d= -9.42533
      n=120 t= 12.000 closed i_coil=+5.8933886201e+00 i_sw=+5.8935e+00 i_d=-9.4107e-05 v_d= -9.41065
      n=121 t= 12.100 open   i_coil=+5.6273371036e+00 i_sw=+1.0747e-04 i_d=+5.6272e+00 v_d= +0.74668
      n=122 t= 12.200 open   i_coil=+5.3722297537e+00 i_sw=+1.0740e-04 i_d=+5.3721e+00 v_d= +0.73959
      n=149 t= 14.900 open   i_coil=+1.3366649576e+00 i_sw=+1.0626e-04 i_d=+1.3366e+00 v_d= +0.62562
      n=150 t= 15.000 open   i_coil=+1.2581799855e+00 i_sw=+1.0622e-04 i_d=+1.2581e+00 v_d= +0.62215

    --- periodicity (same phase in successive periods) ---
      n= 20 t=  2.000 closed i_coil=5.3867611447e+00
      n= 70 t=  7.000 closed i_coil=5.8410884254e+00
      n=120 t= 12.000 closed i_coil=5.8933886201e+00
      n= 50 t=  5.000 open   i_coil=1.1150346928e+00
      n=100 t= 10.000 open   i_coil=1.2433926718e+00
      n=150 t= 15.000 open   i_coil=1.2581799855e+00

    --- freewheel check: is the coil current really carried by the diode when the switch is open? ---
      n=121  i_coil=5.627337  i_diode=5.627230  i_switch=1.075e-04  i_diode/i_coil=0.999981
      n=135  i_coil=2.868638  i_diode=2.868531  i_switch=1.068e-04  i_diode/i_coil=0.999963
      n=150  i_coil=1.258180  i_diode=1.258074  i_switch=1.062e-04  i_diode/i_coil=0.999916

Elmer reproduces this whole 150-step trajectory to 6.7e-09 A on a coil current
that ranges over 1.1..5.9 A; the two asserted steps come out at relative
deviations of 8.8e-11 (step 120) and 1.6e-09 (step 150).
"""
import math

# ---------------------------------------------------------------- parameters
V0     = 10.0                 # DC source [V]
fsw    = 0.2                  # switching frequency [Hz]
duty   = 0.4
phase  = 3.6                  # [deg]
Ron_s  = 0.1                  # switch on  resistance
Roff_s = 1.0e5                # switch off resistance
Vf     = 0.0                  # diode forward voltage
Ron    = 0.1                  # diode on  resistance
Roff   = 1.0e5                # diode off resistance
Vs     = 0.1                  # diode smoothing voltage
L      = 2.3957829546         # calibrated, see docstring
Rcoil  = 1.0
nper   = 50
dt     = (1.0 / fsw) / nper
a      = L / dt


def closed(t):
    """Gate state at time t, exactly as the Fortran evaluates it."""
    return math.fmod(t * fsw - phase / 360.0, 1.0) % 1.0 < duty


def Rs(t):
    return Ron_s if closed(t) else Roff_s


def Rd(v):
    return Ron + 0.5 * (Roff - Ron) * (1.0 - math.tanh((v - Vf) / Vs))


def _grid():
    """Coarse scan over the whole plausible range of the diode voltage plus a
    fine scan across the diode transition (width ~Vs around Vf)."""
    g = [(-3.0 * V0) + (6.0 * V0) * k / 12000.0 for k in range(12001)]
    g += [Vf - 2.0 + 4.0 * k / 40000.0 for k in range(40001)]
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


def step_bdf1(t, i1prev, scan=True):
    """One implicit BDF1 step.  Returns (i_coil, i_switch, i_diode, V, nroots)."""
    R_s = Rs(t)

    def parts(v):
        i1 = (a * i1prev - v) / (a + Rcoil)
        i3 = v / Rd(v)
        i2 = (V0 + v) / R_s
        return i1, i2, i3

    def g(v):
        i1, i2, i3 = parts(v)
        return i2 + i3 - i1

    if not scan:
        lo, hi = -3.0 * V0, 3.0 * V0
        v = _bisect(g, lo, hi, g(lo))
        i1, i2, i3 = parts(v)
        return i1, i2, i3, v, 1

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
        raise RuntimeError('no root at t=%g' % t)
    i1, i2, i3 = parts(roots[0])
    return i1, i2, i3, roots[0], len(roots)


def run_bdf1(nsteps, scan=True):
    i1, out, maxroots = 0.0, [], 1
    for n in range(1, nsteps + 1):
        t = n * dt
        i1, i2, i3, v, nr = step_bdf1(t, i1, scan=scan)
        maxroots = max(maxroots, nr)
        out.append((n, t, closed(t), i1, i2, i3, v))
    return out, maxroots


if __name__ == '__main__':
    print('params: V0=%g fsw=%g duty=%g phase=%gdeg Ron_s=%g Roff_s=%g | '
          'Vf=%g Ron=%g Roff=%g Vs=%g | L=%.10f Rcoil=%g nper=%d dt=%g'
          % (V0, fsw, duty, phase, Ron_s, Roff_s, Vf, Ron, Roff, Vs, L, Rcoil,
             nper, dt))
    sol, maxroots = run_bdf1(150)
    on = [n for n, t, c, i1, i2, i3, v in sol if c and n <= nper]
    print('schedule of one period: closed on steps %d..%d, open on steps %d..%d'
          % (on[0], on[-1], on[-1] + 1, nper))
    print('max number of fixed points found in any step: %d' % maxroots)

    print('\n--- the two ASSERTED instants ---')
    for want in (120, 150):
        n, t, c, i1, i2, i3, v = sol[want - 1]
        print('step %3d  t=%7.3f  %-6s i_coil=%.10e  i_switch=%+.6e  '
              'i_diode=%+.6e  v_diode=%+.6f  R_diode=%.6g'
              % (n, t, 'closed' if c else 'open', i1, i2, i3, v, Rd(v)))

    print('\n--- trajectory around the switching edges of period 3 ---')
    for n, t, c, i1, i2, i3, v in sol:
        if n in (119, 120, 121, 122, 149, 150):
            print('  n=%3d t=%7.3f %-6s i_coil=%+.10e i_sw=%+.4e i_d=%+.4e '
                  'v_d=%+9.5f' % (n, t, 'closed' if c else 'open', i1, i2, i3, v))

    print('\n--- periodicity (same phase in successive periods) ---')
    for want in (20, 70, 120, 50, 100, 150):
        n, t, c, i1 = sol[want - 1][0], sol[want - 1][1], sol[want - 1][2], sol[want - 1][3]
        print('  n=%3d t=%7.3f %-6s i_coil=%.10e' % (n, t, 'closed' if c else 'open', i1))

    print('\n--- freewheel check: is the coil current really carried by the '
          'diode when the switch is open? ---')
    for want in (121, 135, 150):
        n, t, c, i1, i2, i3, v = sol[want - 1]
        print('  n=%3d  i_coil=%.6f  i_diode=%.6f  i_switch=%.3e  '
              'i_diode/i_coil=%.6f' % (n, i1, i3, i2, i3 / i1))
