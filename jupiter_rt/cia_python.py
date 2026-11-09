"""H2-H2, H2-He and H2-CH4 collision-induced absorption."""
from functools import lru_cache
from pathlib import Path
import math
import numpy as np
PI = 3.1415926535898
CLOSCHM = 2.68675484e+19
HBAR = 1.054588757e-27
CLIGHT = 29979250000.0
BOLTZWN = 0.6950304
CALIB_BASE = 2 * PI * CLIGHT * (4 * PI ** 2 / (3 * HBAR * CLIGHT)) * CLOSCHM ** 2

def _poly(x, coeffs):
    r = np.zeros_like(x, dtype=np.float64) if isinstance(x, np.ndarray) else 0.0
    for v in coeffs:
        r = r * x + v
    return r
P1 = (0.0045813, 0.0360768, 0.2659732, 1.2067492, 3.0899424, 3.5156229, 1.0)
P2 = (7.4e-06, 0.0001075, 0.00262698, 0.0348859, 0.23069756, 0.4227842, -0.57721566)
P3 = (0.00032411, 0.00301532, 0.02658733, 0.15084934, 0.51498869, 0.87890594, 0.5)
P4 = (-4.686e-05, -0.00110404, -0.01919402, -0.18156897, -0.67278579, 0.15443144, 1.0)
P5 = (0.00053208, -0.0025154, 0.00587872, -0.01062446, 0.02189568, -0.07832358, 1.25331414)
P6 = (-0.00068245, 0.00325614, -0.00780353, 0.01504268, -0.0365562, 0.23498619, 1.25331414)

def _xk0(x):
    x = np.asarray(x, dtype=np.float64)
    small = x <= 2.0
    out = np.empty_like(x)
    z = x[small]
    if z.size:
        with np.errstate(divide='ignore'):
            out[small] = -np.log(z / 2) * _poly((z / 3.75) ** 2, P1) + _poly((z / 2) ** 2, P2)
    z = x[~small]
    if z.size:
        out[~small] = np.exp(-np.minimum(z, 330.0)) * _poly(2 / z, P5) / np.sqrt(np.minimum(z, 330.0))
    return out

def _xk1(x):
    x = np.asarray(x, dtype=np.float64)
    small = x <= 2.0
    out = np.empty_like(x)
    z = x[small]
    if z.size:
        with np.errstate(divide='ignore', invalid='ignore'):
            out[small] = z * np.log(z / 2) * (z * _poly((z / 3.75) ** 2, P3)) + _poly((z / 2) ** 2, P4)
        out[small] = np.where(z == 0.0, 1.0, out[small])
    z = x[~small]
    if z.size:
        zz = np.minimum(z, 330.0)
        out[~small] = np.sqrt(zz) * np.exp(-zz) * _poly(2 / z, P6)
    return out

def _bgama(fnu, t1, t2, eps, t3, t4, temp, *, methane=False):
    clight = 29979245800.0 if methane else 29979245800.0
    omega = 2 * PI * clight * np.asarray(fnu, dtype=np.float64)
    t0 = 1.0545887e-27 / (2 * 1.380662e-16 * temp)
    z = np.sqrt((1 + (omega * t1) ** 2) * (t2 * t2 + t0 * t0)) / t1
    m = z <= 2
    k1 = np.empty_like(z)
    v = z[m]
    if v.size:
        with np.errstate(divide='ignore'):
            k1[m] = v * v * np.log(v / 2) * _poly((v / 3.75) ** 2, P3) + _poly((v / 2) ** 2, P4)
    v = z[~m]
    if v.size:
        k1[~m] = np.sqrt(v) * np.exp(-v) * _poly(2 / v, P6)
    first = t1 / PI * np.exp(t2 / t1 + t0 * omega) * k1 / (1 + (t1 * omega) ** 2)
    if eps == 0.0:
        return first
    zp = np.sqrt((1 + (omega * t4) ** 2) * (t3 * t3 + t0 * t0)) / t4
    k0 = np.empty_like(zp)
    m = zp <= 2
    v = zp[m]
    if v.size:
        k0[m] = -np.log(v / 2) * _poly((v / 3.75) ** 2, P1) + _poly((v / 2) ** 2, P2)
    v = zp[~m]
    if v.size:
        k0[~m] = np.exp(-v) * _poly(2 / v, P5) / np.sqrt(v)
    return (first + eps * (t3 / PI) * np.exp(t3 / t4 + t0 * omega) * k0) / (1 + eps)

def _bgama0(fnu, t5, t6, temp):
    tau4 = np.sqrt(t5 * t5 + (7.638280918e-12 / (temp * 2)) ** 2)
    omega = 188365156800.0 * np.asarray(fnu, dtype=np.float64)
    x = tau4 * np.sqrt(1 / (t6 * t6) + omega * omega)
    return t5 / 3.141592654 * np.exp(min(t5 / t6, 430.0) + np.asarray(fnu) / (2 * 0.6950304256 * temp)) * _xk0(x)

def _bgama1(fnu, t1, t2, temp):
    tau3 = np.sqrt(t2 * t2 + (7.638280918e-12 / (temp * 2)) ** 2)
    omega = 188365156800.0 * np.asarray(fnu, dtype=np.float64)
    denom = 1 + (omega * t1) ** 2
    x = tau3 / t1 * np.sqrt(denom)
    return t1 / 3.141592654 * np.exp(min(t2 / t1, 430.0) + np.asarray(fnu) / (2 * 0.6950304256 * temp)) * _xk1(x) / denom

def _eh2_vib(n, ii):
    h = n + 0.5
    return 4395.34 * h - 117.905 * h * h + (60.809 - 2.993 * h + 0.025 * h * h) * ii - (0.04648 - 0.00134 * h) * ii * ii

def _ph2_vib(j, temp):
    return (2 * j + 1) * (1 if j % 2 == 0 else 3) * np.exp(-1.4387859 / temp * _eh2_vib(0, j * (j + 1)))

def _parts_vib(temp):
    q = 0.0
    j = 0
    while True:
        dq = _ph2_vib(j, temp)
        q += dq
        j += 1
        if dq <= q / 900.0:
            return (q, j)
        if j > 500:
            raise RuntimeError('The H2 partition-function sum failed to converge')

def _fctl(n):
    if n <= 1:
        return 1.0
    if n <= 15:
        return float(math.factorial(n))
    z = n + 1.0
    p = (((-0.0002294720936 / z - 0.00268132716) / z + 0.003472222222) / z + 0.08333333333) / z + 1.0
    return np.exp(-z) * z ** (z - 0.5) * p * 2.506628274631

@lru_cache(maxsize=4096)
def _clebsqr(l, lam, lp):
    if l + lam < lp or lam + lp < l or l + lp < lam or (l + lp + lam) % 2 or (min(l, lam, lp) < 0):
        return 0.0
    f = 1.0 / (l + lp + 1 - lam)
    if lam:
        for i in range((l + lp - lam) // 2 + 1, (l + lp + lam) // 2 + 1):
            f = f * i / (2 * (2 * i + 1))
    p = (2 * lp + 1) * f * _fctl(lam + l - lp) * _fctl(lam + lp - l)
    return p / (_fctl((lam + l - lp) // 2) * _fctl((lam + lp - l) // 2)) ** 2

def _addspec_h2(freq, temp, q, jrange, s, e, t1, t2, t3, t4, lam1, lam2, lam, nv1=0, nv2=0, like=1):
    calib = CALIB_BASE / (1 + like)
    wki = freq * (1 - np.exp(-freq / (BOLTZWN * temp)))
    out = np.zeros(len(freq))
    if lam1 >= 0 and lam2 >= 0:
        jplus = jrange + max(lam1, lam2)
        for j1 in range(jrange):
            p1 = _ph2_vib(j1, temp) / q
            if p1 < 0.001:
                continue
            for jp1 in range(jplus):
                cg1 = _clebsqr(j1, lam1, jp1)
                if cg1 <= 0:
                    continue
                omega1 = _eh2_vib(nv1, jp1 * (jp1 + 1)) - _eh2_vib(0, j1 * (j1 + 1))
                for j2 in range(jrange):
                    p2 = _ph2_vib(j2, temp) / q
                    if p2 < 0.001:
                        continue
                    for jp2 in range(jplus):
                        cg2 = _clebsqr(j2, lam2, jp2)
                        if cg2 <= 0:
                            continue
                        omega2 = _eh2_vib(nv2, jp2 * (jp2 + 1)) - _eh2_vib(0, j2 * (j2 + 1))
                        fac = calib * p1 * p2 * cg1 * cg2
                        out += s * _bgama(freq - omega1 - omega2, t1, t2, e, t3, t4, temp) * (wki * fac)
    else:
        for j in range(jrange):
            p = _ph2_vib(j, temp) / q
            if p < 0.001:
                continue
            for jp in range(jrange + lam):
                cg = _clebsqr(j, lam, jp)
                if cg <= 0:
                    continue
                omega = _eh2_vib(nv1, jp * (jp + 1)) - _eh2_vib(0, j * (j + 1))
                out += s * _bgama(freq - omega, t1, t2, e, t3, t4, temp) * (wki * calib * p * cg)
    return out

def _y2(x, a, b, c):
    return a * np.exp((c * x + b) * x)

def _y3(x, a, b, c, d):
    return a * np.exp(b * x + c * x * x + d * x * x * x)

def h2_h2(freq, temp):
    freq = np.asarray(freq, dtype=float)
    q, jrange = _parts_vib(temp)
    x = np.log(temp)
    comps = (((2.881e-61, -1.1005, 0.131), (7.3485, -1.3874, 0.166), (7.883e-13, -0.3652, -0.0271), (3.803e-13, -0.4048, -0.0091), (1.0922e-12, -0.481, -0.0127), (5.174e-12, -0.9841, 0.0483), 0, 2, 2, 3), ((2.404e-64, -1.4161, 0.1847), (-0.8033, -0.4474, -0.0235), (3.873e-13, -0.4226, -0.0183), (2.743e-13, -0.3566, -0.014), (4.171e-13, -0.5223, 0.0097), (2.2725e-12, -1.1056, 0.0139), 4, 0, 4, 5), ((6.393e-63, -1.5964, 0.2359), (21.414, -1.2511, 0.1178), (1.876e-13, -0.4615, -0.0012), (4.839e-13, -0.5158, 0.0075), (4.55e-13, -0.5507, 0.0095), (2.045e-12, -0.5266, -0.024), 0, 2, 2, 1), ((5.965e-63, -1.0394, 0.1184), (6.674, -0.9459, 0.1124), (4.764e-13, -0.1725, -0.045), (4.016e-13, -0.3802, -0.0134), (1.0752e-12, -0.4617, -0.0085), (1.1405e-11, -1.2991, 0.0729), 2, 2, 3, 3))
    result = np.zeros(len(freq))
    for comp in comps:
        pars = [_y2(x, *part) for part in comp[:6]]
        result += _addspec_h2(freq, temp, q, jrange, *pars, comp[6], comp[7], comp[8], 0, 0, 1)
    return result

def _param_bc(temp, g0, g1, g2):
    tau0 = 1.05458875e-27 / (2 * 1.38054e-16 * temp)
    tta = g0 * tau0 / g1
    tau1 = np.sqrt((g2 * tta - g0 * (1 + tau0 ** 2 / tta)) / (g0 * (tau0 / tta) ** 2))
    return (tau1, tta / tau1)

def _param_k0(temp, g0, g1, g2):
    tau0 = 1.05458875e-27 / (2 * 1.38054e-16 * temp)
    delta = (tau0 * g1) ** 2 - 4 * (g1 * g1 / g0 + g1 / tau0 - g2) * tau0 * tau0 * g0
    if delta <= 0:
        raise ValueError(f'K0param: discriminante no positivo ({delta:g})')
    tau1 = (-np.sqrt(delta) - tau0 * g1) / (2 * (g1 * g1 / g0 + g1 / tau0 - g2))
    if tau1 <= 0:
        raise ValueError(f'K0param: tau1 no positivo ({tau1:g})')
    tau1 = np.sqrt(tau1)
    tau2 = tau0 * tau1 * g0 / (g1 * tau1 * tau1 - tau0 * g0)
    return (tau1, tau2)

def _addspec_he(freq, temp, q, jrange, g0, t1, t2, lam, lv, ifun):
    out = np.zeros(len(freq))
    wki = freq * (1 - np.exp(-freq / (BOLTZWN * temp)))
    for j in range(jrange):
        p = _ph2_vib(j, temp) / q
        if p < 0.001:
            continue
        for jp in range(jrange + lam):
            c = _clebsqr(j, lam, jp)
            if c <= 0:
                continue
            omega = _eh2_vib(0, jp * (jp + 1)) - _eh2_vib(0, j * (j + 1))
            frq = freq - omega
            bg = _bgama0(frq, t1, t2, temp) if ifun == 0 else _bgama1(frq, t1, t2, temp)
            out += g0 * bg * (wki * CALIB_BASE * p * c)
    return out

def h2_he(freq, temp):
    if not 40 < temp <= 3000:
        raise ValueError('ADDEMHE requires 40 < T <= 3000 K (as in the original Fortran).')
    freq = np.asarray(freq, dtype=float)
    q, jrange = _parts_vib(temp)
    x = np.log(temp)
    pieces = (((2.62e-63, -0.71411, 0.14587, -0.00352), (2.69e-50, -0.98315, 0.21988, -0.00729), (4.06e-36, -2.25664, 0.50098, -0.01925), 2, 3, 1), ((4.24e-65, -0.37075, 0.17473, -0.00489), (1.74e-50, -1.89232, 0.44399, -0.02029), (8.74e-36, -3.27717, 0.69166, -0.02865), 2, 1, 0), ((2.23e-62, -1.89198, 0.45505, -0.02238), (6.25e-49, -1.96486, 0.47059, -0.02402), (3.16e-34, -3.404, 0.72793, -0.03277), 0, 1, 0))
    total = np.zeros(len(freq))
    for a, b, c, lam, lv, ifun in pieces:
        g0, g1, g2 = (_y3(x, *t) for t in (a, b, c))
        t1, t2 = _param_bc(temp, g0, g1, g2) if ifun == 1 else _param_k0(temp, g0, g1, g2)
        total += _addspec_he(freq, temp, q, jrange, g0, t1, t2, lam, lv, ifun)
    return total

def load_corrch4(path):
    lines = [l for l in Path(path).read_text().splitlines() if l.strip()]
    if len(lines) < 42:
        raise ValueError(f'corrch4.txt incomplete: {len(lines)} records; expected 42')
    groups = []
    pos = 0
    for group in range(7):
        values = []
        for ii in range(6):
            line = lines[pos]
            pos += 1
            expected = 10 if ii < 5 else 1
            if len(line) < 8 * expected:
                raise ValueError(f'corrch4.txt: registro corto para grupo {group + 1}')
            values.extend((float(line[k * 8:(k + 1) * 8]) for k in range(expected)))
        groups.append(values)
    a = np.asarray(groups[:3], dtype=float).T
    b = np.asarray(groups[3:], dtype=float).T
    return (a, b)
B01_CH4 = float(np.float32(5.24))
B0_CH4 = float(np.float32(59.3392))
D0_CH4 = float(np.float32(0.04599))

def _parts_ch(temp):
    b0, d0 = (B0_CH4, D0_CH4)
    q = 0.0
    j = 0
    while True:
        ii = j * (j + 1)
        dq = (2 * j + 1) * (1 if j % 2 == 0 else 3) * np.exp(-1.4387859 * (b0 - ii * d0) * ii / temp)
        q += dq
        j += 1
        if dq <= q / 900.0:
            break
        if j > 500:
            raise RuntimeError('PARTSUMCH H2 no converge')
    jr1 = j
    q1 = 0.0
    j = 0
    while True:
        dq = np.exp(-1.4387859 * B01_CH4 * j * (j + 1) / temp) * (2 * j + 1) ** 2
        q1 += dq
        j += 1
        if dq <= q1 / 1000.0:
            break
        if j > 500:
            raise RuntimeError('PARTSUMCH CH4 no converge')
    return (q, jr1, q1, j + 3)

def _profile(rsibb, x, y, wnrmax, dx, ns, nsriup, slit):
    if y <= 0:
        return
    x0 = (nsriup + 1) / 2 + x / dx
    nc = int(x0)
    n1 = nc + 1
    slope = y / slit
    nu = max(1, int(x0 - ns))
    if nu > nsriup:
        return
    no = min(nsriup, int(x0 + ns))
    if no < 1:
        return
    nc = min(nc, nsriup)
    if nc > 1:
        for i in range(nu, nc + 1):
            xi = (i - 1) * dx - wnrmax
            dr = slope * (xi - (x - slit))
            if dr > 0:
                rsibb[i - 1] += dr
    if nc >= nsriup:
        return
    n1 = max(1, n1)
    for i in range(n1, no + 1):
        xi = (i - 1) * dx - wnrmax
        dr = y - slope * (xi - x)
        if dr > 0:
            rsibb[i - 1] += dr
_EB0 = np.array([-26.1308, -24.9681, -22.6532, -19.208, -14.669, -9.0915, -2.564, 4.7286, 12.4947])
_EB1 = np.array([-0.516, -0.1246, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
_BOUND_MATRIX = {'32': {(1, 4): 1.4127e-40, (2, 3): 1.8414e-40, (2, 5): 2.3414e-40, (3, 4): 1.8534e-40, (3, 6): 3.0894e-40, (4, 5): 2.1745e-40, (4, 7): 3.6564e-40, (5, 6): 2.4659e-40, (5, 8): 4.0032e-40, (6, 7): 2.6656e-40, (7, 8): 2.7376e-40, (8, 9): 2.6805e-40, (6, 9): 4.1233e-40}, '43': {(1, 5): 3.3038e-42, (2, 4): 4.5246e-42, (2, 6): 5.2213e-42, (3, 5): 4.2125e-42, (3, 7): 6.5899e-42, (4, 6): 4.6727e-42, (4, 8): 7.4203e-42, (5, 7): 5.0413e-42, (5, 9): 7.7286e-42, (6, 8): 5.1657e-42, (7, 9): 5.0297e-42}, '54': {(1, 6): 9.9184e-43, (2, 5): 1.4168e-42, (2, 7): 1.4969e-42, (3, 4): 1.6093e-42, (3, 6): 1.2455e-42, (4, 5): 1.2587e-42, (4, 7): 1.312e-42, (5, 6): 1.2714e-42, (5, 8): 1.3379e-42, (6, 7): 1.2818e-42, (7, 8): 1.2378e-42, (6, 9): 1.2963e-42, (8, 9): 1.1352e-42, (3, 8): 1.7927e-42, (4, 9): 1.9144e-42}}
_VIB_BOUND = {'32': [(3.9933e-42, 22.5287, 3, 2), (3.079e-42, 18.6921, 4, 1), (4.281e-42, 14.544, 5, 2)], '43': [(1.2153e-43, 19.0835, 4, 2), (8.9592e-44, 14.1528, 5, 1), (1.182e-43, 8.9669, 6, 2)], '54': [(4.5814e-44, 14.5442, 5, 2), (3.2581e-44, 8.5755, 6, 1), (4.1005e-44, 2.4395, 7, 2)]}

def _bound(temp, which, slit=4.3):
    wmax = {'32': 25.0, '43': 30.0, '54': 34.0}[which]
    nsri = 190
    nsriup = 2 * nsri + 1
    dx = wmax / nsri
    ns = int(slit / dx)
    rsibb = np.zeros(401, dtype=float)
    alpha = 1 / (0.69519 * temp)
    rm = 1.7768 * 1.672649e-24
    pf = (rm * 1.380662e-16 * temp * 2 * 3.141592654 / 6.626176e-27 ** 2) ** 1.5
    aa = _BOUND_MATRIX[which]

    def transition(l, j):
        a = aa.get((l, j), 0.0)
        shift = _EB0[j - 1] - _EB0[l - 1]
        pos = a * np.exp(-alpha * _EB0[l - 1]) / pf
        neg = a * np.exp(-alpha * _EB0[j - 1]) / pf
        _profile(rsibb, shift, pos, wmax, dx, ns, nsriup, slit)
        _profile(rsibb, -shift, neg, wmax, dx, ns, nsriup, slit)
    if which == '32':
        for l in range(1, 9):
            if l <= 6:
                transition(l, l + 3)
            transition(l, l + 1)
    elif which == '43':
        for l in range(1, 8):
            if l <= 5:
                transition(l, l + 4)
            transition(l, l + 2)
    else:
        for l in range(1, 9):
            if l <= 4:
                transition(l, l + 5)
            if l <= 6:
                transition(l, l + 3)
            transition(l, l + 1)
    for av, shift, n0, n1 in _VIB_BOUND[which]:
        _profile(rsibb, shift, av * np.exp(-alpha * _EB0[n0 - 1]) / pf, wmax, dx, ns, nsriup, slit)
        _profile(rsibb, -shift, av * np.exp(-alpha * _EB1[n1 - 1]) / pf, wmax, dx, ns, nsriup, slit)
    rsibb[:nsriup] /= 188365183000.0 * slit
    rsi = rsibb[nsri:nsriup].copy()
    initb = 1
    if which == '43':
        nz = np.flatnonzero(rsi)
        if nz.size:
            initb = int(nz[0] + 1)
            rsi = rsi[initb - 1:]
    om = (np.arange(len(rsi)) + initb - 1) * dx
    return (om, rsi, dx, wmax, initb)

def _greville_spline_second(x, y, eps=1e-05):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    if n < 2:
        raise ValueError('SPLINE requiere >=2 puntos')
    s2 = np.zeros(n, float)
    if n == 2:
        return s2
    d = np.diff(y) / np.diff(x)
    b = (x[1:-1] - x[:-2]) * 0.5 / (x[2:] - x[:-2])
    c = 3 * (d[1:] - d[:-1]) / (x[2:] - x[:-2])
    s2[1:-1] = c / 1.5
    for it in range(100000):
        eta = 0.0
        sm = np.max(np.abs(s2))
        epsi = eps * sm
        for i in range(1, n - 1):
            w = (c[i - 1] - b[i - 1] * s2[i - 1] - (0.5 - b[i - 1]) * s2[i + 1] - s2[i]) * 1.0717968
            eta = max(eta, abs(w))
            s2[i] += w
        if eta < epsi:
            break
    else:
        raise RuntimeError('SPLINE de Greville no converge')
    return s2

def _greville_eval(x, y, s2, t):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    s2 = np.asarray(s2, float)
    t = np.asarray(t, float)
    out = np.empty_like(t)
    n = len(x)
    low = t < x[0]
    high = t > x[-1]
    d = np.diff(y) / np.diff(x)
    if low.any():
        ht = t[low] - x[0]
        yp = d[0] + (x[0] - x[1]) * (2 * s2[0] + s2[1]) / 6
        out[low] = y[0] + yp * ht + s2[0] * ht * ht / 2
    if high.any():
        ht = t[high] - x[-1]
        yp = d[-1] + (x[-1] - x[-2]) * (s2[-2] + 2 * s2[-1]) / 6
        out[high] = y[-1] + yp * ht + s2[-1] * ht * ht / 2
    ok = ~(low | high)
    if ok.any():
        tt = t[ok]
        i = np.clip(np.searchsorted(x, tt, side='right') - 1, 0, n - 2)
        ht1 = tt - x[i]
        ht2 = tt - x[i + 1]
        prod = ht1 * ht2
        s3 = (s2[i + 1] - s2[i]) / (x[i + 1] - x[i])
        ss2 = s2[i] + ht1 * s3
        del_sqs = (s2[i] + s2[i + 1] + ss2) / 6
        out[ok] = y[i] + ht1 * d[i] + prod * del_sqs
    return out

class _BoundSpectrum:

    def __init__(self, temp, which):
        om, rsi, dx, wmax, initb = _bound(temp, which)
        self.x = om
        self.wmax = wmax
        if not np.all(rsi > 0):
            raise ValueError(f'BOUND{which} contains non-positive values; cannot take the logarithm')
        self.y = np.log(rsi * 1e+80)
        self.s2 = _greville_spline_second(om, self.y, 1e-05)

    def evaluate(self, frequency, beta):
        f = np.asarray(frequency, float)
        x = np.abs(f)
        tfac = np.where(f < 0, -beta * x, 0.0)
        res = np.empty_like(f)
        beyond = x > self.x[-1]
        if beyond.any():
            xn = x[beyond]
            yn = -(self.y[-2] - self.y[-1]) * (xn - self.x[-1]) / (self.x[-1] - self.x[-2]) + self.y[-1]
            res[beyond] = np.exp(yn + tfac[beyond]) * 1e-80
        if (~beyond).any():
            res[~beyond] = np.exp(tfac[~beyond] + _greville_eval(self.x, self.y, self.s2, x[~beyond])) * 1e-80
        return res

def _bgauch(fnu, g0, delta, omega0, temp):
    d = g0 / (delta * np.sqrt(6.2831853))
    desym = 2 / (1 + np.exp(-fnu / (0.6950304256 * temp)))
    fexp = np.minimum((fnu - omega0) ** 2 / (2 * delta * delta), 300.0)
    return d * desym * np.exp(-fexp)

def _parts_ch_p(j, temp, q):
    ii = j * (j + 1)
    return (2 * j + 1) * (1 if j % 2 == 0 else 3) * np.exp(-1.4387859 * (B0_CH4 - D0_CH4 * ii) * ii / temp) / q

def _addspec_h2ch(freq, temp, q, jrange, s, e, t1, t2, t3, t4, lam, lv, bound, *, scale=1.0, g0bf=0.0, delbf=1.0, om0=0.0):
    wki = freq * (1 - np.exp(-freq / (BOLTZWN * temp)))
    out = np.zeros(len(freq))
    for j in range(jrange):
        p = _parts_ch_p(j, temp, q)
        for jp in range(jrange + lam):
            cg = _clebsqr(j, lam, jp)
            if cg <= 0:
                continue
            ii = j * (j + 1)
            iip = jp * (jp + 1)
            omega = (B0_CH4 - D0_CH4 * iip) * iip - (B0_CH4 - D0_CH4 * ii) * ii
            frq = freq - omega
            xbg = s * _bgama(frq, t1, t2, e, t3, t4, temp, methane=True)
            valid = np.abs(frq) <= bound.wmax
            if valid.any():
                xbg[valid] += scale * bound.evaluate(frq[valid], 1 / (BOLTZWN * temp))
            if lv == 3 and g0bf != 0:
                xbg += _bgauch(frq, g0bf, delbf, om0, temp)
            out += xbg * (wki * CALIB_BASE * p * cg)
    return out

def _addspec_methane(freq, temp, q1, jrange, s, e, t1, t2, t3, t4, lam, bound, qw):
    wki = freq * (1 - np.exp(-freq / (BOLTZWN * temp)))
    beta = 1 / (BOLTZWN * temp)
    out = np.zeros(len(freq))
    for j in range(jrange):
        p = (2 * j + 1) * np.exp(-1.4387859 * B01_CH4 * j * (j + 1) / temp) / q1
        for ip in range(1, lam + 1):
            jp = j + ip
            omega = B01_CH4 * (jp * (jp + 1) - j * (j + 1))
            frq = freq - omega
            cc = 1 + qw[j, ip - 1] / 4.0
            fac = CALIB_BASE * p * cc
            xbg = s * _bgama(frq, t1, t2, e, t3, t4, temp, methane=True)
            valid = np.abs(frq) <= bound.wmax
            if valid.any():
                xbg[valid] += bound.evaluate(frq[valid], beta)
            out += xbg * wki * fac * (2 * jp + 1)
            jp = j - ip
            if jp < 0:
                continue
            omega = B01_CH4 * (jp * (jp + 1) - j * (j + 1))
            frq = freq - omega
            xbg = s * _bgama(frq, t1, t2, e, t3, t4, temp, methane=True)
            valid = np.abs(frq) <= bound.wmax
            if valid.any():
                xbg[valid] += bound.evaluate(frq[valid], beta)
            out += xbg * wki * (CALIB_BASE * p) * (2 * jp + 1)
        frq = freq
        xbg = s * _bgama(frq, t1, t2, e, t3, t4, temp, methane=True)
        valid = np.abs(frq) <= bound.wmax
        if valid.any():
            xbg[valid] += bound.evaluate(frq[valid], beta)
        out += xbg * wki * (CALIB_BASE * p) * (2 * j + 1)
    return out

def h2_ch4(freq, temp, corrch4):
    freq = np.asarray(freq, float)
    qw3, qw4 = load_corrch4(corrch4) if not isinstance(corrch4, tuple) else corrch4
    q, jrange, q1, jrange2 = _parts_ch(temp)
    if jrange2 > 51:
        raise ValueError('PARTSUMCH excede 51 filas de corrch4.txt')
    x = np.log(temp)
    tot = np.zeros(len(freq))
    bound32 = _BoundSpectrum(temp, '32')
    s, e, t1, t2, t3, t4 = [_y2(x, *a) for a in ((2.2916e-60, -1.27134, 0.12423), (-2.019, -0.02259, -0.05891), (4.8254e-14, 0.6421, -0.10109), (9.7826e-13, -0.48654, -0.0146), (3.5247e-13, 0.10164, -0.07879), (1.3961e-14, 1.11146, -0.09561))]
    g0bf = _y2(x, 7.3536e-73, -0.79815, -0.0585)
    delbf = _y2(x, 3.123, -0.00178, 0.00021)
    om0 = _y2(x, 8.6922, 0.00785, -0.00054)
    tot += _addspec_h2ch(freq, temp, q, jrange, s, e, t1, t2, t3, t4, 2, 3, bound32, scale=1.0, g0bf=g0bf, delbf=delbf, om0=om0)
    bound54 = _BoundSpectrum(temp, '54')
    s, e, t1, t2, t3, t4 = [_y2(x, *a) for a in ((2.7321e-60 / 9 * 0.038, -2.32012, 0.23082), (-1.8198, -0.00665, -0.05626), (1.3492e-13, 0.14472, -0.06506), (1.6596e-12, -0.77679, 0.01401), (5.9914e-13, -0.16208, -0.05895), (1.9405e-14, 0.95965, -0.10246))]
    tot += _addspec_h2ch(freq, temp, q, jrange, s, e, t1, t2, t3, t4, 4, 5, bound54, scale=0.038)
    bound43 = _BoundSpectrum(temp, '43')
    s, e, t1, t2, t3, t4 = [_y2(x, *a) for a in ((9.1196e-61 / 7, -1.56529, 0.15284), (-0.5, 0.0, 0.0), (5.1456e-14, 0.60523, -0.10856), (6.2514e-13, -0.51384, 0.00754), (5.5346e-13, -0.40381, 0.00208), (1.3804e-14, 1.9307, -0.27921))]
    tot += _addspec_methane(freq, temp, q1, jrange2, s, e, t1, t2, t3, t4, 3, bound43, qw3)
    s, e, t1, t2, t3, t4 = [_y2(x, *a) for a in ((2.7321e-60 / 9, -2.32012, 0.23082), (-1.8198, -0.00665, -0.05626), (1.3492e-13, 0.14472, -0.06506), (1.6596e-12, -0.77679, 0.01401), (5.9914e-13, -0.16208, -0.05895), (1.9405e-14, 0.95965, -0.10246))]
    tot += _addspec_methane(freq, temp, q1, jrange2, s, e, t1, t2, t3, t4, 4, bound54, qw4)
    return tot

def cia_components(freq, temp, corrch4):
    return (h2_h2(freq, temp), h2_he(freq, temp), h2_ch4(freq, temp, corrch4))
