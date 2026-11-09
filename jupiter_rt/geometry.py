"""Jupiter atmospheric geometry and gravity."""
from __future__ import annotations
import numpy as np
from .constants import C


def _legpol(n: int, z: float) -> float:
    if n == 0:
        return 1.0
    if n == 1:
        return float(z)
    p0, p1 = (1.0, float(z))
    for i in range(2, n + 1):
        p0, p1 = (p1, ((2 * i - 1) * z * p1 - (i - 1) * p0) / i)
    return p1

def gravity(lat_deg: float, h_km: float) -> tuple[float, float]:
    xgm = C.mass_1e24kg * C.grav_constant * 1e+30
    coeff = (C.j2, C.j4, C.j6)
    xradius = C.rjup_m * 100
    e = 1 / (1 - C.ellipticity)
    omega = 2 * np.pi / (C.rotation_days * 24 * 3600)
    lat = np.deg2rad(lat_deg)
    latc = np.arctan(np.tan(lat) / e ** 2)
    ss, cc = (np.sin(latc), np.cos(latc))
    if abs(cc) < 1e-12:
        raise ValueError('NEWGRAV is unstable at exact polar latitudes')
    rr = np.sqrt(cc ** 2 + e ** 2 * ss ** 2)
    r = (xradius + h_km * 100000.0) / rr
    radkm = xradius / rr * 1e-05
    pol = [_legpol(i, ss) for i in range(1, 7)]
    gr = 1 - sum(((2 * i + 1) * rr ** (2 * i) * coeff[i - 1] * pol[2 * i - 1] for i in range(1, 4)))
    gr = gr * xgm / r ** 2 - r * omega ** 2 * cc ** 2
    gt = -sum((4 * i ** 2 * rr ** (2 * i) * coeff[i - 1] * (pol[2 * i - 2] - ss * pol[2 * i - 1]) / cc for i in range(1, 4)))
    gt = gt * xgm / r ** 2 + r * omega ** 2 * cc * ss
    return (float(np.hypot(gr, gt) * 0.01), float(radkm))

def atmosphere_geometry(pressure_bar, temperature_k, lat_deg, mu_emission):
    p = np.asarray(pressure_bar, dtype=float)
    t = np.asarray(temperature_k, dtype=float)
    if p.size != 361 or p.ndim != 1:
        raise ValueError('The original Fortran model requires 361 atmospheric levels')
    if not np.all(p[:-1] > p[1:]):
        raise ValueError('Pressures must be strictly decreasing')
    if not 0 < mu_emission <= 1:
        raise ValueError('Emission-angle cosine must be in (0,1]')
    z = np.zeros(p.size)
    g = np.zeros_like(z)
    for k in range(40, 360):
        g[k], _ = gravity(lat_deg, z[k] / 1000)
        z[k + 1] = z[k] + C.gas_constant * (t[k] + t[k + 1]) / 2 / C.mean_molar_mass_kg_mol / g[k] * np.log(p[k] / p[k + 1])
    g[360], radius = gravity(lat_deg, z[360] / 1000)
    for k in range(40, 0, -1):
        g[k], _ = gravity(lat_deg, z[k] / 1000)
        z[k - 1] = z[k] + C.gas_constant * (t[k] + t[k - 1]) / 2 / C.mean_molar_mass_kg_mol / g[k] * np.log(p[k] / p[k - 1])
    g[0], _ = gravity(lat_deg, z[0] / 1000)
    inside = 1 - (radius / (radius + z / 1000) * np.sqrt(1 - mu_emission ** 2)) ** 2
    if np.any(inside <= 0):
        raise ValueError('Invalid tangent-path geometry; check the emission angle')
    return (z, g, np.sqrt(inside))
