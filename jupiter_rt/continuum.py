"""Collision-induced absorption and continuum optical depths."""
from __future__ import annotations
import numpy as np
from pathlib import Path
from functools import lru_cache
from .constants import C
from .cia_python import cia_components, load_corrch4


@lru_cache(maxsize=768)
def _component_cached(t32, f0, nf, corr_file, mtime_ns, file_size):
    freq = np.float32(f0) + np.arange(nf, dtype=np.float32) * np.float32(2.0)
    factors = load_corrch4(corr_file)
    return cia_components(freq.astype(float), float(t32), factors)

class MissingCIABackend(RuntimeError):
    pass

def _cia_grid(nu):
    nu = np.asarray(nu, dtype=float)
    if nu.ndim != 1 or nu.size < 2 or np.any(np.diff(nu) <= 0):
        raise ValueError('CIA: nu must be a strictly increasing 1D grid with at least two points')
    start = float(nu[0]) - 2.0
    end = float(nu[-1]) + 2.0
    nin = int(np.ceil((end - start) / 2.0)) + 1
    if nin > 801:
        raise ValueError(f'CIA grid exceeds the original limit of 801 samples (nin={nin}). Split the spectrum into segments.')
    ftab = np.float32(start) + np.arange(nin, dtype=np.float32) * np.float32(2.0)
    return ftab

def _layer_scale(p, t, g, k):
    f = np.float32
    if k == 320:
        zn = f(f(f(273.16) * p[k] / t[k]) / f(1.01325))
        aux = f(f(f(f(p[k] * f(100000.0)) / (g[k] * f(C.mean_molar_mass_kg_mol))) * f(C.avogadro) / f(2.6868e+19)) / f(10000.0))
    else:
        zn = f(f(f(f(273.16) * f(2.0)) * f(np.sqrt(f(p[k] * p[k + 1])))) / f(t[k] + t[k + 1]) / f(1.01325))
        aux = f(f(f(f(f(p[k] - p[k + 1]) * f(100000.0)) / (g[k] * f(C.mean_molar_mass_kg_mol))) * f(C.avogadro) / f(2.6868e+19)) / f(10000.0))
    return (zn, aux)

def calculate_cia(case_dir, p, t, g, nu, *, library=None):
    p = np.asarray(p, dtype=np.float32)
    t = np.asarray(t, dtype=np.float32)
    g = np.asarray(g, dtype=np.float32)
    nu = np.asarray(nu, dtype=float)
    if p.shape != (361,) or t.shape != (361,) or g.shape != (361,):
        raise ValueError('The original CIA routine requires 361 levels for Jupiter')
    if np.any(p <= 0) or np.any(t <= 0) or np.any(g <= 0):
        raise ValueError('CIA requires positive pressure, temperature and gravity')
    ftab = _cia_grid(nu)
    corrpath = Path(case_dir) / 'corrch4.txt'
    if not corrpath.exists():
        raise FileNotFoundError(f'Missing corrch4.txt in {case_dir}')
    st = corrpath.stat()
    cache_key = (str(corrpath.resolve()), st.st_mtime_ns, st.st_size)
    h2 = np.float32(C.h2)
    he = np.float32(C.he)
    methane = np.float32(np.float32(1.0) - h2 - he)
    tab = np.zeros((len(ftab), 361), dtype=np.float32)
    for k in range(360):
        hyd, hel, ch4 = _component_cached(float(t[k]), float(ftab[0]), len(ftab), *cache_key)
        zn, aux = _layer_scale(p, t, g, k)
        prefactor = np.float32(np.float32(aux * h2) * zn)
        tab[:, k] = np.asarray(prefactor * (h2 * hyd + he * hel + methane * ch4), dtype=np.float32)
    lo = np.searchsorted(ftab, nu, side='right') - 1
    lo = np.clip(lo, 0, len(ftab) - 2)
    x = (nu - ftab[lo]) / (ftab[lo + 1] - ftab[lo])
    return (1 - x[:, None]) * tab[lo, :].astype(float) + x[:, None] * tab[lo + 1, :].astype(float)
