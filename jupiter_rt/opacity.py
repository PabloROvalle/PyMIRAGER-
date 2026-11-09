"""Binary molecular opacity tables and spectral grids."""
from __future__ import annotations
import numpy as np
from pathlib import Path
from dataclasses import dataclass


@dataclass(frozen=True)
class OpacityGrid:
    origin: float = 500.0
    blocks: int = 150
    nfreq: int = 10000
    dnu: float = 0.001
    ntemp: int = 5
    nlevel: int = 361

    @property
    def block_width(self):
        return self.nfreq * self.dnu

    @property
    def high(self):
        return self.origin + self.blocks * self.block_width

SPECIAL = {'hcn': OpacityGrid(origin=680, blocks=7), 'co2': OpacityGrid(origin=640, blocks=7), 'h2o': OpacityGrid(origin=500, blocks=150), 'c6h6': OpacityGrid(origin=640, blocks=7)}

BASE = OpacityGrid()

def load_temperature_table():
    return np.load(Path(__file__).parent / 'data' / 'temperature_table.npy')

class BinaryOpacity:

    def __init__(self, path, grid: OpacityGrid=BASE, endian='<'):
        self.path = Path(path)
        self.grid = grid
        if not self.path.is_file():
            raise FileNotFoundError(f'Missing opacity file: {self.path}')
        self.dtype = np.dtype(endian + 'f4')
        self.stride = self.dtype.itemsize * grid.nfreq
        expected = grid.nlevel * grid.ntemp * grid.blocks * self.stride
        if self.path.stat().st_size < expected:
            raise ValueError(f'{self.path}: file too small for {grid.blocks} blocks and {grid.nlevel} levels; expected at least {expected:,} B')
        self._mem = np.memmap(self.path, dtype=self.dtype, mode='r')

    def _record(self, number):
        if number < 1 or number * self.grid.nfreq > self._mem.size:
            raise IndexError(f'Opacity record {number} out of range: {self.path.name}')
        ix = (number - 1) * self.grid.nfreq
        return np.asarray(self._mem[ix:ix + self.grid.nfreq], dtype=np.float64)

    def coefficients(self, k, block, table_temp, t):
        g = self.grid
        low = g.origin + block * g.block_width
        if not 0 <= block < g.blocks:
            raise IndexError('Spectral block out of range')
        tt = np.asarray(table_temp, dtype=float)
        if t <= tt[0]:
            it = 0
            frac = 0.0
        elif t >= tt[-1]:
            it = g.ntemp - 2
            frac = 1.0
        else:
            it = int(np.searchsorted(tt, t, side='right') - 1)
            frac = (t - tt[it]) / (tt[it + 1] - tt[it])
        rec = (it + k * g.ntemp) * g.blocks + (block + 1)
        a = self._record(rec)
        b = self._record(rec + g.blocks)
        if np.any(a < 0) or np.any(b < 0):
            raise ValueError(f'{self.path.name}: negative opacities at level k={k}, block={block}')
        return np.power(a, 1 - frac) * np.power(b, frac)

class OpacitySet:

    def __init__(self, gas_names, directory, *, metadata=None, endian='<'):
        directory = Path(directory)
        meta = {}
        if metadata is not None:
            import json
            raw = json.loads(Path(metadata).read_text())
            for n, v in raw.items():
                meta[n.lower()] = OpacityGrid(**v)
        self.tables = {}
        for name in gas_names:
            name = name.lower()
            grid = meta.get(name, SPECIAL.get(name, BASE))
            file = directory / f'jupiter_ch4_{name}.opa'
            self.tables[name] = BinaryOpacity(file, grid, endian)

    def lines(self, gas, k, grid_nu, table_t, t):
        table = self.tables[gas.lower()]
        g = table.grid
        out = np.zeros(len(grid_nu), float)
        starts = np.arange(0, len(grid_nu), g.nfreq)
        if len(grid_nu) % g.nfreq:
            raise ValueError('The internal grid must contain complete blocks')
        for start in starts:
            block_min = float(grid_nu[start])
            block_float = (block_min - g.origin) / g.block_width
            block = int(round(block_float))
            if 0 <= block < g.blocks and abs(block_float - block) < 1e-06:
                out[start:start + g.nfreq] = table.coefficients(k, block, table_t, t)
        return out

    def validate_active_coverage(self, gas_names, wave):
        for name in gas_names:
            g = self.tables[name].grid
            covered = np.logical_and(wave >= g.origin, wave < g.high)
            if not np.any(covered):
                raise ValueError(f'{name.upper()}: no observed wavenumber falls within [{g.origin}, {g.high}) cm^-1')
