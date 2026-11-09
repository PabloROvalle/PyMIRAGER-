"""Physical constants for Jupiter."""
from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class JupiterConstants:
    h2: float = 0.8616
    he: float = 0.136
    ch4: float = 0.0024
    gas_constant: float = 8.31451
    avogadro: float = 6.02214e+23
    boltzmann: float = 1.38066e-23
    planck: float = 6.62608e-34
    cvel_cm_s: float = 29979200000.0
    rjup_m: float = 71492000.0
    j2: float = 14.736 / 1000.0
    j4: float = -587 / 1000000.0
    j6: float = 31 / 1000000.0
    rotation_days: float = 0.413542
    ellipticity: float = 0.0648744
    mass_1e24kg: float = 1898.18722
    grav_constant: float = 6.67408e-11

    @property
    def mean_molar_mass_kg_mol(self) -> float:
        return 0.0020158 * self.h2 + 0.0040026 * self.he + 0.016043 * self.ch4

C = JupiterConstants()
