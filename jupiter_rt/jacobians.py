"""Molecular abundance Jacobian alternatives."""
from __future__ import annotations
import numpy as np


# Imported at call time to avoid a circular dependency with the forward model.

def _molecular_jacobian_exact(transmission, nu, t_layer, kappa, mu, convolution, *, method='analytical', fd_step=0.001):
    from .forward import planck
    nfreq, nlevel = transmission.shape
    if kappa.shape != (nfreq, nlevel) or len(t_layer) != nlevel - 1:
        raise ValueError('Invalid molecular Jacobian array dimensions')
    output = np.zeros((convolution.shape[0], nlevel))
    bprev = planck(nu, t_layer[0])
    partial = -bprev * transmission[:, 0]
    for k in range(nlevel):
        if 0 < k < nlevel - 1:
            b = planck(nu, t_layer[k])
            partial += (bprev - b) * transmission[:, k]
            bprev = b
        elif k == nlevel - 1:
            partial += bprev * transmission[:, k]
        scaled_tau = kappa[:, k] / mu[k]
        if method == 'analytical':
            spectrum_derivative = -partial * scaled_tau
        elif method == 'finite_difference':
            plus = -scaled_tau * np.expm1(fd_step)
            minus = -scaled_tau * np.expm1(-fd_step)
            plus = np.clip(plus, -700.0, 700.0)
            minus = np.clip(minus, -700.0, 700.0)
            spectrum_derivative = partial * (np.expm1(plus) - np.expm1(minus)) / (2 * fd_step)
        else:
            raise ValueError('Unknown molecular Jacobian method: ' + str(method))
        output[:, k] = convolution @ spectrum_derivative
    return output
