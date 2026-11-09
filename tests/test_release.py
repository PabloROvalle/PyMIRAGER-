"""Fast release/regression tests; do not require external Jupiter .opa tables."""
import ast
import unittest
from dataclasses import dataclass
from pathlib import Path
from unittest.mock import patch

import numpy as np

import jupiter_rt
from jupiter_rt.cli import PRESETS, options_for_preset, parser
from jupiter_rt.inverse import run_retrieval, retrieval_step


@dataclass
class Param:
    label: str = 'aerosol'
    kind: str = 'aerosol'


class SyntheticModel:
    """Small deterministic inversion model to test constraints, not RT physics."""
    def __init__(self, nlevel=12):
        self.params = [Param()]
        self.case = type('CaseStub', (), {})()
        self.case.pressure = np.geomspace(1.0, .001, nlevel)
        self.case.taucloud = np.full(nlevel, 0.01)
        self.case.temperature = np.full(nlevel, 150.)
        self.case.profiles = np.empty((nlevel, 0))
        self.case.observed = np.array([1.02, 1.03, 1.04])
        self.case.error = np.ones(3) * 0.1

    def forward(self, with_jacobians=True):
        total = float(self.case.taucloud.sum())
        synthetic = 1 + np.array([.2, .3, .4]) * total
        if not with_jacobians:
            return synthetic
        k = np.ones((3, len(self.case.pressure), 1), dtype=float)
        k[:, :, 0] *= np.array([.2, .3, .4])[:, None]
        return synthetic, k


class ReleaseTests(unittest.TestCase):
    def test_version_is_consistent(self):
        import re
        project_file = Path(__file__).resolve().parents[1] / 'pyproject.toml'
        project_text = project_file.read_text()
        self.assertEqual(jupiter_rt.__version__, '0.4.7')
        self.assertRegex(project_text, r'(?m)^version = "0\.4\.7"$')
        self.assertRegex(project_text, r'(?m)^name = "jupiter-rt-retrieval"$')

    def test_presets_and_bounds(self):
        c = options_for_preset('aerosol_cutoff')
        self.assertEqual((c.aerosol_zero_below_index, c.aerosol_zero_from_index), (85, 120))
        self.assertEqual(c.filter_mode, 'none')
        self.assertEqual(c.c_aerosol, .75)
        self.assertEqual(c.factor, 3)
        self.assertIsNone(options_for_preset('aerosol').aerosol_zero_below_index)
        self.assertIsNone(options_for_preset('aerosol').aerosol_zero_from_index)

    def test_cli_version_without_case(self):
        import contextlib
        import io
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            with self.assertRaises(SystemExit) as err:
                parser().parse_args(['--version'])
        self.assertEqual(err.exception.code, 0)
        self.assertEqual(output.getvalue().strip(), 'PyMIRAGER 0.4.7')

    def test_fixed_zero_constraint_and_jacobian(self):
        model = SyntheticModel()
        options = options_for_preset('aerosol_cutoff', {
            'aerosol_zero_below_index': 3,
            'aerosol_zero_from_index': 8,
            'max_iterations': 2,
        })
        original_step = retrieval_step
        observed_jacobians = []

        def inspect_step(p, kk, *args):
            observed_jacobians.append(kk.copy())
            return original_step(p, kk, *args)

        with patch('jupiter_rt.inverse.retrieval_step', side_effect=inspect_step):
            result = run_retrieval(model, options)
        forbidden = np.r_[0:3, 8:12]
        self.assertTrue(np.array_equal(model.case.taucloud[forbidden], np.zeros(len(forbidden))))
        self.assertEqual(result.sigma.shape, (12, 1))
        self.assertTrue(np.array_equal(result.sigma[forbidden, 0], np.zeros(len(forbidden))))
        self.assertTrue(np.array_equal(result.averaging[forbidden, :, 0], np.zeros((len(forbidden), 12))))
        self.assertTrue(observed_jacobians)
        for kk in observed_jacobians:
            self.assertTrue(np.all(kk[:, forbidden, 0] == 0))
        self.assertTrue(np.all(model.case.taucloud[3:8] >= 0))

    def test_invalid_cutoffs_are_rejected(self):
        for lower, upper in [(5, 5), (8, 5), (-1, 8), (3, 12), (True, 8)]:
            with self.subTest(lower=lower, upper=upper):
                model = SyntheticModel()
                opts = options_for_preset('aerosol_cutoff', {
                    'aerosol_zero_below_index': lower,
                    'aerosol_zero_from_index': upper,
                })
                with self.assertRaises(ValueError):
                    run_retrieval(model, opts)

    def test_cutoffs_require_aerosols(self):
        model = SyntheticModel()
        model.params = [Param(label='T', kind='temperature')]
        opts = options_for_preset('aerosol_cutoff', {
            'aerosol_zero_below_index': 3,
            'aerosol_zero_from_index': 8,
        })
        with self.assertRaisesRegex(ValueError, 'require an aerosol'):
            run_retrieval(model, opts)


if __name__ == '__main__':
    unittest.main()
