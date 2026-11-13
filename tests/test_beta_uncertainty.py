"""Check batched endpoint MC against the previous per-energy calculation."""

from copy import deepcopy
import ctypes
import os
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('CONFLUX_DB', str(ROOT / 'data'))

import conflux.BetaEngine as beta
import conflux.BetaPlusEngine as beta_plus
from conflux.Basic import integrate_trapezoid


def make_branch(module, energies, **kwargs):
    parameters = dict(Z=38, A=90, I=0, Q=3., E0=3., sigma_E0=.03,
                      dJpi='2-', xbins=energies)
    if module is beta:
        parameters.update(frac=1., sigma_frac=0.)
    else:
        parameters.update(bfrac=1., sigma_bfrac=0., ecfrac=0., sigma_ecfrac=0.)
    parameters.update(kwargs)
    return module.BetaBranch(**parameters)


def scalar_reference(module, branch, energy, is_nu, samples=30):
    """Retain the original sampling, endpoint conversion and spectrum call order."""
    if branch.sigma_E0 == 0:
        return 0.
    endpoints = np.random.normal(branch.E0, branch.sigma_E0, samples)
    parameters = deepcopy(branch.Parameters)
    if module is beta:
        parameters['W0'] = endpoints / module.ELECTRON_MASS_MEV + 1
    else:
        parameters['W0'] = endpoints / module.ELECTRON_MASS_MEV - 2 * module.ELECTRON_MASS_MEV
    function = branch.custom_func or (module.neutrino if is_nu else module.electron)
    values = ((1 - branch.mixing) * function(energy, parameters)
              + branch.mixing * function(energy, parameters, numass=branch.numass))
    return np.std(values)


class BatchedUncertaintyTests(unittest.TestCase):
    def setUp(self):
        random_state = np.random.get_state()
        self.addCleanup(np.random.set_state, random_state)
        self.enterContext(patch.object(beta, 'libBSG', None))
        self.enterContext(patch.object(beta_plus, 'libBSG', None))

    def assert_matches_scalar(self, module, branch, energies, is_nu, samples=30):
        np.random.seed(123)
        expected = np.array([scalar_reference(module, branch, x, is_nu, samples)
                             for x in np.asarray(energies).flat]).reshape(np.shape(energies))
        expected_next_random = np.random.random(5)
        np.random.seed(123)
        actual = branch.SpectEndpointUncertMC(energies, is_nu, samples)
        actual_next_random = np.random.random(5)
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-14)
        np.testing.assert_array_equal(actual_next_random, expected_next_random)

    def test_python_batch_matches_scalar_for_both_particles_and_transitions(self):
        energies = np.array([0., .01, .1, 1., 1.9, 2., 2.5, 2.99, 3., 3.03, 4.])
        for module in (beta, beta_plus):
            for dJpi in ('0', '1-', '2-', '3'):
                for is_nu in (False, True):
                    with self.subTest(engine=module.__name__, dJpi=dJpi, neutrino=is_nu):
                        branch = make_branch(module, energies, dJpi=dJpi)
                        self.assert_matches_scalar(module, branch, energies, is_nu)

    def test_mixing_and_multidimensional_energy_arrays(self):
        energies = np.array([[.1, .4, 1.], [1.8, 2., 3.1]])
        for module in (beta, beta_plus):
            for is_nu in (False, True):
                branch = make_branch(module, energies.ravel(), numass=.04, mixing=.3)
                self.assert_matches_scalar(module, branch, energies, is_nu, samples=17)

    def test_each_bin_gets_independent_endpoint_samples(self):
        energies = np.array([.4, .4, .4, .4])
        branch = make_branch(beta, energies)
        np.random.seed(123)
        uncertainties = branch.SpectEndpointUncertMC(energies, True)
        self.assertGreater(np.ptp(uncertainties), 0)
        self.assert_matches_scalar(beta, branch, energies, True)

    def test_zero_endpoint_uncertainty_does_not_draw_random_samples(self):
        energies = np.array([0., .1, 1., 4.])
        for module in (beta, beta_plus):
            branch = make_branch(module, energies, sigma_E0=0.)
            self.assert_matches_scalar(module, branch, energies, True)
            self.assertEqual(branch.SpectEndpointUncertMC(.2, True), 0.)

    def test_scalar_inputs_and_empty_arrays(self):
        for module in (beta, beta_plus):
            branch = make_branch(module, np.array([0., 1.]))
            for is_nu in (False, True):
                self.assert_matches_scalar(module, branch, .3, is_nu)
                self.assert_matches_scalar(module, branch, np.array([]), is_nu)

    def test_custom_functions_keep_the_scalar_call_pattern(self):
        calls = []

        def custom(energy, parameters, numass=0):
            self.assertEqual(np.ndim(energy), 0)
            self.assertEqual(parameters['W0'].shape, (30,))
            calls.append(energy)
            return energy * parameters['W0'] + numass

        energies = np.array([.1, .2, .5])
        for module in (beta, beta_plus):
            branch = make_branch(module, energies, custom_func=custom, mixing=.2, numass=.04)
            self.assert_matches_scalar(module, branch, energies, True)
        self.assertTrue(calls)

    def test_normalized_spectra_and_uncertainties_match_the_previous_path(self):
        for module in (beta, beta_plus):
            for energies in (np.linspace(0, 4, 41), np.linspace(.1, 1., 20)):
                for is_nu in (False, True):
                    branch = make_branch(module, energies)
                    expected_spectrum = branch.BetaSpectrum(energies, is_nu)
                    np.random.seed(123)
                    expected_uncertainty = np.array([
                        scalar_reference(module, branch, x, is_nu) for x in energies])
                    if energies[0] <= 0 and energies[-1] >= branch.E0:
                        norm = integrate_trapezoid(energies, expected_spectrum)
                    else:
                        norm_grid = np.arange(0., branch.E0, branch.E0 / 50.)
                        norm = integrate_trapezoid(norm_grid, branch.BetaSpectrum(norm_grid, is_nu))
                    np.random.seed(123)
                    branch.CalcNormalizedSpectrum(is_nu)
                    np.testing.assert_allclose(branch.spectrum, expected_spectrum / norm,
                                               rtol=1e-12, atol=1e-14)
                    np.testing.assert_allclose(branch.uncertainty, expected_uncertainty / norm,
                                               rtol=1e-12, atol=1e-14)


class CompiledBatchedUncertaintyTests(BatchedUncertaintyTests):
    @classmethod
    def setUpClass(cls):
        path = os.environ.get('CONFLUX_TEST_BSG_LIBRARY')
        if path is None:
            paths = list((ROOT / 'miniBSG/build/lib').glob('libBSG.*'))
            if not paths:
                raise unittest.SkipTest('Compiled BSG library is unavailable')
            path = str(paths[0])
        cls.library = ctypes.CDLL(path)
        cls.library.BSG_beta_spectrum.argtypes = [ctypes.c_double, ctypes.c_double,
            ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_bool]
        cls.library.BSG_beta_spectrum.restype = ctypes.c_double

    def setUp(self):
        super().setUp()
        self.enterContext(patch.object(beta, 'libBSG', self.library))
        self.enterContext(patch.object(beta_plus, 'libBSG', self.library))


if __name__ == '__main__':
    unittest.main()
