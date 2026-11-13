"""Regression coverage for database transition decoding and spectral dispatch."""

import contextlib
import ctypes
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('CONFLUX_DB', str(ROOT / 'data'))

import numpy as np

import conflux.BetaEngine as beta
import conflux.BetaPlusEngine as beta_plus
from conflux.BetaTransition import decode_djpi
from conflux.bsg.CoulombFunctions import lambda_k
from conflux.bsg.SpectralFunctions import shape_factor_unique_forbidden


# dJpi, deltaJ, parity change, forbidden order, unique, spectral L
TRANSITIONS = [
    ('0', 0, False, 0, False, 0),
    ('1', 1, False, 0, False, 0),
    ('0-', 0, True, 1, False, 0),
    ('1-', 1, True, 1, False, 0),
    ('2-', 2, True, 1, True, 2),
    ('2', 2, False, 2, False, 0),
    ('3', 3, False, 2, True, 3),
    ('3-', 3, True, 3, False, 0),
    ('4-', 4, True, 3, True, 4),
    ('4', 4, False, 4, False, 0),
    ('5', 5, False, 4, True, 5),
]


def branch(module, **kwargs):
    args = (38, 90, 0, 3., 3., 0., 1., 0.)
    if module is beta_plus:
        args += (0., 0.)
    return module.BetaBranch(*args, xbins=np.linspace(.001, 1., 100), **kwargs)


def load_fixture(module):
    root = ET.Element('betaDB')
    isotope = ET.SubElement(root, 'isotope', isotope='380900', name='Sr90',
                            Q='3', HL='1', beta_frac='1')
    assignments = [t[0] for t in TRANSITIONS] + ['2-,3', '3,2-', '', '2-,unknown']
    for n, assignment in enumerate(assignments):
        ET.SubElement(isotope, 'branch', end_point_E=str(1 + n / 10),
                      sigma_E0='0', fraction='0.01', sigma_frac='0',
                      IB='0.01', sigma_IB='0', IE='0', sigma_IE='0', dJpi=assignment)
    engine_class = module.BetaEngine if module is beta else module.BetaPlusEngine
    with tempfile.NamedTemporaryFile(suffix='.xml') as db:
        ET.ElementTree(root).write(db.name)
        with contextlib.redirect_stdout(io.StringIO()):
            engine = engine_class(inputlist=[380900], targetDB=db.name,
                                  xbins=np.linspace(.001, 1, 100))
    return engine.istplist[380900]


class TransitionTests(unittest.TestCase):
    def test_both_database_loaders_decode_physical_transition_types(self):
        for module in (beta, beta_plus):
            isotope = load_fixture(module)
            for assignment, deltaJ, parity, order, unique, L in TRANSITIONS:
                with self.subTest(engine=module.__name__, dJpi=assignment):
                    b = next(b for b in isotope.branches.values() if b.dJpi == assignment)
                    self.assertEqual((b.deltaJ, b.parity_change, b.forbidden_order, b.unique),
                                     (deltaJ, parity, order, unique))
                    self.assertEqual(b.Parameters['L'], L)
                    expected = 'allowed' if order == 0 else (
                        'unique_forbidden' if unique else 'allowed_approximation')
                    self.assertEqual(b.shape_model, expected)

    def test_ambiguous_selection_is_independent_of_candidate_order(self):
        a, candidates, ambiguous = decode_djpi('2-,3')
        b, _, _ = decode_djpi('3,2-')
        self.assertEqual(a, b)
        self.assertEqual((a.deltaJ, a.forbidden_order, a.spectral_L), (2, 1, 2))
        self.assertTrue(ambiguous)
        self.assertEqual(len(candidates), 2)
        # At equal order use the nonunique approximation, not an arbitrary
        # unique shape picked by candidate ordering.
        selected, _, _ = decode_djpi('2-,1-')
        self.assertEqual(selected.shape_model, 'allowed_approximation')

    def test_unknown_assignments_never_become_high_order_shapes(self):
        for assignment in ('', 'unknown', '2-,unknown'):
            for module in (beta, beta_plus):
                with self.subTest(engine=module.__name__, assignment=assignment):
                    b = branch(module, dJpi=assignment)
                    self.assertIsNone(b.deltaJ)
                    self.assertIsNone(b.forbidden_order)
                    self.assertEqual(b.Parameters['L'], 0)
                    self.assertEqual(b.shape_model, 'allowed_approximation')

    def test_transition_edits_derive_metadata_and_update_parameters(self):
        for module in (beta, beta_plus):
            isotope = load_fixture(module)
            b = next(b for b in isotope.branches.values() if b.dJpi == '2-')
            isotope.EditBranch(b.E0, dJpi='1-')
            self.assertEqual(b.forbidden_order, 1)
            self.assertFalse(b.unique)
            self.assertEqual(b.Parameters['L'], 0)
            b.dJpi = '3'
            self.assertEqual((b.forbidden_order, b.unique, b.Parameters['L']), (2, True, 3))
            isotope.EditBranch(b.E0, E0=4.)
            replacement = isotope.branches[4.]
            self.assertEqual(replacement.dJpi, '3')
            self.assertEqual(replacement.Parameters['L'], 3)
            with self.assertRaises(TypeError):
                branch(module, forbiddenness=1)
            with self.assertRaises(TypeError):
                isotope.EditBranch(4., forbiddenness=1)
            for value in (-1, 0, 1):
                with self.assertRaises(TypeError):
                    isotope.EditBranch(4., dJpi=value)
            with self.assertRaises(AttributeError):
                replacement.forbidden_order = 1
            self.assertEqual(replacement.dJpi, '3')

    def test_edit_keywords_are_order_independent_and_preserve_other_parameters(self):
        def custom(e, p, **kwargs):
            return np.ones_like(e)

        for module in (beta, beta_plus):
            fraction = 'fraction' if module is beta else 'beta_fraction'
            changes = {'E0': 4., fraction: .25, 'dJpi': '2-', 'sigma_E0': .01,
                       'bAc': 5., 'custom_func': custom, 'numass': .1, 'mixing': .2}
            for edits in (changes, dict(reversed(list(changes.items())))):
                isotope = load_fixture(module)
                isotope.EditBranch(1., **edits)
                b = isotope.branches[4.]
                self.assertNotIn(1., isotope.branches)
                self.assertEqual(b.dJpi, '2-')
                self.assertEqual(b.Parameters['L'], 2)
                self.assertEqual(b.bAc, 5.)
                self.assertEqual(b.sigma_E0, .01)
                self.assertEqual((b.numass, b.mixing), (.1, .2))
                self.assertIs(b.custom_func, custom)
                isotope.EditBranch(4., dJpi='1-')
                self.assertEqual(b.bAc, 5.)
                self.assertIs(b.custom_func, custom)
                self.assertEqual(b.Parameters['L'], 0)
                isotope.EditBranch(4., custom_func=None)
                self.assertEqual(b.shape_model, 'allowed_approximation')
                isotope.EditBranch(5., **{fraction: .3, 'dJpi': '3', 'custom_func': custom})
                added = isotope.branches[5.]
                self.assertEqual(added.Parameters['L'], 3)
                self.assertIs(added.custom_func, custom)
                if module is beta_plus:
                    isotope.EditBranch(5., EC_fraction=.2)
                    self.assertAlmostEqual(added.fraction, .5)

    def test_endpoint_collision_leaves_existing_branches_intact(self):
        for module in (beta, beta_plus):
            isotope = load_fixture(module)
            original = isotope.branches.copy()
            with self.assertRaises(ValueError):
                isotope.EditBranch(1., E0=1.1, dJpi='2-')
            self.assertEqual(isotope.branches, original)
            self.assertEqual(isotope.branches[1.].dJpi, '0')

    def test_conversion_virtual_spectra_use_explicit_assignments(self):
        from conflux.ConversionEngine import VirtualBranch
        virtual = VirtualBranch.__new__(VirtualBranch)
        virtual.Zavg, virtual.Aavg = 38, 90
        energies = np.linspace(.01, 1., 100)
        with patch.object(beta, 'libBSG', None):
            allowed = virtual.BetaSpectrum(energies, 3., 1., dJpi='0', norm=False)
            unique = virtual.BetaSpectrum(energies, 3., 1., dJpi='2-', norm=False)
            expected = branch(beta, dJpi='2-').BetaSpectrum(energies) / (energies[1] - energies[0])
        np.testing.assert_allclose(unique, expected)
        self.assertFalse(np.allclose(unique, allowed))

    def test_nonunique_fallback_and_custom_shape_are_explicit(self):
        energies = np.array([.1, .2, .4])
        for module in (beta, beta_plus):
            with patch.object(module, 'libBSG', None):
                allowed = branch(module, dJpi="0").BetaSpectrum(energies)
                nonunique = branch(module, dJpi='1-')
                np.testing.assert_allclose(nonunique.BetaSpectrum(energies), allowed)
                self.assertEqual(nonunique.Parameters['shape_model'], 'allowed_approximation')
                seen = []

                def custom(e, p, **kwargs):
                    seen.append(p.copy())
                    return np.ones_like(e)

                nonunique.custom_func = custom
                np.testing.assert_array_equal(nonunique.BetaSpectrum(energies), np.ones_like(energies))
                self.assertEqual(nonunique.shape_model, 'custom')
                self.assertEqual(seen[0]['deltaJ'], 1)
                self.assertEqual(seen[0]['shape_model'], 'custom')

    def test_scalar_and_array_dispatch_pass_spin_change_to_compiled_backend(self):
        for module in (beta, beta_plus):
            for assignment, L in (('2-', 2), ('3', 3), ('1-', 0), ('', 0)):
                with self.subTest(engine=module.__name__, assignment=assignment):
                    fake_library = unittest.mock.Mock()
                    fake_library.BSG_beta_spectrum.return_value = 1.
                    with patch.object(module, 'libBSG', fake_library):
                        b = branch(module, dJpi=assignment)
                        b.BetaSpectrum(.2)
                        b.BetaSpectrum(np.array([.1, .2]))
                        b.BetaSpectrum(.2, nu_spectrum=True)
                    self.assertTrue(fake_library.BSG_beta_spectrum.called)
                    for call in fake_library.BSG_beta_spectrum.call_args_list:
                        self.assertEqual(call.args[5], L)
                        self.assertEqual(call.args[3], -38 if module is beta_plus else 38)

    def test_first_unique_shape_matches_neutrino_and_electron_momentum_terms(self):
        b = branch(beta, dJpi='2-')
        W = np.array([1.1, 1.5, 1.9])
        p = b.Parameters
        expected = (p['W0'] - W)**2 + lambda_k(W, p['Z'], p['R'], 2) * (W**2 - 1)
        np.testing.assert_allclose(shape_factor_unique_forbidden(W, **p), expected)
        for invalid in (-1, -2, 0, 1):
            with self.assertRaises(ValueError):
                shape_factor_unique_forbidden(W, invalid, p['W0'], p['Z'], p['R'])

    def test_sr90_database_regression(self):
        with contextlib.redirect_stdout(io.StringIO()):
            engine = beta.BetaEngine(inputlist=[380900], targetDB=str(
                ROOT / 'data/betaDB/ENSDF_betaDB_260707.xml'), xbins=np.linspace(.001, .545, 600))
        b = next(iter(engine.istplist[380900].branches.values()))
        self.assertEqual((b.dJpi, b.forbidden_order, b.unique, b.Parameters['L']), ('2-', 1, True, 2))
        with patch.object(beta, 'libBSG', None):
            spectrum = b.BetaSpectrum(b.xbins)
        mean = np.trapezoid(b.xbins * spectrum, b.xbins) / np.trapezoid(spectrum, b.xbins)
        self.assertAlmostEqual(mean, .1968, delta=.001)


class CompiledBSGTests(unittest.TestCase):
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
        cls.library.shape_factor_unique_forbidden.argtypes = [ctypes.c_double,
            ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_double]
        cls.library.shape_factor_unique_forbidden.restype = ctypes.c_double

    def test_compiled_unique_shape_agrees_with_python(self):
        for Z in (38, -38):
            for L in (2, 3, 4, 5):
                W = np.array([1.1, 1.5, 1.9])
                expected = shape_factor_unique_forbidden(W, L, 3., Z, .014)
                actual = [self.library.shape_factor_unique_forbidden(w, L, 3., Z, .014) for w in W]
                np.testing.assert_allclose(actual, expected, rtol=1e-10)

    def test_compiled_api_rejects_invalid_spin_changes(self):
        for L in (-1, -2, 1):
            self.assertTrue(np.isnan(self.library.BSG_beta_spectrum(1.5, 3., 90, 38, .014, L, False)))
        for L in (-1, 0, 1):
            self.assertTrue(np.isnan(self.library.shape_factor_unique_forbidden(1.5, L, 3., 38, .014)))

    def test_both_engines_use_compiled_fallback_and_unique_shapes(self):
        energies = np.linspace(.01, 1., 100)
        for module in (beta, beta_plus):
            with patch.object(module, 'libBSG', self.library):
                for is_nu in (False, True):
                    allowed = branch(module, dJpi="0").BetaSpectrum(energies, is_nu)
                    nonunique = branch(module, dJpi='1-').BetaSpectrum(energies, is_nu)
                    unique = branch(module, dJpi='2-').BetaSpectrum(energies, is_nu)
                    np.testing.assert_allclose(nonunique, allowed)
                    self.assertTrue(np.all(np.isfinite(unique)))
                    self.assertTrue(np.all(unique >= 0))
                    self.assertFalse(np.allclose(unique, allowed))


if __name__ == '__main__':
    unittest.main()
