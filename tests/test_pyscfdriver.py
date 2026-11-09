import importlib.util
import os
import json
import tempfile
from pathlib import Path
import unittest
from unittest import mock

import numpy as np
from pyscf import gto, dft


ROOT = Path(__file__).resolve().parents[1]


class DriverTests(unittest.TestCase):
    def setUp(self):
        spec = importlib.util.spec_from_file_location('example_driver', ROOT/'gromacs_example/pyscfdriver.py')
        self.driver = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.driver)
        d = self.driver
        d.QM_CHARGE = d.SYSTEM_CHARGE = 0
        d.QM_E_BASIS = 'sto-3g'
        d.QM_E_BASIS_AUX = 'weigend'
        d.DFT_E_XC = 'LDA,VWN'
        d.DFT_DF_NE = d.DFT_DF_NN = True
        d.SCF_CONV_TOL_GRAD = 1e-7
        d.PRINT_TIMING = False
        self.kinds = ['O', 'H', 'H']
        self.coords = [[0., 0., 0.], [0., -.757, .587], [0., .757, .587]]

    def make_mf(self, shift=0.):
        coords = np.array(self.coords)
        coords[1, 0] += shift
        mol = gto.M(atom=list(zip(self.kinds, coords)), basis='sto-3g', verbose=0)
        mf = dft.RKS(mol, xc='LDA,VWN').density_fit(auxbasis='weigend')
        mf.grids.level = 0
        return mf

    def call(self, step, mmcoords=None):
        if mmcoords is None:
            mmcoords = [[3., 0., 0.], [-3., 0., 0.]]
        return self.driver.qmmmCalc(
            step, 'sto-3g', 1, 0, [0, 1, 2], self.kinds, self.coords,
            [3, 4], ['H', 'H'], [.1, -.1], mmcoords, [])

    def test_selection_is_independent_and_excludes_links(self):
        d = self.driver
        for _ in range(3):
            atoms, selected = d.prepare_qmatoms(self.kinds + ['H'], self.coords + [[2., 0., 0.]],
                                               [0, 1, 2, 'L0'])
            self.assertEqual(selected, [1, 2])
            self.assertEqual(d.QM_NUC_INDEX, [])
            self.assertEqual(atoms[-1][0], 'H')
        d.QM_NUC_SELECT = 'custom'
        for selected in ([0], [-1], [3], [1, 1]):
            d.QM_NUC_INDEX = selected
            with self.assertRaises(ValueError):
                d.prepare_qmatoms(self.kinds, self.coords)
        d.QM_METHOD = 'DFT'
        atoms, selected = d.prepare_qmatoms(self.kinds, self.coords)
        self.assertEqual(selected, [])
        self.assertEqual([a[0] for a in atoms], self.kinds)

    def test_history_lifecycle_and_options(self):
        d = self.driver
        d.MD_EXTRAPOLATION = True
        for step, size in ((0, 1), (1, 2), (1, 2), (4, 1), (0, 1)):
            d.run_scf_extrapolation(self.make_mf(), step)
            self.assertEqual(len(d._extrapolator.history), size)
        with self.assertWarnsRegex(RuntimeWarning, 'coordinates changed'):
            d.run_scf_extrapolation(self.make_mf(.001), 0)
        self.assertEqual(len(d._extrapolator.history), 1)
        for option in ('aspc', dict(scheme='fock', history_length=2)):
            d.MD_EXTRAPOLATION = option
            d.run_scf_extrapolation(self.make_mf(), 1)
            d.run_scf_extrapolation(self.make_mf(.001), 2)
            self.assertEqual(len(d._extrapolator.history), 2)
        d.run_scf_extrapolation(self.make_mf(), -1)
        self.assertIsNone(d._extrapolator)
        d.MD_EXTRAPOLATION = 'misspelled'
        with self.assertRaises(ValueError):
            d.run_scf_extrapolation(self.make_mf(), 3)
        self.assertIsNone(d._extrapolator)

    def test_unconverged_scf_fails_before_forces(self):
        d = self.driver
        for option in (False, True):
            d.MD_EXTRAPOLATION = option
            mf = self.make_mf()
            with mock.patch.object(mf, 'kernel', return_value=0.), \
                 mock.patch.object(mf, 'nuc_grad_method') as gradient:
                with self.assertRaisesRegex(RuntimeError, 'SCF failed'):
                    d.run_scf_extrapolation(mf, 0)
                gradient.assert_not_called()
            self.assertIsNone(d._extrapolator)

    def test_cpu_qmmm_and_empty_cutoff(self):
        d = self.driver
        for method in ('cneo', 'DFT'):
            d.QM_METHOD = method
            d.MD_EXTRAPOLATION = False
            reference = self.call(0)
            d.MD_EXTRAPOLATION = True
            for step in (0, 1):
                predicted = self.call(step)
                for actual, expected in zip(predicted, reference):
                    np.testing.assert_allclose(actual, expected, atol=2e-6, rtol=0)
                for force, shape in zip(predicted[1:], ((3, 3), (2, 3))):
                    self.assertEqual(force.shape, shape)
                    self.assertEqual(force.dtype, np.float64)
                    self.assertTrue(force.flags.c_contiguous)
            self.assertEqual(d.QM_NUC_INDEX, [])
            self.coords[1][0] += .001
            moved = self.call(2)
            d.MD_EXTRAPOLATION = False
            fresh = self.call(2)
            for actual, expected in zip(moved, fresh):
                np.testing.assert_allclose(actual, expected, atol=2e-6, rtol=0)
            d.MD_EXTRAPOLATION = True
            d.QMMM_CUT = 5
            with self.assertWarnsRegex(RuntimeWarning, 'forces omit'):
                empty = self.call(3, [[100., 0., 0.], [-100., 0., 0.]])
            d.QMMM_CUT = None
            np.testing.assert_array_equal(empty[2], np.zeros((2, 3)))
            d.MD_EXTRAPOLATION = False
            pure = d.qmCalc(-1, 'sto-3g', 1, 0, self.kinds, self.coords)
            np.testing.assert_allclose(empty[0], pure[0], atol=2e-7, rtol=0)

    def test_point_charges_no_cutoff_and_timing(self):
        d = self.driver
        d.QM_METHOD = 'DFT'
        with tempfile.TemporaryDirectory() as directory:
            d.TIMING_FILE = str(Path(directory)/'timing.jsonl')
            with mock.patch.object(d.itrf, 'mm_charge', wraps=d.itrf.mm_charge) as decorate:
                self.call(0, [[6., 0., 0.], [-6., 0., 0.]])
                self.assertIsNone(decorate.call_args.args[3])
                np.testing.assert_allclose(decorate.call_args.args[1], [[6., 0., 0.], [-6., 0., 0.]])
            record, = [json.loads(line) for line in Path(d.TIMING_FILE).read_text().splitlines()]
            self.assertGreater(record['scf_cycles'], 0)
            self.assertGreater(record['scf_seconds'], 0)
            self.assertGreater(record['grad_seconds'], 0)
            self.assertTrue(record['converged'])

    @unittest.skipUnless(os.environ.get('TEST_DRIVER_GPU') == '1', 'GPU validation is opt-in')
    def test_gpu_qmmm(self):
        d = self.driver
        for method in ('cneo', 'DFT'):
            d.QM_METHOD = method
            d.USE_GPU = False
            d.MD_EXTRAPOLATION = False
            reference = self.call(0)
            d.USE_GPU = True
            d.MD_EXTRAPOLATION = True
            for step in (0, 1):
                actual = self.call(step)
                for x, y in zip(actual, reference):
                    np.testing.assert_allclose(x, y, atol=2e-5, rtol=0)
                self.assertTrue(actual[1].flags.c_contiguous)
                self.assertTrue(actual[2].flags.c_contiguous)


if __name__ == '__main__':
    unittest.main()
