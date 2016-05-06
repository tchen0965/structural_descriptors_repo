from pymatgen.core.structure import Structure
import effective_coordination as econ
import connectivity_from_structure as connectivity
import okeeffe_CN as OKeeffe
import unittest

__author__ = 'Tina'


class TestCoordinationMethods(unittest.TestCase):

    def setUp(self):

        """Loading structures before tests"""
        print "Loading structures from file"
        self.fe_structure = Structure.from_file('Fe.cif', True, True)
        self.caf2_structure = Structure.from_file('CaF2.cif', True, True)
        self.licoo2_structure = Structure.from_file('LiCoO2.cif', True, True)
        self.fe_finder = econ.EffectiveCoordFinder(self.fe_structure)
        self.caf2_finder = econ.EffectiveCoordFinder(self.caf2_structure)
        self.licoo2_finder = econ.EffectiveCoordFinder(self.licoo2_structure)
        print "Structures from file loaded"

    def tearDown(self):
        """Cleaning up after test"""

    def test_econ_Fe(self):
        """Test Routine econ for BCC Fe"""

        cns = self.fe_finder.get_avg_cn(anions=['Fe'])
        for cation in cns:
            self.assertEqual(cation, "Fe", "Only cation in BCC Fe should be Fe")
            self.assertAlmostEqual(cns[cation], 8, 1, "Fe should be, on average, 8-fold coordinated in BCC Fe")

    def test_econ_CaF2(self):
        """Test Routine econ for CaF2"""

        cns = self.caf2_finder.get_avg_cn()
        for cation in cns:
            self.assertEqual(cation, "Ca", "Only cation in CaF2 should be Ca")
            self.assertAlmostEqual(cns[cation], 8, 1, "Ca should be, on average, 8-fold coordinated in CaF2")

    def test_econ_LiCoO2(self):
        """Test Routine econ for LiCoO2"""

        cns = self.licoo2_finder.get_avg_cn()
        for cation in cns:
            self.assertTrue(cation == "Li" or cation == "Co", "Only cations in LiCoO2 should be Li and Co")
            self.assertAlmostEqual(cns[cation], 6, 1,
                                   "Both Li and Co should be, on average, 6-fold coordinated in LiCoO2")

    def test_OKeeffe_Fe(self):
        """Test Routine O'Keeffe CN for BCC Fe"""

        cns = OKeeffe.get_avg_cn(self.fe_structure)
        for cation in cns:
            self.assertEqual(cation, "Fe", "Only atoms in BCC Fe should be Fe")
            self.assertIsInstance(cns[cation], float, "Averaged cation should be a float")
            self.assertGreaterEqual(cns[cation], 8, "Fe coordination should be 8 (OKeeffe tends to overestimate)")

    def test_OKeeffe_CaF2(self):
        """Test Routine O'Keeffe CN for CaF2"""

        cns = OKeeffe.get_avg_cn(self.caf2_structure)
        for cation in cns:
            self.assertTrue(cation == "Ca" or cation == "F", "Only ions in CaF2 should be Ca and F")
            self.assertIsInstance(cns[cation], float, "Averaged cation should be a float")
            if cation == "Ca":
                self.assertGreaterEqual(cns[cation], 8, "Ca coordination should be 8 (OKeeffe tends to overestimate)")
            if cation == "F":
                self.assertGreaterEqual(cns[cation], 4, "F coordination should be 4 (OKeeffe tends to overestimate)")

    def test_OKeeffe_LiCoO2(self):
        """Test Routine O'Keeffe CN for LiCoO2"""

        cns = OKeeffe.get_avg_cn(self.licoo2_structure)
        for cation in cns:
            self.assertTrue(cation == "Li" or cation == "Co" or cation == "O",
                            "Only ions in LiCoO2 should be Li, Co, and O")
            self.assertIsInstance(cns[cation], float, "Averaged cation should be a float")
            if cation == "Li":
                self.assertGreaterEqual(cns[cation], 6, "Li coordination should be 6 (OKeeffe tends to overestimate)")
            if cation == "Co":
                self.assertGreaterEqual(cns[cation], 6, "Co coordination should be 6 (OKeeffe tends to overestimate)")
            if cation == "O":
                self.assertGreaterEqual(cns[cation], 6, "O coordination should be 6 (OKeeffe tends to overestimate)")

    def test_econ_on_polyhedra(self):
        """Test Routine econ for single Polyhedra"""

        licoo2_matrix, licoo2_polyhedra = connectivity.get_connectivity_matrix(self.licoo2_structure, False)
        for polyhedra in licoo2_polyhedra:
            if polyhedra.central_ion_name == "Li":
                test_poly = polyhedra
                break

        cn = econ.get_effective_cn(test_poly)
        self.assertAlmostEqual(cn, 6.0, 2, "Li polyhedra should be 6-fold coordinated")

if __name__ == '__main__':
    
    unittest.main()
