__author__ = 'Tina'

from pymatgen.core.structure import Structure
from structural_descriptors_repo import effective_coordination as ECoN
from structural_descriptors_repo import connectivity_from_structure as connectivity
from structural_descriptors_repo import okeeffe_CN as OKeeffe
import unittest

class TestCoordinationMethods(unittest.TestCase):

    def setUp(self):

        """Loading structures before tests"""
        #print "TestCoordinationMethods:setUp_"
        print "Loading structures from file"
        self.fe_structure = Structure.from_file('Fe.cif', True, True)
        self.caf2_structure = Structure.from_file('CaF2.cif', True, True)
        self.licoo2_structure = Structure.from_file('LiCoO2.cif', True, True)
        self.fe_finder = ECoN.EffectiveCoordFinder(self.fe_structure)
        self.caf2_finder = ECoN.EffectiveCoordFinder(self.caf2_structure)
        self.licoo2_finder = ECoN.EffectiveCoordFinder(self.licoo2_structure)
        print "Structures from file loaded"

    def tearDown(self):
        """Cleaning up after test"""
        #print "TestCoordinationMethods:tearDown_"

    def test_ECoN_Fe(self):
        """Test Routine ECoN for BCC Fe"""

        CNs = self.fe_finder.get_avg_CN(anions=['Fe'])
        for cation in CNs:
            self.assertEqual(cation, "Fe", "Only cation in BCC Fe should be Fe")
            self.assertAlmostEqual(CNs[cation], 8, 1, "Fe should be, on average, 8-fold coordinated in BCC Fe")


    def test_ECoN_CaF2(self):
        """Test Routine ECoN for CaF2"""

        CNs = self.caf2_finder.get_avg_CN()
        for cation in CNs:
            self.assertEqual(cation, "Ca", "Only cation in CaF2 should be Ca")
            self.assertAlmostEqual(CNs[cation], 8, 1, "Ca should be, on average, 8-fold coordinated in CaF2")

    def test_ECoN_LiCoO2(self):
        """Test Routine ECoN for LiCoO2"""

        CNs = self.licoo2_finder.get_avg_CN()
        for cation in CNs:
            self.assertTrue(cation == "Li" or cation == "Co", "Only cations in LiCoO2 should be Li and Co")
            self.assertAlmostEqual(CNs[cation], 6, 1, "Both Li and Co should be, on average, 6-fold coordinated in LiCoO2")


    def test_OKeeffe_Fe(self):
        """Test Routine O'Keeffe CN for BCC Fe"""

        CNs = OKeeffe.get_avg_CN(self.fe_structure)
        for cation in CNs:
            self.assertEqual(cation, "Fe", "Only atoms in BCC Fe should be Fe")
            self.assertIsInstance(CNs[cation], float, "Averaged cation should be a float")
            self.assertGreaterEqual(CNs[cation], 8, "Fe coordination should be 8 (OKeeffe tends to overestimate)")


    def test_OKeeffe_CaF2(self):
        """Test Routine O'Keeffe CN for CaF2"""

        CNs = OKeeffe.get_avg_CN(self.caf2_structure)
        for cation in CNs:
            self.assertTrue(cation == "Ca" or cation == "F", "Only ions in CaF2 should be Ca and F")
            self.assertIsInstance(CNs[cation], float, "Averaged cation should be a float")
            if cation == "Ca":
                self.assertGreaterEqual(CNs[cation], 8, "Ca coordination should be 8 (OKeeffe tends to overestimate)")
            if cation == "F":
                self.assertGreaterEqual(CNs[cation], 4, "F coordination should be 4 (OKeeffe tends to overestimate)")


    def test_OKeeffe_LiCoO2(self):
        """Test Routine O'Keeffe CN for LiCoO2"""

        CNs = OKeeffe.get_avg_CN(self.licoo2_structure)
        for cation in CNs:
            self.assertTrue(cation == "Li" or cation == "Co" or cation == "O", "Only ions in LiCoO2 should be Li, Co, and O")
            self.assertIsInstance(CNs[cation], float, "Averaged cation should be a float")
            if cation == "Li":
                self.assertGreaterEqual(CNs[cation], 6, "Li coordination should be 6 (OKeeffe tends to overestimate)")
            if cation == "Co":
                self.assertGreaterEqual(CNs[cation], 6, "Co coordination should be 6 (OKeeffe tends to overestimate)")
            if cation == "O":
                self.assertGreaterEqual(CNs[cation], 6, "O coordination should be 6 (OKeeffe tends to overestimate)")


    def test_ECoN_on_polyhedra(self):
        """Test Routine ECoN for single Polyhedra"""

        licoo2_matrix, licoo2_polyhedra = connectivity.get_connectivity_matrix(self.licoo2_structure, False)
        for polyhedra in licoo2_polyhedra:
            if polyhedra.central_ion_name == "Li":
                test_poly = polyhedra
                break

        CN = ECoN.get_effective_CN(test_poly)
        self.assertAlmostEqual(CN, 6.0, 2, "Li polyhedra should be 6-fold coordinated")




if __name__ == '__main__':
    
    unittest.main()