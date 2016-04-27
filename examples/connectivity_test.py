__author__ = 'Tina'

from pymatgen.core.structure import Structure
from structural_descriptors_repo import effective_coordination as ECoN
from structural_descriptors_repo import connectivity_from_structure as connectivity
import unittest

# TODO: THROUGHOUT: use assertEqual and verify the exact value (or range using greater/less). Strong tests are better than weak tests. I pointed out a couple places where this needs fixing but there are more...

class TestConnectivity(unittest.TestCase):

    def setUp(self):
        """Loading structures before tests"""
        #print "TestConnectivityMethods:setUp_"
        print "Loading structures from file"
        self.FeStructure = Structure.from_file('Fe.cif', True, True)  # TODO: use Pythonic naming conventions. e.g. self.fe_structure. Note that the unittest class itself doesn't use Pythonic names, but we shouldn't copy them in this case.
        self.CaF2Structure = Structure.from_file('CaF2.cif', True, True)
        self.LiCoO2Structure = Structure.from_file('LiCoO2.cif', True, True)
        print "Structures from file loaded"

    def tearDown(self):
        """Cleaning up after test"""
        #print "TestConnectivityMethods:tearDown_"


    def test_connectivity_matrix_Fe(self):
        """Test Routine Connectivity Matrix for BCC Fe"""

        central_species = ['Fe']
        peripheral_species = ['Fe']
        FeMatrix, FePolyhedra = connectivity.get_connectivity_matrix(self.FeStructure, False, 2.8, peripheral_species, central_species)
        self.assertIn('Fe', FeMatrix.keys(), "Fe polyhedra not found in BCC Fe matrix")
        self.assertNotIn('Fe1', FeMatrix.keys(), "Found Fe1 polyhedra instead of Fe polyhedra")
        self.assertTrue(FeMatrix['Fe']['Fe']['point'] != 0, "Fe should be point-sharing")
        self.assertTrue(FeMatrix['Fe']['Fe']['edge'] != 0, "Fe should be edge-sharing")
        self.assertTrue(FeMatrix['Fe']['Fe']['face'] != 0, "Fe should be face-sharing")

        self.assertLessEqual(FeMatrix['Fe']['Fe']['face'], 6, "Face-sharing instances exceeds number of faces")

        for poly in FePolyhedra:
            self.assertIsInstance(poly, connectivity.Polyhedra, "List of polyhedra includes a non-polyhedra element")


    def test_connectivity_matrix_sites_diff_Fe(self):
        """Test Routine Connectivity Matrix Differentiating Sites for BCC Fe"""

        central_species = ['Fe']
        peripheral_species = ['Fe']
        FeMatrix, FePolyhedra = connectivity.get_connectivity_matrix(self.FeStructure, True, 2.8, peripheral_species, central_species)
        self.assertIn('Fe1', FeMatrix.keys(), "Fe1 polyhedra not found in BCC Fe matrix")
        self.assertNotIn('Fe', FeMatrix.keys(), "Found Fe polyhedra instead of Fe1 polyhedra")
        self.assertTrue(FeMatrix['Fe1']['Fe1']['point'] != 0, "Fe1 should be point-sharing")  # TODO: use assertEqual and verify the exact value (or range using greater/less). Strong tests are better than weak tests.
        self.assertTrue(FeMatrix['Fe1']['Fe1']['edge'] != 0, "Fe1 should be edge-sharing") # TODO: use assertEqual and verify the exact value (or range using greater/less). Strong tests are better than weak tests.
        self.assertTrue(FeMatrix['Fe1']['Fe1']['face'] != 0, "Fe1 should be face-sharing") # TODO: use assertEqual and verify the exact value (or range using greater/less). Strong tests are better than weak tests.

        self.assertLessEqual(FeMatrix['Fe1']['Fe1']['face'], 6, "Face-sharing instances exceeds number of faces")

        for poly in FePolyhedra:
            self.assertIsInstance(poly, connectivity.Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_connectivity_matrix_CaF2(self):
        """Test Routine Connectivity Matrix for CaF2"""

        central_species = ['Ca']
        peripheral_species = ['F']
        CaF2Matrix, CaF2Polyhedra = connectivity.get_connectivity_matrix(self.CaF2Structure, False, 2.8, peripheral_species, central_species)
        self.assertIn('Ca', CaF2Matrix.keys(), "Ca polyhedra not found in CaF2 matrix")
        self.assertNotIn('F', CaF2Matrix.keys(), "Found F polyhedra instead of Ca polyhedra")
        self.assertTrue(CaF2Matrix['Ca']['Ca']['point'] == 0, "Ca should not be point-sharing")
        self.assertTrue(CaF2Matrix['Ca']['Ca']['edge'] != 0, "Ca should be edge-sharing")
        self.assertTrue(CaF2Matrix['Ca']['Ca']['face'] == 0, "Ca should not be face-sharing")

        self.assertLessEqual(CaF2Matrix['Ca']['Ca']['face'], 6, "Face-sharing instances exceeds number of faces")
        central_species = ['F']
        peripheral_species = ['Ca']

        F2CaMatrix, F2CaPolyhedra = connectivity.get_connectivity_matrix(self.CaF2Structure, False, 2.8, peripheral_species, central_species)
        self.assertIn('F', F2CaMatrix.keys(), "F polyhedra not found in CaF2 matrix")
        self.assertNotIn('Ca', F2CaMatrix.keys(), "Found Ca polyhedra instead of F polyhedra")
        self.assertTrue(F2CaMatrix['F']['F']['point'] != 0, "F should point-sharing")
        self.assertTrue(F2CaMatrix['F']['F']['edge'] != 0, "F should be edge-sharing")
        self.assertTrue(F2CaMatrix['F']['F']['face'] == 0, "F should not be face-sharing")

        self.assertLessEqual(F2CaMatrix['F']['F']['face'], 4, "Face-sharing instances exceeds number of faces")


    def test_connectivity_matrix_sites_diff_CaF2(self):
        """Test Routine Connectivity Matrix Differentiating Sites for CaF2"""

        central_species = ['Ca']
        peripheral_species = ['F']
        CaF2Matrix, CaF2Polyhedra = connectivity.get_connectivity_matrix(self.CaF2Structure, True, 2.8, peripheral_species, central_species)
        self.assertIn('Ca1', CaF2Matrix.keys(), "Ca1 polyhedra not found in CaF2 matrix")
        self.assertNotIn('F1', CaF2Matrix.keys(), "Found F1 polyhedra instead of Ca polyhedra")
        self.assertTrue(CaF2Matrix['Ca1']['Ca1']['point'] == 0, "Ca should not be point-sharing")
        self.assertTrue(CaF2Matrix['Ca1']['Ca1']['edge'] != 0, "Ca should be edge-sharing")
        self.assertTrue(CaF2Matrix['Ca1']['Ca1']['face'] == 0, "Ca should not be face-sharing")

        central_species = ['F']
        peripheral_species = ['Ca']
        F2CaMatrix, F2CaPolyhedra = connectivity.get_connectivity_matrix(self.CaF2Structure, True, 2.8, peripheral_species, central_species)
        print F2CaMatrix
        self.assertIn('F1', F2CaMatrix.keys(), "F1 polyhedra not found in CaF2 matrix")
        self.assertNotIn('Ca1', F2CaMatrix.keys(), "Found Ca1 polyhedra instead of F polyhedra")
        self.assertTrue(F2CaMatrix['F1']['F1']['point'] != 0, "F1 should be point-sharing with F1")
        self.assertTrue(F2CaMatrix['F1']['F1']['edge'] == 0, "F1 should not be edge-sharing with F1")
        self.assertTrue(F2CaMatrix['F1']['F1']['face'] == 0, "F1 sholud not be face-sharing with F1")
        self.assertTrue(F2CaMatrix['F2']['F2']['point'] != 0, "F2 should be point-sharing with F2")
        self.assertTrue(F2CaMatrix['F2']['F2']['edge'] == 0, "F2 should not be edge-sharing with F2")
        self.assertTrue(F2CaMatrix['F2']['F2']['face'] == 0, "F2 should not beface-sharing with F2")
        self.assertTrue(F2CaMatrix['F1']['F2']['point'] != 0, "F1 should be point-sharing with F2")
        self.assertTrue(F2CaMatrix['F1']['F2']['edge'] != 0, "F1 should be edge-sharing with F2")
        self.assertTrue(F2CaMatrix['F1']['F2']['face'] == 0, "F1 should not be face-sharing with F2")
        self.assertEqual(F2CaMatrix['F1']['F2'], F2CaMatrix['F2']['F1'],
                         "F2Ca connectivity matrix not symmetric across diagonal")

    def test_connectivity_matrix_LiCoO2(self):
        """Test Routine Connectivity Matrix for LiCoO2"""

        LiCoO2Matrix, LiCoO2Polyhedra = connectivity.get_connectivity_matrix(self.LiCoO2Structure, False)
        self.assertIn('Li', LiCoO2Matrix.keys(), "Li not found in LiCoO2 matrix")
        self.assertIn('Co', LiCoO2Matrix.keys(), "Co not found in LiCoO2 matrix")
        self.assertTrue(LiCoO2Matrix['Li']['Li']['point'] == 0, "Li should not be point-sharing")
        self.assertTrue(LiCoO2Matrix['Li']['Li']['edge'] != 0, "Li should be edge-sharing")
        self.assertTrue(LiCoO2Matrix['Li']['Li']['face'] == 0, "Li should not be face-sharing")
        self.assertTrue(LiCoO2Matrix['Co']['Co']['point'] == 0, "Co should not be point-sharing")
        self.assertTrue(LiCoO2Matrix['Co']['Co']['edge'] != 0, "Co should be edge-sharing")
        self.assertTrue(LiCoO2Matrix['Co']['Co']['face'] == 0, "Co should not be face-sharing")
        self.assertTrue(LiCoO2Matrix['Li']['Co']['point'] != 0, "Li and Co should be point-sharing")
        self.assertTrue(LiCoO2Matrix['Li']['Co']['edge'] != 0, "Li and Co should be edge-sharing")
        self.assertTrue(LiCoO2Matrix['Li']['Co']['face'] == 0, "Li and Co should not be face-sharing")
        self.assertEqual(LiCoO2Matrix['Li']['Co'], LiCoO2Matrix['Co']['Li'],
                         "Connectivity matrix should be symmetric")

        for poly in LiCoO2Polyhedra:
            self.assertIsInstance(poly, connectivity.Polyhedra, "List of polyhedra includes a non-polyhedra element")


    def test_connectivity_description(self):
        """Test Routine Connectivity Description on BCC Fe"""

        central_species = ['Fe']
        peripheral_species = ['Fe']
        FeMatrix, FePolyhedra = connectivity.get_connectivity_matrix(self.FeStructure, False, 2.8, peripheral_species, central_species)
        FeDescriptions = connectivity.get_connectivity_description(FeMatrix, FePolyhedra, self.FeStructure, False)
        for cation in FeDescriptions.keys():
            self.assertIn(cation, FeMatrix.keys())
            self.assertTrue(isinstance(FeDescriptions[cation], str)
                            or isinstance(FeDescriptions[cation], unicode),
                            "Descriptions are not type str or unicode")

    def test_connectivity_description_sites_diff(self):
        """Test Routine Connectivity Description on LiCoO2"""


        LiCoO2Matrix, LiCoO2Polyhedra = connectivity.get_connectivity_matrix(self.LiCoO2Structure, True)
        LiCoO2Descriptions = connectivity.get_connectivity_description(LiCoO2Matrix, LiCoO2Polyhedra, self.LiCoO2Structure, True)
        for cation in LiCoO2Descriptions.keys():
            inMatrixKeys = False
            for keys in LiCoO2Matrix.keys():
                if cation in keys:
                    inMatrixKeys = True
            self.assertTrue(inMatrixKeys, "Cation being described not in connectivity matrix")
            self.assertTrue(isinstance(LiCoO2Descriptions[cation], str)
                            or isinstance(LiCoO2Descriptions[cation], unicode),
                            "Descriptions are not type str or unicode")


    def test_surrounding_connectivity(self):
        """Test Routine Surrounding Connectivity on Li polyhedra of LiCoO2"""

        # Testing polyhedra that are not specified by site number
        LiCoO2Matrix, LiCoO2Polyhedra = connectivity.get_connectivity_matrix(self.LiCoO2Structure, False)
        for polyhedra in LiCoO2Polyhedra:
            if polyhedra.central_ion_name == "Li":
                testPoly = polyhedra
                break
        connectedPolyhedra = connectivity.get_surrounding_connectivity(self.LiCoO2Structure, testPoly)

        for poly, connections in connectedPolyhedra:
            self.assertIn(poly.central_ion_name, LiCoO2Matrix.keys(), "Polyhedra connected to species not in matrix")
            self.assertTrue(connections <= 2, "Should not be any face-sharing")

        # Testing polyhedra that are specified by site number
        """
        LiCoO2Matrix2, LiCoO2Polyhedra2 = connectivity.get_connectivity_matrix(self.LiCoO2Structure, False)
        for polyhedra in LiCoO2Polyhedra2:
            if polyhedra.central_ion_name == "Li1":
                testPoly2 = polyhedra
                break
        connectedPolyhedra2 = connectivity.get_surrounding_connectivity(self.LiCoO2Structure, testPoly2)
        for poly, connections in connectedPolyhedra2:
            self.assertIn(poly.central_ion_name, LiCoO2Matrix.keys(), "Polyhedra connected to species not in matrix")
            self.assertTrue(connections <= 2, "Should not be any face-sharing")
        """



if __name__ == '__main__':

    unittest.main()