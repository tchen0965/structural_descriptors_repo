__author__ = 'Tina'

from pymatgen.core.structure import Structure
from structural_descriptors_repo import connectivity_from_structure as connectivity
import unittest


class TestConnectivity(unittest.TestCase):

    def setUp(self):
        """Loading structures before tests"""
        #print "TestConnectivityMethods:setUp_"
        print "Loading structures from file"
        self.fe_structure = Structure.from_file('Fe.cif', True, True)
        self.caf2_structure = Structure.from_file('CaF2.cif', True, True)
        self.licoo2_structure = Structure.from_file('LiCoO2.cif', True, True)
        print "Structures from file loaded"

    def tearDown(self):
        """Cleaning up after test"""
        #print "TestConnectivityMethods:tearDown_"


    def test_connectivity_matrix_Fe(self):
        """Test Routine Connectivity Matrix for BCC Fe"""

        central_species = ['Fe']
        peripheral_species = ['Fe']
        fe_matrix, fe_polyhedra = connectivity.get_connectivity_matrix(self.fe_structure, False, 2.8, peripheral_species, central_species)
        self.assertIn('Fe', fe_matrix.keys(), "Fe polyhedra not found in BCC Fe matrix")
        self.assertNotIn('Fe1', fe_matrix.keys(), "Found Fe1 polyhedra instead of Fe polyhedra")
        self.assertEqual(fe_matrix['Fe']['Fe']['point'],  8, "Fe should be point-sharing")
        self.assertEqual(fe_matrix['Fe']['Fe']['edge'], 12, "Fe should be edge-sharing")
        self.assertEqual(fe_matrix['Fe']['Fe']['face'], 6, "Fe should be face-sharing")

        self.assertLessEqual(fe_matrix['Fe']['Fe']['face'], 6, "Face-sharing instances exceeds number of faces")

        for poly in fe_polyhedra:
            self.assertIsInstance(poly, connectivity.Polyhedra, "List of polyhedra includes a non-polyhedra element")


    def test_connectivity_matrix_sites_diff_Fe(self):
        """Test Routine Connectivity Matrix Differentiating Sites for BCC Fe"""

        central_species = ['Fe']
        peripheral_species = ['Fe']
        fe_matrix, fe_polyhedra = connectivity.get_connectivity_matrix(self.fe_structure, True, 2.8, peripheral_species, central_species)
        self.assertIn('Fe1', fe_matrix.keys(), "Fe1 polyhedra not found in BCC Fe matrix")
        self.assertNotIn('Fe', fe_matrix.keys(), "Found Fe polyhedra instead of Fe1 polyhedra")
        self.assertEqual(fe_matrix['Fe1']['Fe1']['point'], 8, "Fe1 should be point-sharing")
        self.assertEqual(fe_matrix['Fe1']['Fe1']['edge'], 12, "Fe1 should be edge-sharing")
        self.assertEqual(fe_matrix['Fe1']['Fe1']['face'], 6, "Fe1 should be face-sharing")

        self.assertLessEqual(fe_matrix['Fe1']['Fe1']['face'], 6, "Face-sharing instances exceeds number of faces")

        for poly in fe_polyhedra:
            self.assertIsInstance(poly, connectivity.Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_connectivity_matrix_CaF2(self):
        """Test Routine Connectivity Matrix for CaF2"""

        central_species = ['Ca']
        peripheral_species = ['F']
        caf2_matrix, caf2_polyhedra = connectivity.get_connectivity_matrix(self.caf2_structure, False, 2.8, peripheral_species, central_species)
        self.assertIn('Ca', caf2_matrix.keys(), "Ca polyhedra not found in CaF2 matrix")
        self.assertNotIn('F', caf2_matrix.keys(), "Found F polyhedra instead of Ca polyhedra")
        self.assertEqual(caf2_matrix['Ca']['Ca']['point'], 0, "Ca should not be point-sharing")
        self.assertEqual(caf2_matrix['Ca']['Ca']['edge'], 12, "Ca should be edge-sharing")
        self.assertEqual(caf2_matrix['Ca']['Ca']['face'], 0, "Ca should not be face-sharing")

        self.assertLessEqual(caf2_matrix['Ca']['Ca']['face'], 6, "Face-sharing instances exceeds number of faces")
        central_species = ['F']
        peripheral_species = ['Ca']

        f2ca_matrix, f2ca_polyhedra = connectivity.get_connectivity_matrix(self.caf2_structure, False, 2.8, peripheral_species, central_species)
        self.assertIn('F', f2ca_matrix.keys(), "F polyhedra not found in CaF2 matrix")
        self.assertNotIn('Ca', f2ca_matrix.keys(), "Found Ca polyhedra instead of F polyhedra")
        self.assertEqual(f2ca_matrix['F']['F']['point'], 32, "F should point-sharing")
        self.assertEqual(f2ca_matrix['F']['F']['edge'], 12, "F should be edge-sharing")
        self.assertEqual(f2ca_matrix['F']['F']['face'], 0, "F should not be face-sharing")


    def test_connectivity_matrix_sites_diff_CaF2(self):
        """Test Routine Connectivity Matrix Differentiating Sites for CaF2"""

        central_species = ['Ca']
        peripheral_species = ['F']
        caf2_matrix, caf2_polyhedra = connectivity.get_connectivity_matrix(self.caf2_structure, True, 2.8, peripheral_species, central_species)
        self.assertIn('Ca1', caf2_matrix.keys(), "Ca1 polyhedra not found in CaF2 matrix")
        self.assertNotIn('F1', caf2_matrix.keys(), "Found F1 polyhedra instead of Ca polyhedra")
        self.assertEqual(caf2_matrix['Ca1']['Ca1']['point'], 0, "Ca should not be point-sharing")
        self.assertEqual(caf2_matrix['Ca1']['Ca1']['edge'], 12, "Ca should be edge-sharing")
        self.assertEqual(caf2_matrix['Ca1']['Ca1']['face'], 0, "Ca should not be face-sharing")

        central_species = ['F']
        peripheral_species = ['Ca']
        f2ca_matrix, f2ca_polyhedra = connectivity.get_connectivity_matrix(self.caf2_structure, True, 2.8, peripheral_species, central_species)
        self.assertIn('F1', f2ca_matrix.keys(), "F1 polyhedra not found in CaF2 matrix")
        self.assertNotIn('Ca1', f2ca_matrix.keys(), "Found Ca1 polyhedra instead of F polyhedra")
        self.assertEqual(f2ca_matrix['F1']['F1']['point'], 12, "F1 should be point-sharing with F1")
        self.assertEqual(f2ca_matrix['F1']['F1']['edge'], 0, "F1 should not be edge-sharing with F1")
        self.assertEqual(f2ca_matrix['F1']['F1']['face'], 0, "F1 sholud not be face-sharing with F1")
        self.assertEqual(f2ca_matrix['F2']['F2']['point'], 12, "F2 should be point-sharing with F2")
        self.assertEqual(f2ca_matrix['F2']['F2']['edge'], 0, "F2 should not be edge-sharing with F2")
        self.assertEqual(f2ca_matrix['F2']['F2']['face'], 0, "F2 should not beface-sharing with F2")
        self.assertEqual(f2ca_matrix['F1']['F2']['point'], 4, "F1 should be point-sharing with F2")
        self.assertEqual(f2ca_matrix['F1']['F2']['edge'], 6, "F1 should be edge-sharing with F2")
        self.assertEqual(f2ca_matrix['F1']['F2']['face'], 0, "F1 should not be face-sharing with F2")
        self.assertEqual(f2ca_matrix['F1']['F2'], f2ca_matrix['F2']['F1'],
                         "F2Ca connectivity matrix not symmetric across diagonal")

    def test_connectivity_matrix_LiCoO2(self):
        """Test Routine Connectivity Matrix for LiCoO2"""

        licoo2_matrix, licoo2_polyhedra = connectivity.get_connectivity_matrix(self.licoo2_structure, False)
        self.assertIn('Li', licoo2_matrix.keys(), "Li not found in LiCoO2 matrix")
        self.assertIn('Co', licoo2_matrix.keys(), "Co not found in LiCoO2 matrix")
        self.assertEqual(licoo2_matrix['Li']['Li']['point'], 0, "Li should not be point-sharing")
        self.assertEqual(licoo2_matrix['Li']['Li']['edge'], 6, "Li should be edge-sharing")
        self.assertEqual(licoo2_matrix['Li']['Li']['face'], 0, "Li should not be face-sharing")
        self.assertEqual(licoo2_matrix['Co']['Co']['point'], 0, "Co should not be point-sharing")
        self.assertEqual(licoo2_matrix['Co']['Co']['edge'], 6, "Co should be edge-sharing")
        self.assertEqual(licoo2_matrix['Co']['Co']['face'], 0, "Co should not be face-sharing")
        self.assertEqual(licoo2_matrix['Li']['Co']['point'], 6, "Li and Co should be point-sharing")
        self.assertEqual(licoo2_matrix['Li']['Co']['edge'], 6, "Li and Co should be edge-sharing")
        self.assertEqual(licoo2_matrix['Li']['Co']['face'], 0, "Li and Co should not be face-sharing")
        self.assertEqual(licoo2_matrix['Li']['Co'], licoo2_matrix['Co']['Li'],
                         "Connectivity matrix should be symmetric")

        for poly in licoo2_polyhedra:
            self.assertIsInstance(poly, connectivity.Polyhedra, "List of polyhedra includes a non-polyhedra element")


    def test_connectivity_description(self):
        """Test Routine Connectivity Description on BCC Fe"""

        central_species = ['Fe']
        peripheral_species = ['Fe']
        fe_matrix, fe_polyhedra = connectivity.get_connectivity_matrix(self.fe_structure, False, 2.8, peripheral_species, central_species)
        fe_descriptions = connectivity.get_connectivity_description(fe_matrix, fe_polyhedra, self.fe_structure, False)
        for cation in fe_descriptions.keys():
            self.assertIn(cation, fe_matrix.keys())
            self.assertTrue(isinstance(fe_descriptions[cation], str)
                            or isinstance(fe_descriptions[cation], unicode),
                            "Descriptions are not type str or unicode")

    def test_connectivity_description_sites_diff(self):
        """Test Routine Connectivity Description on LiCoO2"""


        licoo2_matrix, licoo2_polyhedra = connectivity.get_connectivity_matrix(self.licoo2_structure, True)
        licoo2_descriptions = connectivity.get_connectivity_description(licoo2_matrix, licoo2_polyhedra, self.licoo2_structure, True)
        for cation in licoo2_descriptions.keys():
            in_matrix_keys = False
            for keys in licoo2_matrix.keys():
                if cation in keys:
                    in_matrix_keys = True
            self.assertTrue(in_matrix_keys, "Cation being described not in connectivity matrix")
            self.assertTrue(isinstance(licoo2_descriptions[cation], str)
                            or isinstance(licoo2_descriptions[cation], unicode),
                            "Descriptions are not type str or unicode")


    def test_surrounding_connectivity(self):
        """Test Routine Surrounding Connectivity on Li polyhedra of LiCoO2"""

        # Testing polyhedra that are not specified by site number
        licoo2_matrix, licoo2_polyhedra = connectivity.get_connectivity_matrix(self.licoo2_structure, False)
        for polyhedra in licoo2_polyhedra:
            if polyhedra.central_ion_name == "Li":
                test_poly = polyhedra
                break
        connected_polyhedra = connectivity.get_surrounding_connectivity(self.licoo2_structure, test_poly)

        for poly, connections in connected_polyhedra:
            self.assertIn(poly.central_ion_name, licoo2_matrix.keys(), "Polyhedra connected to species not in matrix")
            self.assertLessEqual(connections, 2, "Should not be any face-sharing")

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