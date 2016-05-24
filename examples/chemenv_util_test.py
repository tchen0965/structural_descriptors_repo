from pymatgen import Structure
from pymatgen.analysis.chemenv.coordination_environments import coordination_geometries_files as cg_files
from pymatgen.analysis.chemenv.coordination_environments import coordination_geometry_finder as polyfinder
from pymatgen.analysis.chemenv.coordination_environments import chemenv_strategies as strategies
from pymatgen.analysis.chemenv.coordination_environments import structure_environments as se
from chemenv_util import find_site_ce, find_species_string_ce, find_species_ce_from_light_se
import unittest
import json
import os


class TestChemEnvUtil(unittest.TestCase):

    def setUp(self):
        """Set up Structure from Cif file and path to jsons folder in Chemenv module in Pymatgen"""
        self.structure = Structure.from_file('LiCoO2.cif', True, False)
        self.structure.make_supercell([3, 3, 3])
        for isite, site in enumerate(self.structure._sites):
            if site.species_string == 'Li':
                self.first_site = site
                self.ifirst_site = isite
                break

        self.path_to_jsons = os.path.dirname(cg_files.__file__)

    def test_site_chemenv(self):
        """Test finding the coordination environment (mp symbol) of a single site in a structure using Chemenv"""
        ce = find_site_ce(self.structure, self.ifirst_site)
        with open(self.path_to_jsons+"/%s.json" % ce) as json_file:
            data = json.load(json_file)
        self.assertEqual(data['mp_symbol'], "O:6", "Li polyhedra should be 6-coordinated octahedron")
        self.assertEqual(data['coordination'], 6, "Li polyhedra should be 6-coordinated octahedron")
        self.assertEqual(data['name'], "Octahedron", "Li polyhedra should be 6-coordinated octahedron")

    def test_species_chemenv(self):
        """Test finding the coordination environments (mp symbols) of a specific species in a structure using Chemenv"""
        ces = find_species_string_ce(self.structure, self.first_site.species_string)

        for ce in ces:
            with open(self.path_to_jsons+"/%s.json" % ce) as json_file:
                data = json.load(json_file)
            self.assertEqual(data['mp_symbol'], "O:6", "Li polyhedra should be 6-coordinated octahedron")
            self.assertEqual(data['coordination'], 6, "Li polyhedra should be 6-coordinated octahedron")
            self.assertEqual(data['name'], "Octahedron", "Li polyhedra should be 6-coordinated octahedron")

    def test_species_from_light_structure_chemenv(self):
        """Test finding the ce's (mp symbols) of a species but starting from a Light Structure Environment object"""
        s1_finder = polyfinder.LocalGeometryFinder()
        s1_finder.setup_structure(self.structure)
        s1_finder.setup_parameters(centering_type='standard', structure_refinement='none')
        environments = s1_finder.compute_structure_environments_detailed_voronoi(maximum_distance_factor=1.5)

        light_se = se.LightStructureEnvironments(strategies.SimplestChemenvStrategy(), environments)

        ces = find_species_ce_from_light_se(light_se, self.first_site.species_string)

        for ce in ces:
            with open(self.path_to_jsons+"/%s.json" % ce) as json_file:
                data = json.load(json_file)
            self.assertEqual(data['mp_symbol'], "O:6", "Li polyhedra should be 6-coordinated octahedron")
            self.assertEqual(data['coordination'], 6, "Li polyhedra should be 6-coordinated octahedron")
            self.assertEqual(data['name'], "Octahedron", "Li polyhedra should be 6-coordinated octahedron")


if __name__ == '__main__':
    unittest.main()
