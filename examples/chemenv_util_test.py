from pymatgen import Structure
from pymatgen.analysis.chemenv.coordination_environments import coordination_geometries_files as cg_files
from chemenv_util import find_site_ce, find_species_string_ce
import unittest
import json
import os


class TestChemEnvUtil(unittest.TestCase):

    def setUp(self):
        self.structure = Structure.from_file('LiCoO2.cif', True, False)
        self.structure.make_supercell([3, 3, 3])
        for isite, site in enumerate(self.structure._sites):
            if site.species_string == 'Co':
                self.first_site = site
                self.ifirst_site = isite
                break

    def test_site_chemenv(self):
        ce = find_site_ce(self.structure, self.ifirst_site)
        path_to_jsons = os.path.dirname(cg_files.__file__)
        with open(path_to_jsons+"/%s.json" % ce) as json_file:
            data = json.load(json_file)
        self.assertEqual(data['mp_symbol'], "O:6", "Li polyhedra should be 6-coordinated octahedron")
        self.assertEqual(data['coordination'], 6, "Li polyhedra should be 6-coordinated octahedron")
        self.assertEqual(data['name'], "Octahedron", "Li polyhedra should be 6-coordinated octahedron")

    def test_species_chemenv(self):
        ces = find_species_string_ce(self.structure, self.first_site.species_string)

        path_to_jsons = os.path.dirname(cg_files.__file__)
        for ce in ces:
            with open(path_to_jsons+"/%s.json" % ce) as json_file:
                data = json.load(json_file)
            self.assertEqual(data['mp_symbol'], "O:6", "Li polyhedra should be 6-coordinated octahedron")
            self.assertEqual(data['coordination'], 6, "Li polyhedra should be 6-coordinated octahedron")
            self.assertEqual(data['name'], "Octahedron", "Li polyhedra should be 6-coordinated octahedron")

if __name__ == '__main__':
    unittest.main()
