from pymatgen.core.structure import Structure
from connectivity_from_structure import get_connectivity_matrix, Polyhedra
from effective_coordination import EffectiveCoordFinder
import unittest

__author__ = 'Tina'


"""

Each unit test checks whether the connectivities between polyhedra as well as the coordination number match those in
the literature ("The Major Ternary Structural Families", O. Muller, R. Roy; "Transition Metal Oxides", Rao Raveau)

"""


class TestVariousStructures(unittest.TestCase):

    def test_barite(self):

        """Testing coordination and connectivity matrix for barite structure BaSO4"""

        barite_structure = Structure.from_file('test_structures/barite.cif', True, True)
        cn_finder = EffectiveCoordFinder(barite_structure)
        cns = cn_finder.get_avg_cn(anions=['O'])
        for cation in cns:
            self.assertTrue(cation == "Ba" or cation == "S",
                            "Ba and SO4 polyanions should be the only ions in barite")
            if cation == 'Ba':
                self.assertEqual(round(cns[cation]), 8, "Ba should be 8-fold coordinated")
            if cation == 'S':
                self.assertEqual(round(cns[cation]), 4, "S should be 4-fold coordinated")

        central_species = ['Ba', 'S']
        peripheral_species = ['O']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(barite_structure, False, 3.0, peripheral_species, central_species)

        self.assertIn('Ba', connectivity_matrix.keys(), "Ba polyhedra not found in barite matrix")
        self.assertIn('S', connectivity_matrix.keys(), "SO4 polyanions not found in barite matrix")
        self.assertEqual(connectivity_matrix['Ba']['Ba']['point'], 4, "Ba should be point-sharing")
        self.assertEqual(connectivity_matrix['Ba']['Ba']['edge'], 4, "Ba should be edge--sharing")
        self.assertEqual(connectivity_matrix['Ba']['Ba']['face'], 0, "Ba should be not be face-sharing")
        self.assertEqual(connectivity_matrix['S']['S']['point'], 0, "S should be isolated")
        self.assertEqual(connectivity_matrix['S']['S']['edge'], 0, "S should be isolated")
        self.assertEqual(connectivity_matrix['S']['S']['face'], 0, "S should be isolated")
        self.assertEqual(connectivity_matrix['Ba']['S']['point'], 6, "Ba and S should be point-sharing")
        self.assertEqual(connectivity_matrix['Ba']['S']['edge'], 1, "Ba and S should be edge-sharing")
        self.assertEqual(connectivity_matrix['Ba']['S']['face'], 0, "Ba and S should not be face-sharing")

        for connectivity_type in connectivity_matrix['Ba']['S'].keys():
            self.assertEqual(connectivity_matrix['Ba']['S'][connectivity_type]*barite_structure.composition['Ba'],
                             connectivity_matrix['S']['Ba'][connectivity_type]*barite_structure.composition['S'],
                             "Total number of sharing instances between 'Ba' and 'S' "
                             "should be same as total number of sharing"
                             "instances between 'S' and 'Ba'")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_k2so4(self):
        """Testing coordination and connectivity matrix for beta-K2SO4 structure"""

        k2so4_structure = Structure.from_file('test_structures/beta-K2SO4.cif', True, True)
        cn_finder = EffectiveCoordFinder(k2so4_structure)
        cns = cn_finder.get_avg_cn(radius=3.2, anions=['O'])
        for cation in cns:
            self.assertTrue(cation == "K" or cation == "S", "Ba and S should be only ions in beta-K2SO4")
            if cation == 'K':
                self.assertEqual(round(cns[cation]), 8, "K should be 8-fold coordinated")
            if cation == 'S':
                self.assertEqual(round(cns[cation]), 4, "S should be 4-fold coordinated")

        central_species = ['K', 'S']
        peripheral_species = ['O']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(k2so4_structure, False, 3.2, peripheral_species, central_species)
        self.assertIn('K', connectivity_matrix.keys(), "K polyhedra not found in beta-K2SO4 matrix")
        self.assertIn('S', connectivity_matrix.keys(), "SO4 polyanions not found in beta-K2SO4 matrix")
        self.assertEqual(connectivity_matrix['K']['K']['point'], 6, "K should be point-sharing")
        self.assertEqual(connectivity_matrix['K']['K']['edge'], 4, "K should be edge-sharing")
        self.assertEqual(connectivity_matrix['K']['K']['face'], 6, "K should be face-sharing")
        self.assertEqual(connectivity_matrix['S']['S']['point'], 0, "S should be isolated")
        self.assertEqual(connectivity_matrix['S']['S']['edge'], 0, "S should be isolated")
        self.assertEqual(connectivity_matrix['S']['S']['face'], 0, "S should be isolated")
        self.assertEqual(connectivity_matrix['S']['K']['point'], 4, "S and K should be point-sharing")
        self.assertEqual(connectivity_matrix['S']['K']['edge'], 7, "S and K should be edge-sharing")
        self.assertEqual(connectivity_matrix['S']['K']['face'], 0, "S and K should not be face-sharing")

        for connectivity_type in connectivity_matrix['K']['S'].keys():
            self.assertEqual(connectivity_matrix['K']['S'][connectivity_type]*k2so4_structure.composition['K'],
                             connectivity_matrix['S']['K'][connectivity_type]*k2so4_structure.composition['S'],
                             "Total number of sharing instances between 'K' and 'S' "
                             "should be same as total number of sharing"
                             "instances between 'S' and 'K'")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_fluorite(self):
        """Testing coordination and connectivity matrix for fluorite structure CaF2"""

        caf2_structure = Structure.from_file('test_structures/fluorite.cif', True, True)

        cn_finder = EffectiveCoordFinder(caf2_structure)
        cns = cn_finder.get_avg_cn(radius=2.6, anions=['F'])
        for cation in cns:
            self.assertEqual(cation, "Ca", "Ca should be the only ions in CaF2 fluorite")
            if cation == 'Ca':
                self.assertEqual(round(cns[cation]), 8, "Ca should be 8-fold coordinated")

        central_species = ['Ca']
        peripheral_species = ['F']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(caf2_structure, False, 2.6, peripheral_species, central_species)

        self.assertIn('Ca', connectivity_matrix.keys(), "Ca polyhedra not found in CaF2 fluorite matrix")
        self.assertEqual(connectivity_matrix['Ca']['Ca']['point'], 0, "Ca should not be point-sharing")
        self.assertEqual(connectivity_matrix['Ca']['Ca']['edge'], 12, "Ca should be edge-sharing")
        self.assertEqual(connectivity_matrix['Ca']['Ca']['face'], 0, "Ca should not be face-sharing")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_cafe2o4(self):
        """Testing coordination and connectivity matrix for CaFe2O4 structure"""

        cafe2o4_structure = Structure.from_file('test_structures/CaFe2O4.cif', True, True)

        cn_finder = EffectiveCoordFinder(cafe2o4_structure)
        cns = cn_finder.get_avg_cn(radius=2.6, anions=['O'])
        for cation in cns:
            self.assertTrue(cation == "Ca" or cation == "Fe",
                            "Ca and Fe should be the only cations in CaFe2O4 structure")
            if cation == 'Ca':
                self.assertLessEqual(round(cns[cation]), 8,
                                     "Ca should be 8-fold coordinated (effective coordination underestimates)")
            if cation == 'Fe':
                self.assertEqual(round(cns[cation]), 6, "Fe should be 6-fold coordinated")

        central_species = ['Ca', 'Fe']
        peripheral_species = ['O']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(cafe2o4_structure, False, 2.6, peripheral_species, central_species)

        self.assertIn('Ca', connectivity_matrix.keys(), "Ca cations not found in matrix")
        self.assertIn('Fe', connectivity_matrix.keys(), "Fe cations not found in matrix")
        self.assertEqual(connectivity_matrix['Ca']['Ca']['point'], 2, "Ca should be point-sharing")
        self.assertEqual(connectivity_matrix['Ca']['Ca']['edge'], 0, "Ca should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Ca']['Ca']['face'], 2, "Ca should be face-sharing")
        self.assertEqual(connectivity_matrix['Fe']['Fe']['point'], 4, "Fe should be point-sharing")
        self.assertEqual(connectivity_matrix['Fe']['Fe']['edge'], 4, "Fe should be edge-sharing")
        self.assertEqual(connectivity_matrix['Fe']['Fe']['face'], 0, "Fe should not be face-sharing")
        self.assertEqual(connectivity_matrix['Ca']['Fe']['point'], 8, "Fe and Ca should be point-sharing")
        self.assertEqual(connectivity_matrix['Ca']['Fe']['edge'], 5, "Fe and Ca should be edge-sharing")
        self.assertEqual(connectivity_matrix['Ca']['Fe']['face'], 2, "Fe and Ca should be face-sharing")

        for connectivity_type in connectivity_matrix['Ca']['Fe'].keys():
            self.assertEqual(connectivity_matrix['Ca']['Fe'][connectivity_type]*cafe2o4_structure.composition['Ca'],
                             connectivity_matrix['Fe']['Ca'][connectivity_type]*cafe2o4_structure.composition['Fe'],
                             "Total number of sharing instances between 'Ca' and 'Fe' "
                             "should be same as total number of sharing"
                             "instances between 'Fe' and 'Ca'")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_calcite(self):
        """Testing coordination and connectivity matrix for calcite structure CaCO3"""

        calcite_structure = Structure.from_file('test_structures/calcite.cif', True, True)
        cn_finder = EffectiveCoordFinder(calcite_structure)
        cns = cn_finder.get_avg_cn(radius=2.5, anions=['O'])
        for cation in cns:
            self.assertTrue(cation == "Ca" or cation == "C", "Ca and C should be the only ions in CaCO3")
            if cation == 'Ca':
                self.assertEqual(round(cns[cation]), 6, "K should be 8-fold coordinated")
            if cation == 'C':
                self.assertEqual(round(cns[cation]), 3, "S should be 4-fold coordinated")

        central_species = ['Ca', 'C']
        peripheral_species = ['O']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(calcite_structure, False, 2.6, peripheral_species, central_species)

        self.assertIn('Ca', connectivity_matrix.keys(), "Ca polyhedra not found in CaCO3 calcite matrix")
        self.assertIn('C', connectivity_matrix.keys(), "CO3 polyanions not found in CaCO3 calcite matrix")
        self.assertEqual(connectivity_matrix['Ca']['Ca']['point'], 6, "Ca should be point-sharing")
        self.assertEqual(connectivity_matrix['Ca']['Ca']['edge'], 0, "Ca should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Ca']['Ca']['face'], 0, "Ca should not be face-sharing")
        self.assertEqual(connectivity_matrix['C']['C']['point'], 0, "C should be isolated")
        self.assertEqual(connectivity_matrix['C']['C']['edge'], 0, "C should be isolated")
        self.assertEqual(connectivity_matrix['C']['C']['face'], 0, "C should be isolated")
        self.assertEqual(connectivity_matrix['Ca']['C']['point'], 6, "Ca and C should be point-sharing")
        self.assertEqual(connectivity_matrix['Ca']['C']['edge'], 0, "Ca and C should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Ca']['C']['face'], 0, "Ca and C should not be face-sharing")

        for connectivity_type in connectivity_matrix['Ca']['C'].keys():
            self.assertEqual(connectivity_matrix['Ca']['C'][connectivity_type]*calcite_structure.composition['Ca'],
                             connectivity_matrix['C']['Ca'][connectivity_type]*calcite_structure.composition['C'],
                             "Total number of sharing instances between 'Ca' and 'C' "
                             "should be same as total number of sharing"
                             "instances between 'C' and 'Ca'")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_corundum(self):
        """Testing coordination and connectivity matrix for corundum structure Cr2O3"""

        corundum_structure = Structure.from_file('test_structures/corundum.cif', True, True)

        cn_finder = EffectiveCoordFinder(corundum_structure)
        cns = cn_finder.get_avg_cn(radius=2.6, anions=['O'])
        for cation in cns:
            self.assertTrue(cation == "Cr", "Cr should be the only ions in corundum Cr2O3")
            if cation == 'Cr':
                self.assertEqual(round(cns[cation]), 6, "Cr should be 6-fold coordinated")

        central_species = ['Cr']
        peripheral_species = ['O']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(corundum_structure, False, 2.6, peripheral_species, central_species)

        self.assertIn('Cr', connectivity_matrix.keys(), "Cr polyhedra not found in corundum Cr2O3 matrix")
        self.assertEqual(connectivity_matrix['Cr']['Cr']['point'], 9, "Cr should be point-sharing")
        self.assertEqual(connectivity_matrix['Cr']['Cr']['edge'], 3, "Cr should be edge-sharing")
        self.assertEqual(connectivity_matrix['Cr']['Cr']['face'], 1, "Cr should be face-sharing")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_hexagonal(self):
        """Testing coordination and connectivity matrix for hexagonal structure BaNiO3"""

        hexagonal_structure = Structure.from_file('test_structures/HexagonalABX3.cif', True, True)

        cn_finder = EffectiveCoordFinder(hexagonal_structure)
        cns = cn_finder.get_avg_cn(radius=3.2, anions=['O'])
        for cation in cns:
            self.assertTrue(cation == "Ba" or cation == "Ni",
                            "Ba and Ni should be the only ions in hexagonal BaNiO3")
            if cation == 'Ba':
                self.assertEqual(round(cns[cation]), 12,
                                 "Ba should be 12-fold coordinated (effective coordination underestimates)")
            if cation == 'Ni':
                self.assertEqual(round(cns[cation]), 6, "Si should be 6-fold coordinated")

        central_species = ['Ba', 'Ni']
        peripheral_species = ['O']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(hexagonal_structure, False, 3.2, peripheral_species, central_species)

        self.assertIn('Ba', connectivity_matrix.keys(), "Ba cations not found in matrix")
        self.assertIn('Ni', connectivity_matrix.keys(), "Ni cations not found in matrix")
        self.assertEqual(connectivity_matrix['Ba']['Ba']['point'], 6, "Ba should be point-sharing")
        self.assertEqual(connectivity_matrix['Ba']['Ba']['edge'], 0, "Ba should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Ba']['Ba']['face'], 8, "Ba should be face-sharing")
        self.assertEqual(connectivity_matrix['Ni']['Ni']['point'], 0, "Ni should not be point-sharing")
        self.assertEqual(connectivity_matrix['Ni']['Ni']['edge'], 0, "Ni should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Ni']['Ni']['face'], 2, "Ni should only be face-sharing")
        self.assertEqual(connectivity_matrix['Ni']['Ba']['point'], 6, "Ni and Ba should be point-sharing")
        self.assertEqual(connectivity_matrix['Ni']['Ba']['edge'], 0, "Ni and Ba should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Ni']['Ba']['face'], 6, "Ni and Ba should be face-sharing")

        for connectivity_type in connectivity_matrix['Ba']['Ni'].keys():
            self.assertEqual(connectivity_matrix['Ba']['Ni'][connectivity_type]*hexagonal_structure.composition['Ba'],
                             connectivity_matrix['Ni']['Ba'][connectivity_type]*hexagonal_structure.composition['Ni'],
                             "Total number of sharing instances between 'Ba' and 'Ni' "
                             "should be same as total number of sharing"
                             "instances between 'Ni' and 'Ba'")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_k2nif4(self):
        """Testing coordination and connectivity matrix for K2NiF4 structure"""

        k2nif4_structure = Structure.from_file('test_structures/K2NiF4.cif', True, True)

        cn_finder = EffectiveCoordFinder(k2nif4_structure)
        cns = cn_finder.get_avg_cn(radius=3.0, anions=['F'])
        for cation in cns:
            self.assertTrue(cation == "K" or cation == "Ni",
                            "K and Ni should be the only cations in K2NiF4 structure")
            if cation == 'K':
                self.assertLessEqual(round(cns[cation]), 9,
                                     "Zr should be 9-fold coordinated (effective coordination underestimates)")
            if cation == 'Ni':
                self.assertEqual(round(cns[cation]), 6, "Ni should be 6-fold coordinated")

        central_species = ['K', 'Ni']
        peripheral_species = ['F']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(k2nif4_structure, False, 3.0, peripheral_species, central_species)

        self.assertIn('K', connectivity_matrix.keys(), "K cations not found in matrix")
        self.assertIn('Ni', connectivity_matrix.keys(), "Ni cations not found in matrix")
        self.assertEqual(connectivity_matrix['K']['K']['point'], 8, "K should be point-sharing")
        self.assertEqual(connectivity_matrix['K']['K']['edge'], 4, "K should not be edge-sharing")
        self.assertEqual(connectivity_matrix['K']['K']['face'], 5, "K should not be face-sharing")
        self.assertEqual(connectivity_matrix['Ni']['Ni']['point'], 4, "Ni should be point-sharing")
        self.assertEqual(connectivity_matrix['Ni']['Ni']['edge'], 0, "Ni should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Ni']['Ni']['face'], 0, "Ni should not be face-sharing")
        self.assertEqual(connectivity_matrix['Ni']['K']['point'], 2, "Ni and K should not be point-sharing")
        self.assertEqual(connectivity_matrix['Ni']['K']['edge'], 0, "Ni and K should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Ni']['K']['face'], 8, "Ni and K should only be face-sharing")

        for connectivity_type in connectivity_matrix['K']['Ni'].keys():
            self.assertEqual(connectivity_matrix['K']['Ni'][connectivity_type]*k2nif4_structure.composition['K'],
                             connectivity_matrix['Ni']['K'][connectivity_type]*k2nif4_structure.composition['Ni'],
                             "Total number of sharing instances between 'K' and 'Ni' "
                             "should be same as total number of sharing"
                             "instances between 'Ni' and 'K'")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_olivine(self):
        """Testing coordination and connectivity matrix for olivine structure Fe2SiO4"""

        olivine_structure = Structure.from_file('test_structures/olivine.cif', True, True)

        cn_finder = EffectiveCoordFinder(olivine_structure)
        cns = cn_finder.get_avg_cn(radius=2.6, anions=['O'])
        for cation in cns:
            self.assertTrue(cation == "Fe" or cation == "Si",
                            "Fe and Si cations should be the only cations in olivine Fe2SiO4")
            if cation == 'Fe':
                self.assertEqual(round(cns[cation]), 6, "Fe should be 6-fold coordinated")
            if cation == 'Si':
                self.assertEqual(round(cns[cation]), 4, "Si should be 4-fold coordinated")

        central_species = ['Fe', 'Si']
        peripheral_species = ['O']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(olivine_structure, False, 2.6, peripheral_species, central_species)

        self.assertIn('Fe', connectivity_matrix.keys(), "Fe polyhedra not found in matrix")
        self.assertIn('Si', connectivity_matrix.keys(), "Si polyhedra not found in matrix")
        self.assertEqual(connectivity_matrix['Fe']['Fe']['point'], 6, "Fe should be point-sharing")
        self.assertEqual(connectivity_matrix['Fe']['Fe']['edge'], 3, "Fe should be edge-sharing")
        self.assertEqual(connectivity_matrix['Fe']['Fe']['face'], 0, "Fe should be face-sharing")
        self.assertEqual(connectivity_matrix['Si']['Si']['point'], 0, "Si should be isolated")
        self.assertEqual(connectivity_matrix['Si']['Si']['edge'], 0, "Si should be isolated")
        self.assertEqual(connectivity_matrix['Si']['Si']['face'], 0, "Si should be isolated")
        self.assertEqual(connectivity_matrix['Si']['Fe']['point'], 6, "Fe and Si should be point-sharing")
        self.assertEqual(connectivity_matrix['Si']['Fe']['edge'], 3, "Fe and Si should be edge-sharing")
        self.assertEqual(connectivity_matrix['Si']['Fe']['face'], 0, "Fe and Si should not be face-sharing")

        for connectivity_type in connectivity_matrix['Fe']['Si'].keys():
            self.assertEqual(connectivity_matrix['Fe']['Si'][connectivity_type]*olivine_structure.composition['K'],
                             connectivity_matrix['Si']['Fe'][connectivity_type]*olivine_structure.composition['S'],
                             "Total number of sharing instances between 'Fe' and 'Si' "
                             "should be same as total number of sharing"
                             "instances between 'Si' and 'Fe'")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_phenacite(self):
        """Testing coordination and connectivity matrix for phenacite structure Be2SiO4"""

        phenacite_structure = Structure.from_file('test_structures/phenacite.cif', True, True)
        cn_finder = EffectiveCoordFinder(phenacite_structure)
        cns = cn_finder.get_avg_cn(radius=2.6, anions=['O'])
        for cation in cns:
            self.assertTrue(cation == "Be" or cation == "Si",
                            "Be and Si polyanions should be the only ions in phenacite Be2SiO4")
            if cation == 'Be':
                self.assertEqual(round(cns[cation]), 4, "Be should be 4-fold coordinated")
            if cation == 'Si':
                self.assertEqual(round(cns[cation]), 4, "Si should be 4-fold coordinated")

        central_species = ['Be', 'Si']
        peripheral_species = ['O']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(phenacite_structure, False, 2.6, peripheral_species, central_species)
        self.assertIn('Be', connectivity_matrix.keys(), "Be cations not found in matrix")
        self.assertIn('Si', connectivity_matrix.keys(), "Si cations not found in matrix")
        self.assertEqual(connectivity_matrix['Be']['Be']['point'], 4, "Be should be point-sharing")
        self.assertEqual(connectivity_matrix['Be']['Be']['edge'], 0, "Be should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Be']['Be']['face'], 0, "Be should not be face-sharing")
        self.assertEqual(connectivity_matrix['Si']['Si']['point'], 0, "Si should be isolated")
        self.assertEqual(connectivity_matrix['Si']['Si']['edge'], 0, "Si should be isolated")
        self.assertEqual(connectivity_matrix['Si']['Si']['face'], 0, "Si should be isolated")
        self.assertEqual(connectivity_matrix['Be']['Si']['point'], 4, "Be and S should be point-sharing")
        self.assertEqual(connectivity_matrix['Be']['Si']['edge'], 0, "Be and S should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Be']['Si']['face'], 0, "Be and S should not be face-sharing")

        for connectivity_type in connectivity_matrix['Be']['Si'].keys():
            self.assertEqual(connectivity_matrix['Be']['Si'][connectivity_type]*phenacite_structure.composition['Be'],
                             connectivity_matrix['Si']['Be'][connectivity_type]*phenacite_structure.composition['Si'],
                             "Total number of sharing instances between 'Be' and 'Si' "
                             "should be same as total number of sharing"
                             "instances between 'Si' and 'Be'")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_pyroxene(self):
        """Testing coordination and connectivity matrix for pyroxene structure Mg16Si16O48"""

        pyroxene_structure = Structure.from_file('test_structures/pyroxene.cif', True, True)

        cn_finder = EffectiveCoordFinder(pyroxene_structure)
        cns = cn_finder.get_avg_cn(radius=2.6, anions=['O'])
        for cation in cns:
            self.assertTrue(cation == "Mg" or cation == "Si", "Mg and Si should be the only ions in pyroxene MgSiO4")
            if cation == 'Mg':
                self.assertLessEqual(round(cns[cation]), 6,
                                     "Mg should be 6-fold coordinated (effective coordination underestimates)")
            if cation == 'Si':
                self.assertEqual(round(cns[cation]), 4, "Si should be 4-fold coordinated")

        central_species = ['Mg', 'Si']
        peripheral_species = ['O']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(pyroxene_structure, False, 2.6, peripheral_species, central_species)

        self.assertIn('Mg', connectivity_matrix.keys(), "Mg cations not found in matrix")
        self.assertIn('Si', connectivity_matrix.keys(), "Si cations not found in matrix")
        self.assertEqual(connectivity_matrix['Mg']['Mg']['point'], 0, "Mg should be point-sharing")
        self.assertEqual(connectivity_matrix['Mg']['Mg']['edge'], 4, "Mg should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Mg']['Mg']['face'], 0, "Mg should not be face-sharing")
        self.assertEqual(connectivity_matrix['Si']['Si']['point'], 2, "Ti should be point-sharing")
        self.assertEqual(connectivity_matrix['Si']['Si']['edge'], 0, "Ti should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Si']['Si']['face'], 0, "Ti should not be face-sharing")
        self.assertEqual(connectivity_matrix['Si']['Mg']['point'], 6, "Ti and Mg should not be point-sharing")
        self.assertEqual(connectivity_matrix['Si']['Mg']['edge'], 0.5, "Ti and Mg should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Si']['Mg']['face'], 0, "Ti and Mg should only be face-sharing")

        for connectivity_type in connectivity_matrix['Mg']['Si'].keys():
            self.assertEqual(connectivity_matrix['Mg']['Si'][connectivity_type]*pyroxene_structure.composition['Mg'],
                             connectivity_matrix['Si']['Mg'][connectivity_type]*pyroxene_structure.composition['Si'],
                             "Total number of sharing instances between 'Mg' and 'Si' "
                             "should be same as total number of sharing"
                             "instances between 'Si' and 'Mg'")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_reo3(self):
        """Testing coordination and connectivity matrix for ReO3 structure"""

        reo3_structure = Structure.from_file('test_structures/ReO3.cif', True, True)

        cn_finder = EffectiveCoordFinder(reo3_structure)
        cns = cn_finder.get_avg_cn(radius=2.6, anions=['O'])
        for cation in cns:
            self.assertTrue(cation == "Re", "Re should be the only ions in ReO3 structure")
            if cation == 'Re':
                self.assertEqual(round(cns[cation]), 6, "Re should be 6-fold coordinated")

        central_species = ['Re']
        peripheral_species = ['O']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(reo3_structure, False, 2.6, peripheral_species, central_species)
        self.assertIn('Re', connectivity_matrix.keys(), "Re polyhedra not found in ReO3 matrix")
        self.assertEqual(connectivity_matrix['Re']['Re']['point'], 6, "Re should be point-sharing")
        self.assertEqual(connectivity_matrix['Re']['Re']['edge'], 0, "Re should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Re']['Re']['face'], 0, "Re should not be face-sharing")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_rocksalt(self):
        """Testing coordination and connectivity matrix for rock salt structure NaCl"""

        rocksalt_structure = Structure.from_file('test_structures/rocksalt.cif', True, True)

        cn_finder = EffectiveCoordFinder(rocksalt_structure)
        cns = cn_finder.get_avg_cn(radius=3.0, anions=['Cl'])
        for cation in cns:
            self.assertTrue(cation == "Na", "Na should be the only ions in rocksalt NaCl structure")
            if cation == 'Na':
                self.assertEqual(round(cns[cation]), 6, "Na should be 6-fold coordinated")

        central_species = ['Na']
        peripheral_species = ['Cl']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(rocksalt_structure, False, 3.0, peripheral_species, central_species)
        self.assertIn('Na', connectivity_matrix.keys(), "Na polyhedra not found in NaCl matrix")
        self.assertEqual(connectivity_matrix['Na']['Na']['point'], 6, "Na should not be point-sharing")
        self.assertEqual(connectivity_matrix['Na']['Na']['edge'], 12, "Na should be edge-sharing")
        self.assertEqual(connectivity_matrix['Na']['Na']['face'], 0, "Na should not be face-sharing")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_rutile(self):
        """Testing coordination and connectivity matrix for rutile structure TiO2"""

        rutile_structure = Structure.from_file('test_structures/rutile.cif', True, True)

        cn_finder = EffectiveCoordFinder(rutile_structure)
        cns = cn_finder.get_avg_cn(radius=2.6, anions=['O'])
        for cation in cns:
            self.assertTrue(cation == "Ti", "Ti should be the only ions in rocksalt NaCl structure")
            if cation == 'Ti':
                self.assertEqual(round(cns[cation]), 6, "Ti should be 6-fold coordinated")

        central_species = ['Ti']
        peripheral_species = ['O']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(rutile_structure, False, 2.6, peripheral_species, central_species)
        self.assertIn('Ti', connectivity_matrix.keys(), "Ti polyhedra not found in NaCl matrix")
        self.assertEqual(connectivity_matrix['Ti']['Ti']['point'], 4, "Ti should not be point-sharing")
        self.assertEqual(connectivity_matrix['Ti']['Ti']['edge'], 4, "Ti should be edge-sharing")
        self.assertEqual(connectivity_matrix['Ti']['Ti']['face'], 0, "Ti should not be face-sharing")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_scheelite(self):
        """Testing coordination and connectivity matrix for scheelite structure CaWO4"""

        scheelite_structure = Structure.from_file('test_structures/scheelite.cif', True, True)

        cn_finder = EffectiveCoordFinder(scheelite_structure)
        cns = cn_finder.get_avg_cn(radius=2.6, anions=['O'])
        for cation in cns:
            self.assertTrue(cation == "Ca" or cation == "W", "Ca and W should be the only ions in scheelite CaWO4")
            if cation == 'Ca':
                self.assertEqual(round(cns[cation]), 8,
                                 "Ca should be 8-fold coordinated (effective coordination underestimates)")
            if cation == 'W':
                self.assertEqual(round(cns[cation]), 4, "W should be 4-fold coordinated")

        central_species = ['Ca', 'W']
        peripheral_species = ['O']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(scheelite_structure, False, 2.6, peripheral_species, central_species)

        self.assertIn('Ca', connectivity_matrix.keys(), "Ca cations not found in matrix")
        self.assertIn('W', connectivity_matrix.keys(), "W cations not found in matrix")
        self.assertEqual(connectivity_matrix['Ca']['Ca']['point'], 0, "Ca should not be point-sharing")
        self.assertEqual(connectivity_matrix['Ca']['Ca']['edge'], 4, "Ca should be be edge-sharing")
        self.assertEqual(connectivity_matrix['Ca']['Ca']['face'], 0, "Ca should not be face-sharing")
        self.assertEqual(connectivity_matrix['W']['W']['point'], 0, "W should not be point-sharing")
        self.assertEqual(connectivity_matrix['W']['W']['edge'], 0, "W should not be edge-sharing")
        self.assertEqual(connectivity_matrix['W']['W']['face'], 0, "W should not be face-sharing")
        self.assertEqual(connectivity_matrix['Ca']['W']['point'], 8, "Ca and W should only be point-sharing")
        self.assertEqual(connectivity_matrix['Ca']['W']['edge'], 0, "Ca and W should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Ca']['W']['face'], 0, "Ca and W should not be face-sharing")

        for connectivity_type in connectivity_matrix['Ca']['W'].keys():
            self.assertEqual(connectivity_matrix['Ca']['W'][connectivity_type]*scheelite_structure.composition['Ca'],
                             connectivity_matrix['W']['Ca'][connectivity_type]*scheelite_structure.composition['W'],
                             "Total number of sharing instances between 'Ca' and 'W' "
                             "should be same as total number of sharing"
                             "instances between 'W' and 'Ca'")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_sio2(self):
        """Testing coordination and connectivity matrix for SiO2 structure"""

        sio2_structure = Structure.from_file('test_structures/SiO2.cif', True, True)

        cn_finder = EffectiveCoordFinder(sio2_structure)
        cns = cn_finder.get_avg_cn(radius=2.6, anions=['O'])
        for cation in cns:
            self.assertTrue(cation == "Si", "Si should be the only ions in SiO2 structure")
            if cation == 'Si':
                self.assertEqual(round(cns[cation]), 4, "Si should be 6-fold coordinated")

        central_species = ['Si']
        peripheral_species = ['O']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(sio2_structure, False, 2.6, peripheral_species, central_species)
        self.assertIn('Si', connectivity_matrix.keys(), "Si polyhedra not found in SiO2 matrix")
        self.assertEqual(connectivity_matrix['Si']['Si']['point'], 4, "Si should not be point-sharing")
        self.assertEqual(connectivity_matrix['Si']['Si']['edge'], 0, "Si should be edge-sharing")
        self.assertEqual(connectivity_matrix['Si']['Si']['face'], 0, "Si should not be face-sharing")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_spinel(self):
        """Testing coordination and connectivity matrix for spinel structure MgAl2O4"""

        spinel_structure = Structure.from_file('test_structures/spinel.cif', True, True)
        cn_finder = EffectiveCoordFinder(spinel_structure)
        cns = cn_finder.get_avg_cn(radius=2.6, anions=['O'])
        for cation in cns:
            self.assertTrue(cation == "Mg" or cation == "Al",
                            "Mg and Al polyanions should be the only ions in spinel MgAl2O4")
            if cation == 'Mg':
                self.assertEqual(round(cns[cation]), 4, "Mg should be 4-fold coordinated")
            if cation == 'Al':
                self.assertEqual(round(cns[cation]), 6, "Al should be 4-fold coordinated")

        central_species = ['Mg', 'Al']
        peripheral_species = ['O']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(spinel_structure, False, 2.6, peripheral_species, central_species)
        self.assertIn('Mg', connectivity_matrix.keys(), "Mg cations not found in matrix")
        self.assertIn('Al', connectivity_matrix.keys(), "Al cations not found in matrix")

        self.assertEqual(connectivity_matrix['Mg']['Mg']['point'], 0, "Mg should be point-sharing")
        self.assertEqual(connectivity_matrix['Mg']['Mg']['edge'], 0, "Mg should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Mg']['Mg']['face'], 0, "Mg should not be face-sharing")
        self.assertEqual(connectivity_matrix['Al']['Al']['point'], 0, "Al should be isolated")
        self.assertEqual(connectivity_matrix['Al']['Al']['edge'], 6, "Al should be isolated")
        self.assertEqual(connectivity_matrix['Al']['Al']['face'], 0, "Al should be isolated")
        self.assertEqual(connectivity_matrix['Mg']['Al']['point'], 12, "Mg and Al should be point-sharing")
        self.assertEqual(connectivity_matrix['Mg']['Al']['edge'], 0, "Mg and Al should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Mg']['Al']['face'], 0, "Mg and Al should not be face-sharing")

        for connectivity_type in connectivity_matrix['Mg']['Al'].keys():
            self.assertEqual(connectivity_matrix['Mg']['Al'][connectivity_type]*spinel_structure.composition['Mg'],
                             connectivity_matrix['Al']['Mg'][connectivity_type]*spinel_structure.composition['Al'],
                             "Total number of sharing instances between 'Mg' and 'Al' "
                             "should be same as total number of sharing"
                             "instances between 'Al' and 'Mg'")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_perovskite(self):
        """Testing coordination and connectivity matrix for perovskite structure SrTiO3"""

        perovskite_structure = Structure.from_file('test_structures/perovskite.cif', True, True)
        cn_finder = EffectiveCoordFinder(perovskite_structure)
        cns = cn_finder.get_avg_cn(radius=3.2, anions=['O'])
        for cation in cns:
            self.assertTrue(cation == "Sr" or cation == "Ti",
                            "Be and Si polyanions should be the only ions in phenacite Be2SiO4")
            if cation == 'Sr':
                self.assertLessEqual(round(cns[cation]), 12,
                                     "Sr should be 12-fold coordinated (effective coordination underestimates)")
            if cation == 'Ti':
                self.assertEqual(round(cns[cation]), 6, "Ti should be 6-fold coordinated")

        central_species = ['Sr', 'Ti']
        peripheral_species = ['O']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(perovskite_structure, False, 3.2, peripheral_species, central_species)
        self.assertIn('Sr', connectivity_matrix.keys(), "Sr cations not found in matrix")
        self.assertIn('Ti', connectivity_matrix.keys(), "Ti cations not found in matrix")
        self.assertEqual(connectivity_matrix['Sr']['Sr']['point'], 12, "Sr should be point-sharing")
        self.assertEqual(connectivity_matrix['Sr']['Sr']['edge'], 0, "Sr should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Sr']['Sr']['face'], 6, "Sr should not be face-sharing")
        self.assertEqual(connectivity_matrix['Ti']['Ti']['point'], 6, "Ti should be point-sharing")
        self.assertEqual(connectivity_matrix['Ti']['Ti']['edge'], 0, "Ti should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Ti']['Ti']['face'], 0, "Ti should not be face-sharing")
        self.assertEqual(connectivity_matrix['Ti']['Sr']['point'], 0, "Sr and Ti should not be point-sharing")
        self.assertEqual(connectivity_matrix['Ti']['Sr']['edge'], 0, "Sr and Ti should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Ti']['Sr']['face'], 8, "Sr and Ti should only be face-sharing")

        for connectivity_type in connectivity_matrix['Sr']['Ti'].keys():
            self.assertEqual(connectivity_matrix['Sr']['Ti'][connectivity_type]*perovskite_structure.composition['Sr'],
                             connectivity_matrix['Ti']['Sr'][connectivity_type]*perovskite_structure.composition['Ti'],
                             "Total number of sharing instances between 'Sr' and 'Ti' "
                             "should be same as total number of sharing"
                             "instances between 'Ti' and 'Sr'")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_zircon(self):
        """Testing coordination and connectivity matrix for zircon structure ZrSiO4"""

        zircon_structure = Structure.from_file('test_structures/zircon.cif', True, True)
        cn_finder = EffectiveCoordFinder(zircon_structure)
        cns = cn_finder.get_avg_cn(radius=2.6, anions=['O'])
        for cation in cns:
            self.assertTrue(cation == "Zr" or cation == "Si",
                            "Zr and Si polyanions should be the only ions in zircon ZrSiO4")
            if cation == 'Zr':
                self.assertLessEqual(round(cns[cation]), 8,
                                     "Zr should be 8-fold coordinated (effective coordination underestimates)")
            if cation == 'Si':
                self.assertEqual(round(cns[cation]), 4, "Si should be 4-fold coordinated")

        central_species = ['Zr', 'Si']
        peripheral_species = ['O']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(zircon_structure, False, 2.6, peripheral_species, central_species)
        self.assertIn('Zr', connectivity_matrix.keys(), "Zr cations not found in matrix")
        self.assertIn('Si', connectivity_matrix.keys(), "Si cations not found in matrix")
        self.assertEqual(connectivity_matrix['Zr']['Zr']['point'], 0, "Zr should not be point-sharing")
        self.assertEqual(connectivity_matrix['Zr']['Zr']['edge'], 4, "Zr should be edge-sharing")
        self.assertEqual(connectivity_matrix['Zr']['Zr']['face'], 0, "Zr should not be face-sharing")
        self.assertEqual(connectivity_matrix['Si']['Si']['point'], 0, "Si should be point-sharing")
        self.assertEqual(connectivity_matrix['Si']['Si']['edge'], 0, "Si should not be edge-sharing")
        self.assertEqual(connectivity_matrix['Si']['Si']['face'], 0, "Si should not be face-sharing")
        self.assertEqual(connectivity_matrix['Si']['Zr']['point'], 4, "Zr and Si should be point-sharing")
        self.assertEqual(connectivity_matrix['Si']['Zr']['edge'], 2, "Zr and Si should be edge-sharing")
        self.assertEqual(connectivity_matrix['Si']['Zr']['face'], 0, "Zr and Si should not be face-sharing")

        for connectivity_type in connectivity_matrix['Zr']['Si'].keys():
            self.assertEqual(connectivity_matrix['Zr']['Si'][connectivity_type]*zircon_structure.composition['Zr'],
                             connectivity_matrix['Si']['Zr'][connectivity_type]*zircon_structure.composition['Si'],
                             "Total number of sharing instances between 'Zr' and 'Si' "
                             "should be same as total number of sharing"
                             "instances between 'Si' and 'Zr'")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")

    def test_znso4(self):
        """Testing coordination and connectivity matrix for ZnSO4 structure"""

        znso4_structure = Structure.from_file('test_structures/ZnSO4.cif', True, True)
        cn_finder = EffectiveCoordFinder(znso4_structure)
        cns = cn_finder.get_avg_cn(radius=2.6, anions=['O'])
        for cation in cns:
            self.assertTrue(cation == "Zn" or cation == "S", "Zn and S should be the only ions in zircon ZrSiO4")
            if cation == 'Zn':
                self.assertLessEqual(round(cns[cation]), 6,
                                     "Zn should be 6-fold coordinated (effective coordination underestimates)")
            if cation == 'S':
                self.assertEqual(round(cns[cation]), 4, "S should be 4-fold coordinated")

        central_species = ['Zn', 'S']
        peripheral_species = ['O']

        connectivity_matrix, connectivity_polyhedra = \
            get_connectivity_matrix(znso4_structure, False, 2.6, peripheral_species, central_species)
        self.assertIn('Zn', connectivity_matrix.keys(), "Zn cations not found in matrix")
        self.assertIn('S', connectivity_matrix.keys(), "S cations not found in matrix")
        self.assertEqual(connectivity_matrix['Zn']['Zn']['point'], 0, "Zn should not be point-sharing")
        self.assertEqual(connectivity_matrix['Zn']['Zn']['edge'], 2, "Zn should be edge-sharing")
        self.assertEqual(connectivity_matrix['Zn']['Zn']['face'], 0, "Zn should not be face-sharing")
        self.assertEqual(connectivity_matrix['S']['S']['point'], 0, "S should be isolated")
        self.assertEqual(connectivity_matrix['S']['S']['edge'], 0, "S should be isolated")
        self.assertEqual(connectivity_matrix['S']['S']['face'], 0, "S should be isolated")
        self.assertEqual(connectivity_matrix['S']['Zn']['point'], 6, "S and Zn should only be point-sharing")
        self.assertEqual(connectivity_matrix['S']['Zn']['edge'], 0, "S and Zn should not be edge-sharing")
        self.assertEqual(connectivity_matrix['S']['Zn']['face'], 0, "S and Zn should be face-sharing")

        for connectivity_type in connectivity_matrix['Zn']['S'].keys():
            self.assertEqual(connectivity_matrix['Zn']['S'][connectivity_type]*znso4_structure.composition['Zn'],
                             connectivity_matrix['S']['Zn'][connectivity_type]*znso4_structure.composition['S'],
                             "Total number of sharing instances between 'Zn' and 'S' "
                             "should be same as total number of sharing"
                             "instances between 'S' and 'Zn'")

        for poly in connectivity_polyhedra:
            self.assertIsInstance(poly, Polyhedra, "List of polyhedra includes a non-polyhedra element")


if __name__ == '__main__':
    unittest.main()