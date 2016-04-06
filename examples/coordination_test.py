from pymatgen.core.structure import Structure
from structural_descriptors_repo import effective_coordination as ECoN
from structural_descriptors_repo import okeeffe_coordination_number_from_structure as OKeeffe
import unittest

# TODO: make real unit tests!!
class TestCoordinationMethods(unittest.TestCase):

    def setUp(self):
        self.assertEqual(1,1)

    def test_1(self):
        self.assertEqual(1, 1)

if __name__ == '__main__':
    s = Structure.from_file('LiCoO2.cif', True, True)

    print "Example structure: LiCoO2"
    print s
    print ""

    print "Average effective coordination number given by 3.2 angstrom radius for each cation"
    finder = ECoN.EffectiveCoordFinder(s)
    print finder.getAvgCN(3.2)
    print ""

    print "Average O'Keeffe coordination number for each cation"
    print OKeeffe.getAvgCN(s)
    print ""