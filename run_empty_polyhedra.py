from pymatgen import Structure
from Voronoi_sites import add_voronoi
from empty_polyhedra_util import split_sites, structure_from_sites
from empty_polyhedra import step1, step2, step3


__author__ = 'Tina'


class EmptyPolyhedraFinder():

    def __init__(self, input_structure, voronoi_species='Rn', max_radius=2.6, min_neighbors=4, anions_list=None,
                 num_matches=4):
        self._structure = input_structure
        self._voronoi_species = voronoi_species
        self._radius = max_radius
        self._min_neighbors = min_neighbors
        if anions_list is None:
            self._anions = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S', 'N', 'N3-']
        else:
            self._anions = anions_list
        self._num_matches = num_matches

    def run_empty_polyhedra(self):

        dup_structure = self._structure.copy()

        structure_with_voronoi = add_voronoi(dup_structure)

        # find all voronoi sites in the structure
        print "splitting sites into voronoi and non-voronoi"
        allvoronois, othersites = split_sites(structure_with_voronoi, self._voronoi_species)

        # process potential empty polyhedral center sites
        print "processing structure, step 1"
        remaining_voronois = step1(allvoronois, structure_with_voronoi, self._radius, self._anions,
                                   self._min_neighbors)

        # create intermediate structure
        print "creating structure1"
        structure_with_less_voronoi = structure_from_sites(remaining_voronois, othersites, structure_with_voronoi)

        print "processing structure, step 2 (may take a few minutes)"
        # find polyhedral of voronoi, organizing them by pairing those with the same set of peripheral anions
        # match each voronoi point in voronois to a number - dictionary
        # create a method like getCationPolyhedral that finds the polyhedral for the specific site
        # check polyhedral in the same way
        # use the number from 2nd step to choose which sites to remove
        remaining_voronois = step2(remaining_voronois, structure_with_less_voronoi, self._radius, self._anions,
                                   self._voronoi_species, self._min_neighbors)

        print "creating structure3"
        structure_with_less_voronoi = structure_from_sites(remaining_voronois, othersites, structure_with_less_voronoi)

        print "processing structure, step 3"
        # add criterion based on C-V and C-O distances
        # if C-O < 2.5, if C-V distance < C-O distance, then remove
        # for cation closest to V with all nearby oxygen
        final_voronois = step3(remaining_voronois, structure_with_less_voronoi, self._radius, self._anions,
                               self._num_matches)


        print "creating structure4"
        structure_with_final_voronoi = structure_from_sites(final_voronois, othersites, structure_with_less_voronoi)

        return structure_with_final_voronoi

if __name__ == '__main__':
    path = "CaF2"

    structure = Structure.from_file("%s.cif" % path, True, False)
    epf = EmptyPolyhedraFinder(structure)
    final_structure = epf.run_empty_polyhedra()
    final_structure.to(filename="%s_with_voronoi.cif" % path)

    # Add one more step which is if the Voronoi points are within 1 angstrom of each other, then remove one