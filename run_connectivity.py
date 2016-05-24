from pymatgen.core.structure import Structure
from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
import connectivity_from_structure as connectivity

__author__ = 'Tina_Chen'
__contributor__ = 'Anubhav Jain'

if __name__ == '__main__':

    # path to structure file
    path = 'LiCoO2.cif'

    # optional list of strings with central/target species abbreviations
    # (i.e. we aim to find connectivity between polyhedra with these species at the center)
    central_species = ['Li']

    # optional list of strings with peripheral species abbreviations
    # (i.e. we aim to find connectivity where these species are surrounding the central/target species)
    peripheral_species = None

    # radius from the central/target sites to which we look for peripheral species
    # only needs to be changed if central species is very large (eg: Ba, Sr)
    radius = 2.6

    structure = Structure.from_file(path, True, True)
    print structure
    connectivity_matrix, all_polyhedra, supercell = \
        connectivity.get_connectivity_matrix(structure, True, radius, peripheral_species, central_species)
    descriptions = connectivity.get_connectivity_description(connectivity_matrix, all_polyhedra, structure, True, anions=peripheral_species)

    for cation in descriptions:
        print ""
        print descriptions[cation]
