__author__ = 'Tina_Chen'
__contributor__ = 'Stephen_Dacek'


from pymatgen import Structure
from pymatgen import PeriodicSite
import effective_coordination
from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
from itertools import combinations
import math
import numpy


"""
Finds the connectivities between polyhedra within a structure, where we define connectivity by the number of shared
peripheral ions between two polyhedra. 0 shared peripheral ions indicates that there is no connectivity between the
two polyhedra; 1 shared peripheral ion indicates that there is point-sharing connectivity between the polyhedra; 2
shared peripheral ions indicate that there is edge-sharing connectivity between the polyhedra; 3 or more shared
peripheral ions indicate that there is face-sharing connectivity between the polyhedra.

For the purpose of this code, the peripheral ions will be defined by the ions surrounding a central site as given by
structure's get_neighbors function, with bond weight contributions (given by Hoppe, 1979) less than a set cut-off value.

"""

class Polyhedra(object):

    """
    Object representing a polyhedra in a structure with central ion site "cation" and surrounding ions site
    "peripheralIons"

    """

    def __init__(self, cation, peripheralIons):
        """

        :param cation: Site object representing the central cation
        :param peripheralIons: list of Site objects representing the central cation's surrounding ions
        """

        self.central_ion = cation
        self.central_ion_name = cation.species_string
        self.peripheral_ions = peripheralIons
        self.cation_num = 1     # Number to identify unique cation site in center cell from other sites with same species

        self.composition = self.central_ion._species
        for site in self.peripheral_ions:
            self.composition += site._species


    def set_site_number(self, cation_number):
        self.cation_num = cation_number

    def get_num_connections(self, other):
        """
        Gives the connectivity between the given polyhedra and another polyhedra by counting the number of atoms shared
        between the two polyhedra

        :param other: (Polyhedra) target Polyhedra against which we are checking for connectivity with the current
        polyhedra
        :return: (int) Integer giving the number of peripheral ions shared by the current polyhedra and the target
        polyhedra; returns -1 if the same polyhedra is being compared to itself
        """
        return len(self.get_connections(other))

    def get_connections(self, other):
        """
        Gives the shared sites between the current Polyhedra and the given Polyhedra

        :param polyhedra2: (Polyhedra) target Polyhedra against which we are checking for connectivity with the current
        polyhedra
        :return: (list) list of sites shared between the current Polyhedra and the given Polyhedra
        """

        if self == other:
            raise ValueError('Checking connections between exact same polyhedra')

        shared_sites = []

        for site in other.peripheral_ions:
            for c_site in self.peripheral_ions:
                if c_site == site:
                    shared_sites.append(site)
        return shared_sites

    def __eq__(self, other):
        """

        :param other: another Polyhedra
        :return: true is the central ion site of the current polyhedra is the same as the other polyhedra
        """
        if self.central_ion == other.central_ion:
            return True
        else:
            return False


    def __str__(self):
        """
        List of properties of the polyhedra: central species name, central species site, list of peripheral ions of
        the polyhedra
        """

        return str([self.central_ion_name,
                    [self.central_ion.frac_coords[0], self.central_ion.frac_coords[1], self.central_ion.frac_coords[2]],
                    self.peripheral_ions])


def get_surrounding_connectivity(structure, polyhedra, radius, peripheral_species):
    """
    Gives the surrounding connectivity of a specific polyhedra in a given structure

    :param structure: (Structure) target structure
    :param polyhedra: (Polyhedra) target polyhedra
    :return: (list) list of polyhedra and their connections to the target polyhedra in the target structure, elements
    of the list are of the form: [[polyhedra, connection], ...]
    """

    polyhedraList = get_supercell_polyhedra(structure, radius, peripheral_species)

    connectedPolyhedra = []
    for otherPolyhedra in polyhedraList:
        connection = polyhedra.get_num_connections(otherPolyhedra)
        if connection > 0:
            connectedPolyhedra.append([otherPolyhedra, connection])

    return connectedPolyhedra

def get_polyhedra(structure, site, peripheral_species, radius=2.8):
    """
    Gives a polyhedra in a structure based on the site, where nearby atoms are considered peripheral ions to the site
    if their bond weights (as defined by Hoppe, 1979) are greater than given value (i.e. the atoms contribute to the
    ECoN)

    :param structure: (Structure) target structure
    :param site: (PeriodicSite) target site
    :param peripheral_species: (list of Strings) List of strings with species names of the ions that can be peripheral
    ions
    :param radius: (float) largest possible distance within which peripheral ions can be obtained
    :return: (Polyhedra) polyhedra object representing the polyhedra around the target site in the target structure
    """
    #anions = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S']  # this should be a tuneable parameter

    anionSites = []
    bondlengths = []
    peripheralSites = []

    for entry in structure.get_neighbors(site, radius):
        if entry[1] < radius:
            if entry[0].species_string in peripheral_species:
                anionSites.append(entry[0])
                bondlengths.append(entry[1])
    for potentialPeripheralIons in anionSites:
        #do not count nearby ions that do not contribute
        if effective_coordination.calculate_bond_weight(potentialPeripheralIons[1], bondlengths) > 0.5:
            peripheralSites.append(potentialPeripheralIons)

    return Polyhedra(site, peripheralSites)



def count_connectivity(structure, radius = 2.8):

    """

    Gives the total counts of each type of connectivity (point-sharing, edge-sharing, face-sharing, as well as polyhedra
    that do not have sharing with any other polyhedra) within a given structure for each cation

    :param structure: (Structure) target structure
    :return: (dict) dictionary of all connectivity counts for each cation, where connectivity of 0 implies no
    connectivity, 1 implies point connectivity, 2 implies edge connectivity, and 3 or above implies face connectivity;
    the keys of the dict contain a cation while the values of the dict are the counts of each type of
    connectivity in a list, dict(cation) = [no-connectivity counts, point-connectivity counts, edge-connectivity counts,
    face-connectivity counts]
    """

    polyhedra_list = get_supercell_polyhedra(structure, radius)

    # Get the polyhedra in the central cell of the supercell from which to calculate connectivities
    center_cell_polyhedra = []
    for polyhedra in polyhedra_list:
        site_coords = polyhedra.central_ion.frac_coords
        if site_coords[0] >= (1.0/3)-(1e-5) and site_coords[0] < (2.0/3):
            if site_coords[1] >= (1.0/3)-(1e-5) and site_coords[1] < (2.0/3):
                if site_coords[2] >= (1.0/3)-(1e-5) and site_coords[2] < (2.0/3):
                    center_cell_polyhedra.append(polyhedra)


    # Count the number of connections between the
    connections = {}
    for polyhedra1 in center_cell_polyhedra:
        cation = polyhedra1.central_ion_name
        connected = False #keep track of whether a polyhedra is not connected to any other polyhedra
        if not cation in connections.keys():
            connections[cation] = {"none": 0, "point": 0,"edge": 0, "face": 0}
        for polyhedra2 in polyhedra_list:
            connection = polyhedra1.get_num_connections(polyhedra2)
            if connection < 0:
                continue
            if connection == 1:
                connections[cation]["point"] += 1
                connected = True
            if connection == 2:
                connections[cation]["edge"] += 1
                connected = True
            if connection >= 3:
                connections[cation]["face"] += 1
                connected = True
        if connected == False:
            connections[cation]["none"] += 1

    return connections

def get_supercell_polyhedra(structure, radius, peripheral_species, central_species = []):
    """
    Obtains all cation polyhedra needed to describe connectivity by creating a supercell and converting sites in the
    structure to cation polyhedra (peripheral ions of the polyhedra given by the

    :param structure: (Structure) target structure
    :param radius: (float) radius within which to determine whether a nearby atom is a peripheral ion
    :param peripheral_species: (list of Strings) List of strings with species names of the ions that can be peripheral
    ions
    :param central_species: (list of Strings) List of strings with species names of the ions that can be central ions;
    if there are no central species specified, then we consider all ions not in the list of peripheral ions to be
    central species
    :return: (list) list of Polyhedra objects in supercell large enough to determine connectivities between polyhedra
    """

    structure.make_supercell((6, 6, 6))

    polyhedra = []
    if len(central_species) == 0:
        for site_index in range(structure.num_sites):
            if structure[site_index].species_string not in peripheral_species:
                polyhedra.append(get_polyhedra(structure, structure[site_index], peripheral_species, radius))
    else:
        for site_index in range(structure.num_sites):
            if structure[site_index].species_string in central_species:
                polyhedra.append(get_polyhedra(structure, structure[site_index], peripheral_species, radius))

    return polyhedra


def get_connectivity_matrix(structure, radius = 2.8, peripheral_species=['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S'], central_species = []):
    """
    Creates a connectivity matrix to describe the connectivity between cations in a structure; connections between
    cation polyhedra and reflections of itself (and reflections of other cation polyhedra) are also counted; different
    sites of the same cation species are NOT differentiated (that is, all instances of each cation are counted all
    together)

    :param structure: (Structure) target structure
    :param radius: (float) radius within which to determine whether a nearby atom is a peripheral ion
    :param peripheral_species: (list of Strings) List of strings with species names of the ions that can be peripheral
    ions
    :param central_species: (list of Strings) List of strings with species names of the ions that can be central ions;
    if there are no central species specified, then we consider all ions not in the list of peripheral ions to be
    central species
    :return: (dict) dictionary of dictionaries, with the first set of keys being the specified (unreflected) Polyhedra
    sites, the second (inner) set of keys being the specified (unreflected) Polyhedra sites which we are comparing
    connectivity to, and the values being a list [x, y, z] which give the numbers of [point-sharing, edge-sharing, and
    face-sharing] instances between the first Polyhedra and all of the reflections of the second Polyhedra
    """

    polyhedra_list = get_supercell_polyhedra(structure, radius, peripheral_species, central_species)

    # get the polyhedra in the central cell of the supercell from which to calculate connectivities
    center_cell_polyhedra = []
    cation_names = []
    for polyhedra in polyhedra_list:
        site_coords = polyhedra.central_ion.frac_coords
        if site_coords[0] >= (1.0/2)-(1e-5) and site_coords[0] < (2.0/3)-1e-5:
            if site_coords[1] >= (1.0/2)-(1e-5) and site_coords[1] < (2.0/3)-1e-5:
                if site_coords[2] >= (1.0/2)-(1e-5) and site_coords[2] < (2.0/3)-1e-5:
                    center_cell_polyhedra.append(polyhedra)
                    cation_names.append(polyhedra.central_ion_name)

    # instantiate nested dictionaries for connections
    connections = {}
    for cation_1 in cation_names:
        connections[cation_1] = {}
        for cation_2 in cation_names:
            connections[cation_1][cation_2] = {"point": 0, "edge": 0, "face": 0}

    # count connections between specific polyhedra to create adjacency matrix
    # matrix should be symmetric, i.e. connections between cation1 and cation2 should be the same as connections
    # between cation2 and cation1
    for inner_polyhedra in center_cell_polyhedra:
        inner_cation = inner_polyhedra.central_ion_name
        #print inner_cation
        for outer_polyhedra in polyhedra_list:
            outer_cation = outer_polyhedra.central_ion_name
            if inner_polyhedra == outer_polyhedra:
                continue
            connection = inner_polyhedra.get_num_connections(outer_polyhedra)
            if connection == 1:
                connections[inner_cation][outer_cation]["point"] += 1
            if connection == 2:
                connections[inner_cation][outer_cation]["edge"] += 1
            if connection >= 3:
                connections[inner_cation][outer_cation]["face"] += 1

    return connections

# TODO: using list as peripheral_species is a common Python "Gotcha!" Please look it up and fix.
# TODO: same thing with central_species

def get_connectivity_matrix_2(structure, radius = 2.8, peripheral_species=['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S'], central_species=[]):

    """
    Creates a connectivity matrix to describe the connectivity between cations in a structure; connections between
    cation polyhedra and reflections of itself (and reflections of other cation polyhedra) are also counted; different
    sites of the same species ARE differentiated (that is, connectivity for cations with different sites, not including
    reflections, are counted separately)

    :param structure: (Structure) target structure
    :param radius: (float) radius within which to determine whether a nearby atom is a peripheral ion
    :param peripheral_species: (list of Strings) List of strings with species names of the ions that can be peripheral
    ions
    :param central_species: (list of Strings) List of strings with species names of the ions that can be central ions;
    if there are no central species specified, then we consider all ions not in the list of peripheral ions to be
    central species
    :return: (dict) dictionary of dictionaries, with the first set of keys being the specified (unreflected) Polyhedra
    sites, the second (inner) set of keys being the specified (unreflected) Polyhedra sites which we are comparing
    connectivity to, and the values being a list [x, y, z] which give the numbers of [point-sharing, edge-sharing, and
    face-sharing] instances between the first Polyhedra and all of the reflections of the second Polyhedra
    """

    polyhedra_list = get_supercell_polyhedra(structure, radius, peripheral_species, central_species)

    # for all polyhedra in list
        # for all polyhedra in list
            # find connection between the polyhedra
            # check both sites (keep track of polyhedra of same species using site)
                # if either site is an image, record the connection is the original site's entry

    # get the polyhedra in the central cell of the supercell from which to calculate connectivities, relabeling
    # to differentiate between different sites with same cation species
    center_cell_polyhedra = []
    cation_iterators = {}
    for polyhedra in polyhedra_list:
        site_coords = polyhedra.central_ion.frac_coords
        if site_coords[0] >= (1.0/2)-1e-5 and site_coords[0] < (2.0/3)-1e-5:
            if site_coords[1] >= (1.0/2)-1e-5 and site_coords[1] < (2.0/3)-1e-5:
                if site_coords[2] >= (1.0/2)-1e-5 and site_coords[2] < (2.0/3)-1e-5  :
                    if polyhedra.central_ion_name not in cation_iterators.keys():
                        cation_iterators[polyhedra.central_ion_name] = 1
                    else:
                        cation_iterators[polyhedra.central_ion_name] += 1
                    polyhedra.set_site_number(cation_iterators[polyhedra.central_ion_name])
                    center_cell_polyhedra.append(polyhedra)


    # classify supercell polyhedra into center-cell polyhedra by checking whether they are images
    for outer_polyhedra in polyhedra_list:
        for inner_polyhedra in center_cell_polyhedra:
            if check_image_in_supercell(outer_polyhedra.central_ion, inner_polyhedra.central_ion):
                outer_polyhedra.set_site_number(inner_polyhedra.cation_num)
                break

    # instantiate nested dictionaries for connections
    connections = {}
    for cation_1 in center_cell_polyhedra:
        cation_1_unique_site_name = cation_1.central_ion_name + str(cation_1.cation_num)
        connections[cation_1_unique_site_name] = {}
        for cation_2 in center_cell_polyhedra:
            cation_2_unique_site_name = cation_2.central_ion_name + str(cation_2.cation_num)
            connections[cation_1_unique_site_name][cation_2_unique_site_name] = {"point": 0, "edge": 0, "face": 0}

    # count connections between specific polyhedra to create adjacency matrix
    # matrix should be symmetric, i.e. connections between cation1 and cation2 should be the same as connections
    # between cation2 and cation1
    for inner_polyhedra in center_cell_polyhedra:
        inner_cation = inner_polyhedra.central_ion_name + str(inner_polyhedra.cation_num)
        for outer_polyhedra in polyhedra_list:
            outer_cation = outer_polyhedra.central_ion_name + str(outer_polyhedra.cation_num)
            if inner_polyhedra == outer_polyhedra:
                continue
            connection = inner_polyhedra.get_num_connections(outer_polyhedra)
            if connection == 1:
                connections[inner_cation][outer_cation]["point"] += 1
            if connection == 2:
                connections[inner_cation][outer_cation]["edge"] += 1
            if connection >= 3:
                connections[inner_cation][outer_cation]["face"] += 1

    return connections


def check_image_in_supercell(site1, site2):
    is_image = False
    x1 = site1.frac_coords[0]
    x2 = site2.frac_coords[0]
    y1 = site1.frac_coords[1]
    y2 = site2.frac_coords[1]
    z1 = site1.frac_coords[2]
    z2 = site2.frac_coords[2]
    if round((x1 - x2) * 3, 5).is_integer() and round((y1 - y2) * 3, 5).is_integer() and round((z1 - z2) * 3, 5).is_integer():
        is_image = True

    return is_image


if __name__ == '__main__':

    # For some reason, can only get one connectivity matrix at a time or else the next connectivity matrix takes
    # a very long time to obtain (don't know if it's for it to finish)
    # TODO: regarding the comment above, see my note on default mutable args for peripheral_species and central_species

    # TODO: make these into unit tests
    print "Testing on BCC Fe"
    print "Note: for this situation, where the central ion's peripheral ions is the same species as the ion itself, " \
          "we need to specify both the central species and the peripheral species\n"
    s1 = Structure.from_file('Fe.cif', True, False)
    print s1
    central_species = ['Fe']
    peripheral_species = ['Fe']
    print get_connectivity_matrix(s1, 2.8, peripheral_species, central_species)
    print ""

    """
    print "Testing on CaF2"
    s2 = Structure.from_file('CaF2.cif', True, False)
    print s2
    print ""
    print "Identifying connectivity between Ca-centered polyhedra"
    central_species = ['Ca']
    peripheral_species = ['F']
    print get_connectivity_matrix(s2, 2.8, peripheral_species, central_species)
    print ""

    central_species = ['F']
    peripheral_species = ['Ca']
    print "Identifying connectivity between F-centered polyhedra"
    print get_connectivity_matrix(s2, 2.8, peripheral_species, central_species)
    print ""
    print "Identifying connectivity between F-centered polyhedra, distinguishing between the two different sites in the" \
          " structure"
    print get_connectivity_matrix_2(s2, 2.8, peripheral_species, central_species)
    print ""


    print "Testing on LiCoO2"
    print "Note: default value for peripheral species is " \
          "['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S']"
    print "Default value for central species is all species in the material that are not a peripheral species"
    s3 = Structure.from_file('LiCoO2.cif', True, False)
    print s3
    print get_connectivity_matrix(s3, 2.8)
    """

