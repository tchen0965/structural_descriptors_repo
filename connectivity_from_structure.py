__author__ = 'Tina_Chen'
__contributor__ = 'Anubhav Jain'


from pymatgen import Structure
from pymatgen import PeriodicSite
import effective_coordination as ECoN
from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
from itertools import combinations
import math
import numpy
import re


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
        """

        :param cation_number: number differentiating sites with the same species sitting on the site
        """
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
        :return: true if the central ion site of the current polyhedra is the same as the other polyhedra
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
        if ECoN.calculate_bond_weight(potentialPeripheralIons[1], bondlengths) > 0.5:
            peripheralSites.append(potentialPeripheralIons)

    return Polyhedra(site, peripheralSites)


def get_supercell_polyhedra(structure, radius, peripheral_species, central_species = [], n=3):

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
    supercell = structure.copy()

    supercell.make_supercell((n, n, n))

    polyhedra = []
    if len(central_species) == 0:
        for site_index in range(supercell.num_sites):
            if supercell[site_index].species_string not in peripheral_species:
                polyhedra.append(get_polyhedra(supercell, supercell[site_index], peripheral_species, radius))
    else:
        for site_index in range(supercell.num_sites):
            if supercell[site_index].species_string in central_species:
                polyhedra.append(get_polyhedra(supercell, supercell[site_index], peripheral_species, radius))

    return polyhedra

def get_supercell_size(structure, radius, peripheral_species, central_species):
    """
    Gets the side length of supercell to build in order to obtain all necessary polyhedra; gives a 6x6x6 supercell only
    if the polyhedra extend across multiple unit cells, or else uses a 3x3x3 supercell. We determine what supercell size
    to use by looking at an example polyhedra. If the test polyhedra extends beyond an entire unit cell, then it needs a
    larger super cell. This specifically addresses cases such as BCC Fe in which the polyhedra spans multiple unit cells
    and thus needs a larger supercell to consider connectivities to all surrounding polyhedra

    :param structure: (Structure) target structure
    :param radius: (float) radius within which to determine whether a nearby atom is a peripheral ion
    :param peripheral_species: (list of Strings) List of strings with species names of the ions that can be peripheral
    ions
    :param central_species: (list of Strings) List of strings with species names of the ions that can be central ions;
    if there are no central species specified, then we consider all ions not in the list of peripheral ions to be
    central species
    :return: (integer) side length of cell to build (in unit cells)
    """

    # obtain a test polyhedra based on the given structure and given peripheral and central species
    if len(central_species) == 0:
        for site_index in range(structure.num_sites):
            if structure[site_index].species_string not in peripheral_species:
                test_polyhedra = get_polyhedra(structure, structure[site_index], peripheral_species, radius)
                break
    else:
        for site_index in range(structure.num_sites):
            if structure[site_index].species_string in central_species:
                test_polyhedra = get_polyhedra(structure, structure[site_index], peripheral_species, radius)
                break

    # find the bond lengths between the central ion site and each of the peripheral ion sites
    frac_bond_lengths = []
    for peripheral_ion in test_polyhedra.peripheral_ions:
        frac_bond_lengths.append(numpy.linalg.norm(test_polyhedra.central_ion.frac_coords - peripheral_ion.frac_coords))

    # if the maximum bond length is greater than half a unit cell, use a larger supercell size; or else, use a
    # normal 3x3x3 supercell which should give necessary polyhedra surrounding the central cell
    if (max(frac_bond_lengths) >= 0.5-1e-5):
        supercell_size = 6
    else:
        supercell_size = 3

    return supercell_size

def get_central_cell_polyhedra(polyhedra_list, supercell_size, sites_diff):
    """
    Obtains all Polyhedra in the center unit cell of the supercell

    :param polyhedra_list: (list of Polyhedra) list of all Polyhedra in the super cell
    :param supercell_size: (integer) side length of super cell
    :return: (list of Polyhedra) list of Polyhedra that are in the center cell of the super cell
    :return: (list of String) list of all cations in the structure
    """
    center_cell_polyhedra = []
    cation_iterators = {}
    lower_boundary = (supercell_size/2)/float(supercell_size)
    upper_boundary = ((supercell_size/2)+1)/float(supercell_size)
    for polyhedra in polyhedra_list:
        site_coords = polyhedra.central_ion.frac_coords
        if site_coords[0] >= lower_boundary-(1e-5) and site_coords[0] < upper_boundary-(1e-5):
            if site_coords[1] >= lower_boundary-(1e-5) and site_coords[1] < upper_boundary-(1e-5):
                if site_coords[2] >= lower_boundary-(1e-5) and site_coords[2] < upper_boundary-(1e-5):
                    if polyhedra.central_ion_name not in cation_iterators.keys():
                        cation_iterators[polyhedra.central_ion_name] = 1
                    else:
                        cation_iterators[polyhedra.central_ion_name] += 1
                    if sites_diff:
                        polyhedra.set_site_number(cation_iterators[polyhedra.central_ion_name])
                    center_cell_polyhedra.append(polyhedra)

    return center_cell_polyhedra, cation_iterators.keys()

"""
def get_central_cell_polyhedra_sites_diff(polyhedra_list, supercell_size):

    Obtains all Polyhedra in the center unit cell of the supercell, keeping track of different sites with the same
    species occupying the site

    :param polyhedra_list: (list of Polyhedra) list of all Polyhedra in the super cell
    :param supercell_size: (integer) side length of super cell
    :return: (list) list of Polyhedra that are in the center cell of the super cell, with the cation site number set to
    identify different sites with the same occupying species
    :return: (dictionary) dict of keys = cation name, values = number of unique sites of that cation found in the
    central unit cell

    center_cell_polyhedra = []
    cation_iterators = {}
    lower_boundary = (supercell_size/2)/float(supercell_size)
    upper_boundary = ((supercell_size/2)+1)/float(supercell_size)
    for polyhedra in polyhedra_list:
        site_coords = polyhedra.central_ion.frac_coords
        if site_coords[0] >= lower_boundary-(1e-5) and site_coords[0] < upper_boundary-(1e-5):
            if site_coords[1] >= lower_boundary-(1e-5) and site_coords[1] < upper_boundary-(1e-5):
                if site_coords[2] >= lower_boundary-(1e-5) and site_coords[2] < upper_boundary-(1e-5):
                    if polyhedra.central_ion_name not in cation_iterators.keys():
                        cation_iterators[polyhedra.central_ion_name] = 1
                    else:
                        cation_iterators[polyhedra.central_ion_name] += 1
                    polyhedra.set_site_number(cation_iterators[polyhedra.central_ion_name])
                    center_cell_polyhedra.append(polyhedra)
    return center_cell_polyhedra, cation_iterators
"""

def get_ex_poly(polyhedra_list, species_tag):
    """

    :param polyhedra_list: (list of Polyhedra) list of all Polyhedra in the super cell structure
    :param species_tag: (String) name of species given by a connectivity matrix; may contain the cation site number
    appended to the species name
    :return: (Polyhedra) first Polyhedra found in polyhedra list that matches the species_tag
    """
    species_name = re.sub(r"[^A-Za-z]+", '', species_tag)
    if (len(re.findall(r'\d+', species_tag))==0):
        print species_tag
        print re.findall(r'\d+', species_tag)
    else:
        species_site_num = int(re.findall(r'\d+', species_tag)[0])
    for polyhedra in polyhedra_list:
        if polyhedra.central_ion_name == species_name:
            if polyhedra.cation_num == species_site_num:
                return polyhedra


def get_CN_description(CN):
    """
    Writes a verbal description of the coordination number/polyhedra based on a given coordination number
    :param CN: (integer) rounded coordination number
    :return: (String) verbal description of coordination number/polyhedra
    """
    roundedCN = round(CN)
    description = str(int(roundedCN)) + "-fold coordinated"
    return description

def check_image_in_supercell(site1, site2, supercell_size):
    """
    Checks whether site1 and site2 are periodic images of each other in the super cell structure given the size of the
    super cell
    :param site1: (Site) site in super cell
    :param site2: (Site) site in super cell
    :param supercell_size: (integer) side length of super cell (in unit cells)
    :return: (boolean) whether site1 and site2 are periodic images of each other in the super cell
    """
    is_image = False
    x1 = site1.frac_coords[0]
    x2 = site2.frac_coords[0]
    y1 = site1.frac_coords[1]
    y2 = site2.frac_coords[1]
    z1 = site1.frac_coords[2]
    z2 = site2.frac_coords[2]
    if round((x1 - x2) * supercell_size, 5).is_integer() and \
            round((y1 - y2) * supercell_size, 5).is_integer() and \
            round((z1 - z2) * supercell_size, 5).is_integer():
        is_image = True

    return is_image


def get_connectivity_matrix(structure, sites_diff, radius = 2.8, peripheral_species=['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S'], central_species = []):
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
    :return: (list of Polyhedra) list of all Polyhedra in the super cell
    """

    # identify supercell size (3x3x3 vs. 6x6x6) needed to calculate connectivities based on size of polyhedra
    supercell_size = get_supercell_size(structure, radius, peripheral_species, central_species)

    # obtain list of all polyhedra in the supercell
    polyhedra_list = get_supercell_polyhedra(structure, radius, peripheral_species, central_species, supercell_size)

    # get the polyhedra in the central cell of the supercell from which to calculate connectivities
    center_cell_polyhedra, cation_names = get_central_cell_polyhedra(polyhedra_list, supercell_size, sites_diff)

    # for matrix with differentiated cation sites, classify outer polyhedra sites
    if sites_diff:
        for outer_polyhedra in polyhedra_list:
            for inner_polyhedra in center_cell_polyhedra:
                if check_image_in_supercell(outer_polyhedra.central_ion, inner_polyhedra.central_ion, supercell_size):
                    outer_polyhedra.set_site_number(inner_polyhedra.cation_num)
                    break

    # instantiate nested dictionaries for connections
    connections = {}
    if sites_diff:
        for cation_1 in center_cell_polyhedra:
            cation_1_unique_site_name = cation_1.central_ion_name + str(cation_1.cation_num)
            connections[cation_1_unique_site_name] = {}
            for cation_2 in center_cell_polyhedra:
                cation_2_unique_site_name = cation_2.central_ion_name + str(cation_2.cation_num)
                connections[cation_1_unique_site_name][cation_2_unique_site_name] = {"point": 0, "edge": 0, "face": 0}
    else:
        for cation_1 in cation_names:
            connections[cation_1] = {}
            for cation_2 in cation_names:
                connections[cation_1][cation_2] = {"point": 0, "edge": 0, "face": 0}

    # count connections between specific polyhedra to create adjacency matrix
    # matrix should be symmetric, i.e. connections between cation1 and cation2 should be the same as connections
    # between cation2 and cation1
    for inner_polyhedra in center_cell_polyhedra:
        if sites_diff:
            inner_cation = inner_polyhedra.central_ion_name + str(inner_polyhedra.cation_num)
        else:
            inner_cation = inner_polyhedra.central_ion_name
        for outer_polyhedra in polyhedra_list:
            if sites_diff:
                outer_cation = outer_polyhedra.central_ion_name + str(outer_polyhedra.cation_num)
            else:
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

    return connections, polyhedra_list

"""
def get_connectivity_matrix_sites_diff(structure, sites_diff, radius = 2.8, peripheral_species=['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S'], central_species=[]):

<<<<<<< HEAD
=======
# TODO: using list as peripheral_species is a common Python "Gotcha!" Please look it up and fix.
# TODO: same thing with central_species

def get_connectivity_matrix_2(structure, radius = 2.8, peripheral_species=['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S'], central_species=[]):
>>>>>>> origin/master

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
    :return: (List) list of all Polyhedra in the super cell, with cation site numbers assigned



    supercell_size = get_supercell_size(structure, radius, peripheral_species, central_species)

    polyhedra_list = get_supercell_polyhedra(structure, radius, peripheral_species, central_species, supercell_size)

    # get the polyhedra in the central cell of the supercell from which to calculate connectivities, relabeling
    # to differentiate between different sites with same cation species
    center_cell_polyhedra, iterators = get_central_cell_polyhedra(polyhedra_list, supercell_size, sites_diff)


    # classify supercell polyhedra into center-cell polyhedra by checking whether they are images
    for outer_polyhedra in polyhedra_list:
        for inner_polyhedra in center_cell_polyhedra:
            if check_image_in_supercell(outer_polyhedra.central_ion, inner_polyhedra.central_ion, supercell_size):
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

    return connections, polyhedra_list
"""

def get_connectivity_description(connectivity_matrix, polyhedra, structure, sites_diff, radius=2.8, anions=['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S']):
    """
    Writes a verbal description of the connectivity between cations in a structure; connections between cation
    polyhedra and reflections of itself (and reflections of other cation polyhedra) are also counted; different
    sites of the same species are not differentiated (that is, connectivity for cations with different sites, not
    including reflections, are not counted separately)

    :param connectivity_matrix: (dict) dicitonary of dictionaries containing connectivities between different cation
    polyhedra (should be output from get_connectivity_matrix)
    :param structure: (Structure) target structure
    :param radius: (float) radius within which to determine whether a nearby atom is a peripheral ion
    :return: (dict) dictionary of strings containing verbal descriptions of connectivites between cations in the given
    structure
    """

    CN_list = ECoN.EffectiveCoordFinder(structure).get_avg_CN(radius, anions)
    descriptions = {}
    for cation_1 in connectivity_matrix.keys():
        descriptions[cation_1] = ""
    for cation_1 in connectivity_matrix.keys():
        if sites_diff:
            CN1 = ECoN.get_effective_CN(get_ex_poly(polyhedra, cation_1))
        else:
            CN1 = CN_list[cation_1]
        descriptions[cation_1] += cation_1 + " are " + get_CN_description(CN1) + ". \n"
        for cation_2 in connectivity_matrix[cation_1].keys():
            connected = False
            for connectivity_type in connectivity_matrix[cation_1][cation_2].keys():
                if connectivity_matrix[cation_1][cation_2][connectivity_type] != 0:
                    connected = True
            if connected == True:
                descriptions[cation_1] += "They are "
                first = True
                for connectivity_type in connectivity_matrix[cation_1][cation_2].keys():
                    if connectivity_matrix[cation_1][cation_2][connectivity_type] != 0:
                        if first:
                            descriptions[cation_1] += connectivity_type + "-connected "
                            first = False
                        else:
                            descriptions[cation_1] += "and " + connectivity_type + "-connected "
                if sites_diff:
                    CN2 = ECoN.get_effective_CN(get_ex_poly(polyhedra, cation_2))
                else:
                    CN2 = CN_list[cation_2]
                descriptions[cation_1] += "to " + get_CN_description(CN2) + " " + cation_2 + ". "
    return descriptions

"""
def get_connectivity_description_sites_diff(connectivity_matrix, polyhedra):

    Writes a verbal description of the connectivity between cations in a structure; connections between cation
    polyhedra and reflections of itself (and reflections of other cation polyhedra) are also counted; different
    sites of the same species are differentiated (that is, connectivity for cations with different sites, not
    including reflections, are counted separately)

    :param connectivity_matrix: (dict) dicitonary of dictionaries containing connectivities between different cation
    polyhedra (should be output from get_connectivity_matrix)
    :param polyhedra: (list of Polyhedra) list of all Polyhedra in the supercell
    :return: (dict) dictionary of strings containing verbal descriptions of connectivites between cations in the given
    structure, differentiating between different sites with the same occupying species

    descriptions = {}
    for cation_1 in connectivity_matrix.keys():
        descriptions[cation_1] = ""
    for cation_1 in connectivity_matrix.keys():
        descriptions[cation_1] += cation_1 + " are " + \
                                  get_CN_description(ECoN.get_effective_CN(get_ex_poly(polyhedra, cation_1))) + ". \n"
        for cation_2 in connectivity_matrix[cation_1].keys():
            connected = False
            for connectivity_type in connectivity_matrix[cation_1][cation_2].keys():
                if connectivity_matrix[cation_1][cation_2][connectivity_type] != 0:
                    connected = True
            if connected:
                descriptions[cation_1] += "They are "
                first = True
                for connectivity_type in connectivity_matrix[cation_1][cation_2].keys():
                    if connectivity_matrix[cation_1][cation_2][connectivity_type] != 0:
                        if first:
                            descriptions[cation_1] += connectivity_type + "-connected "
                            first = False
                        else:
                            descriptions[cation_1] += "and " + connectivity_type + "-connected "
                descriptions[cation_1] += "to " + \
                                          get_CN_description(ECoN.get_effective_CN(get_ex_poly(polyhedra, cation_2))) \
                                          + " " + cation_2 + ". "
    return descriptions
"""


def get_surrounding_connectivity(structure, polyhedra, radius=2.8, peripheral_species=['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S'], central_species = []):
    """
    Gives the surrounding connectivity of a specific polyhedra in a given structure

    :param structure: (Structure) target structure
    :param polyhedra: (Polyhedra) target polyhedra
    :return: (list) list of polyhedra and their connections to the target polyhedra in the target structure, elements
    of the list are of the form: [[polyhedra, connection], ...]
    """
    targetSite = polyhedra.central_ion

    polyhedraList = get_supercell_polyhedra(structure, radius, peripheral_species, central_species, 6)

    targetPolyhedra = polyhedra

    # identify equivalent polyhedra in one of center unit cells
    for poly in polyhedraList:
        site = poly.central_ion
        if (targetSite.species_string == site.species_string):
            if (abs(site.frac_coords[0] - 0.5-targetSite.frac_coords[0]) < 0.01):
                if (abs(site.frac_coords[1] - 0.5-targetSite.frac_coords[1]) < 0.01):
                    if (abs(site.frac_coords[2] - 0.5-targetSite.frac_coords[2]) < 0.01):
                        targetPolyhedra = poly

    connectedPolyhedra = []
    for otherPolyhedra in polyhedraList:
        if otherPolyhedra.central_ion != targetPolyhedra.central_ion:
            connection = targetPolyhedra.get_num_connections(otherPolyhedra)
            if connection > 0:
                connectedPolyhedra.append([otherPolyhedra, connection])

    return connectedPolyhedra



if __name__ == '__main__':

<<<<<<< HEAD

=======
    # For some reason, can only get one connectivity matrix at a time or else the next connectivity matrix takes
    # a very long time to obtain (don't know if it's for it to finish)
    # TODO: regarding the comment above, see my note on default mutable args for peripheral_species and central_species

    # TODO: make these into unit tests
>>>>>>> origin/master
    print "Testing on BCC Fe"
    print "Note: for this situation, where the central ion's peripheral ions is the same species as the ion itself, " \
          "we need to specify both the central species and the peripheral species\n"
    s1 = Structure.from_file('Fe.cif', True, False)
    print s1
    central_species = ['Fe']
    peripheral_species = ['Fe']
    Fematrix, Fepolyhedra = get_connectivity_matrix(s1, True, 2.8, peripheral_species, central_species)
    print ""
    print Fematrix
    print ""
    Fedescriptions = get_connectivity_description(Fematrix, Fepolyhedra, s1, True)
    for cation_1 in Fedescriptions.keys():
        print Fedescriptions[cation_1]
        print ""

    print "Testing on CaF2"

    s2 = Structure.from_file('CaF2.cif', True, False)

    print s2
    print ""
    print "Identifying connectivity between Ca-centered polyhedra"
    central_species = ['Ca']
    peripheral_species = ['F']
    CaFmatrix, CaFpolyhedra = get_connectivity_matrix(s2, True, 2.8, peripheral_species, central_species)
    print ""
    print CaFmatrix
    print ""
    CaFdescriptions = get_connectivity_description(CaFmatrix, CaFpolyhedra, s2, True)
    for cation_1 in CaFdescriptions.keys():
        print CaFdescriptions[cation_1]
        print ""


    s2 = Structure.from_file('CaF2.cif', True, False)
    central_species = ['F']
    peripheral_species = ['Ca']
    #print "Identifying connectivity between F-centered polyhedra"
    #print get_connectivity_matrix(s2, 2.8, peripheral_species, central_species)
    #print ""
    print "Identifying connectivity between F-centered polyhedra, distinguishing between the two different sites in the" \
          " structure"
    print ""
    FCamatrix, FCapolyhedra = get_connectivity_matrix(s2, True, 2.8, peripheral_species, central_species)
    print FCamatrix
    print ""
    FCadescriptions = get_connectivity_description(FCamatrix, FCapolyhedra, s2, True)
    for cation_1 in FCadescriptions.keys():
        print FCadescriptions[cation_1]
        print ""


    print "Testing on LiCoO2"
    print "Note: default value for peripheral species is " \
          "['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S']"
    print "Default value for central species is all species in the material that are not a peripheral species"
    s3 = Structure.from_file('LiCoO2.cif', True, False)
    print s3
    print ""
    LiCoO2matrix, LiCoO2polyhedra = get_connectivity_matrix(s3, True, 2.8)
    print LiCoO2matrix
    print ""
    LiCoO2descriptions = get_connectivity_description(LiCoO2matrix, LiCoO2polyhedra, s3, True)
    for cation_1 in LiCoO2descriptions.keys():
        print LiCoO2descriptions[cation_1]
        print ""



