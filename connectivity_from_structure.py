from pymatgen.analysis.chemenv.coordination_environments import coordination_geometry_finder as polyfinder
from pymatgen.analysis.chemenv.coordination_environments import chemenv_strategies as strategies
from pymatgen.util.coord_utils import find_in_coord_list_pbc
from pymatgen.analysis.chemenv.coordination_environments import structure_environments as se
from effective_coordination import calculate_bond_weight, EffectiveCoordFinder, get_effective_cn
from pymatgen.analysis.chemenv.coordination_environments import coordination_geometries_files as cg_files
from polyhedra import Polyhedra
from chemenv_util import find_species_ce_from_light_se
import json
import os
import numpy
import re


__author__ = 'Tina_Chen'
__contributor__ = 'Anubhav Jain'

# TODO: please talk to me about overall code organization
# TODO: please talk to me about adding a README file that shows you how to use the library

"""
Finds the connectivities between polyhedra within a structure, where we define connectivity by the number of shared
peripheral ions between two polyhedra. 0 shared peripheral ions indicates that there is no connectivity between the
two polyhedra; 1 shared peripheral ion indicates that there is point-sharing connectivity between the polyhedra; 2
shared peripheral ions indicate that there is edge-sharing connectivity between the polyhedra; 3 or more shared
peripheral ions indicate that there is face-sharing connectivity between the polyhedra.

For the purpose of this code, the peripheral ions will be defined by the ions surrounding a central site as given by
structure's get_neighbors function, with bond weight contributions (given by Hoppe, 1979) less than a set cut-off value.

"""


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

    anion_sites = []
    bondlengths = []
    peripheral_sites = []

    for entry in structure.get_neighbors(site, radius):
        if entry[1] < radius:
            if entry[0].species_string in peripheral_species:
                anion_sites.append(entry[0])
                bondlengths.append(entry[1])
    for potentialPeripheralIons in anion_sites:
        # do not count nearby ions that do not contribute
        if calculate_bond_weight(potentialPeripheralIons[1], bondlengths) > 0:
            peripheral_sites.append(potentialPeripheralIons)

    return Polyhedra(site, peripheral_sites)


def get_polyhedra_in_structure(structure, radius, peripheral_species, central_species):

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

    # obtain structure's polyhedra based on the given structure and given peripheral and central species
    test_polyhedra = []

    if len(central_species) == 0:
        for site_index in range(structure.num_sites):
            if structure[site_index].species_string not in peripheral_species:
                test_polyhedra.append(get_polyhedra(structure, structure[site_index], peripheral_species, radius))
    else:
        for site_index in range(structure.num_sites):
            if structure[site_index].species_string in central_species:
                test_polyhedra.append(get_polyhedra(structure, structure[site_index], peripheral_species, radius))

    # find the bond lengths between the central ion site and each of the peripheral ion sites
    frac_bond_lengths = []
    for polyhedra in test_polyhedra:
        for peripheral_ion in polyhedra.peripheral_ions:
            frac_bond_lengths.append(numpy.linalg.norm(polyhedra.central_ion.frac_coords - peripheral_ion.frac_coords))

    # if the maximum bond length is greater than half a unit cell, use a larger supercell size; or else, use a
    # normal 3x3x3 supercell which should give necessary polyhedra surrounding the central cell
    if max(frac_bond_lengths) >= 0.5-1e-5:
        supercell_size = 6
    else:
        supercell_size = 3

    return supercell_size


def get_central_cell_polyhedra(polyhedra_list, supercell_size, sites_diff):
    """
    Obtains all Polyhedra in the center unit cell of the supercell

    :param polyhedra_list: (list of Polyhedra) list of all Polyhedra in the super cell
    :param supercell_size: (integer) side length of super cell
    :param sites_diff: (Boolean) whether sites with same cation species should be differentiated
    :return: (list of Polyhedra) list of Polyhedra that are in the center cell of the super cell
    :return: (list of String) list of all cations in the structure
    """
    center_cell_polyhedra = []
    cation_iterators = {}
    lower_boundary = (supercell_size/2)/float(supercell_size)
    upper_boundary = ((supercell_size/2)+1)/float(supercell_size)
    for polyhedra in polyhedra_list:
        site_coords = polyhedra.central_ion.frac_coords
        if upper_boundary-1e-5 > site_coords[0] >= lower_boundary-1e-5:
            if upper_boundary-1e-5 > site_coords[1] >= lower_boundary-1e-5:
                if upper_boundary-1e-5 > site_coords[2] >= lower_boundary-1e-5:
                    if polyhedra.central_ion_name not in cation_iterators.keys():
                        cation_iterators[polyhedra.central_ion_name] = 1
                    else:
                        cation_iterators[polyhedra.central_ion_name] += 1
                    if sites_diff:
                        polyhedra.set_site_number(cation_iterators[polyhedra.central_ion_name])
                    center_cell_polyhedra.append(polyhedra)

    return center_cell_polyhedra, cation_iterators.keys()


def get_ex_poly(polyhedra_list, species_tag):
    """

    :param polyhedra_list: (list of Polyhedra) list of all Polyhedra in the super cell structure
    :param species_tag: (String) name of species given by a connectivity matrix; may contain the cation site number
    appended to the species name
    :return: (Polyhedra) first Polyhedra found in polyhedra list that matches the species_tag
    """
    species_name = re.sub(r"[^A-Za-z]+", '', species_tag)
    if len(re.findall(r'\d+', species_tag)) == 0:
        raise ValueError('Site-differentiated species tag has no site number')
    else:
        species_site_num = int(re.findall(r'\d+', species_tag)[0])
    for polyhedra in polyhedra_list:
        if polyhedra.central_ion_name == species_name:
            if polyhedra.cation_num == species_site_num:
                return polyhedra


def get_cn_description(cn):
    """
    Writes a verbal description of the coordination number/polyhedra based on a given coordination number
    :param cn: (integer) rounded coordination number
    :return: (String) verbal description of coordination number/polyhedra
    """
    rounded_cn = round(cn)
    description = str(int(rounded_cn)) + "-fold coordinated"
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


def get_connectivity_matrix(structure, sites_diff=True, radius=2.8, peripheral_species=None, central_species=None):
    """
    Creates a connectivity matrix to describe the connectivity between cations in a structure; connections between
    cation polyhedra and reflections of itself (and reflections of other cation polyhedra) are also counted; different
    sites of the same cation species are NOT differentiated (that is, all instances of each cation are counted all
    together)

    :param structure: (Structure) target structure
    :param sites_diff: (Boolean) whether sites with the same cation species should be differentiated
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

    if peripheral_species is None:
        peripheral_species = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S', 'N', 'N3-']
    if central_species is None:
        central_species = []

    # identify supercell size (3x3x3 vs. 6x6x6) needed to calculate connectivities based on size of polyhedra
    supercell_size = get_supercell_size(structure, radius, peripheral_species, central_species)
    supercell = structure.copy()
    supercell.make_supercell([supercell_size, supercell_size, supercell_size])

    # obtain list of all polyhedra in the supercell
    polyhedra_list = get_polyhedra_in_structure(supercell, radius, peripheral_species, central_species)

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

    if not sites_diff:
        comp = structure.composition
        for inner_cation in connections.keys():
            for outer_cation in connections[inner_cation].keys():
                connections[inner_cation][outer_cation]['face'] = \
                    connections[inner_cation][outer_cation]['face']/comp[inner_cation]
                connections[inner_cation][outer_cation]['edge'] = \
                    connections[inner_cation][outer_cation]['edge']/comp[inner_cation]
                connections[inner_cation][outer_cation]['point'] = \
                    connections[inner_cation][outer_cation]['point']/comp[inner_cation]

    return connections, polyhedra_list, supercell


# TODO: please discuss the method below with me. It should not require setting a radius or anions.
# TODO: The user should be able to provide an appropriate object and this should just take care of describing that
# TODO: object. It should not start analyzing the structure using its own choice of algorithms...
def get_connectivity_description(connectivity_matrix, polyhedra, structure, sites_diff, radius=2.6, anions=None):
    """
    Writes a verbal description of the connectivity between cations in a structure; connections between cation
    polyhedra and reflections of itself (and reflections of other cation polyhedra) are also counted; different
    sites of the same species are not differentiated (that is, connectivity for cations with different sites, not
    including reflections, are not counted separately)

    :param connectivity_matrix: (dict) dictionary of dictionaries containing connectivities between different cation
    polyhedra (should be output from get_connectivity_matrix)
    :param polyhedra: (List of Polyhedra) list of all Polyhedra in the supercell structure
    :param structure: (Structure) target structure
    :param sites_diff: (Boolean) whether sites with the same cation species should be differentiated
    :param radius: (float) radius within which to determine whether a nearby atom is a peripheral ion
    :param anions: (List of Strings) list of species strings of what we consider anions in the structure
    :return: (dict) dictionary of strings containing verbal descriptions of connectivites between cations in the given
    structure; keys = cation tag (eg: Li1, Na2, where the letters represent the species and the number distinguishes
    the species from the other instances of that same species in the unit cell); dict values = String descriptions of
    connectivity
    """

    if anions is None:
        anions = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S', 'N', 'N3-']

    cn_list = EffectiveCoordFinder(structure).get_avg_cn(radius, anions)
    descriptions = {}
    for cation_1 in connectivity_matrix.keys():
        descriptions[cation_1] = ""
    for cation_1 in connectivity_matrix.keys():
        if sites_diff:
            cn1 = get_effective_cn(get_ex_poly(polyhedra, cation_1))
        else:
            cn1 = cn_list[cation_1]
        descriptions[cation_1] += cation_1 + " are " + get_cn_description(cn1) + ". \n"
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
                if sites_diff:
                    cn2 = get_effective_cn(get_ex_poly(polyhedra, cation_2))
                else:
                    cn2 = cn_list[cation_2]
                descriptions[cation_1] += "to " + get_cn_description(cn2) + " " + cation_2 + ". "
    return descriptions


def get_connectivity_description_1(structure, radius=2.6, peripheral_species=None, central_species=None):
    """
    Writes a verbal description of the connectivity between cations in a structure and does differentiate between
    species on different sites in the unit cell. Uses the Chemenv module in Pymatgen to calculate the coordination
    environment for specific sites.

    :param structure: (Structure) target Structure
    :param radius: (float) radius within which to determine whether a nearby atom is a peripheral ion
    :param peripheral_species: (List of Strings) list of species strings of what we consider the peripheral species of
    the polyhedra in the structure
    :param central_species: (List of Strings) list of species strings of what we consider the central species of the
    polyhedra in the structure
    :return: (dict) dictionary of strings containing verbal descriptions of connectivites between cations in the given
    structure; keys = cation tag (eg: Li1, Na2, where the letters represent the species and the number distinguishes
    the species from the other instances of that same species in the unit cell); dict values = String descriptions of
    connectivity
    """

    if peripheral_species is None:
        peripheral_species = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S', 'N', 'N3-']
    if central_species is None:
        central_species = []

    connectivity_matrix, all_polyhedra, supercell = \
        get_connectivity_matrix(structure, True, radius, peripheral_species, central_species)

    lgf = polyfinder.LocalGeometryFinder()
    lgf.setup_parameters(centering_type='standard', structure_refinement='none')
    lgf.setup_structure(supercell)
    environments = lgf.compute_structure_environments_detailed_voronoi(maximum_distance_factor=1.5)
    light_se = se.LightStructureEnvironments(strategies.SimplestChemenvStrategy(), environments)
    path_to_jsons = os.path.dirname(cg_files.__file__)

    descriptions = {}
    for cation_1 in connectivity_matrix.keys():
        descriptions[cation_1] = ""
    for cation_1 in connectivity_matrix.keys():
        example_poly1 = get_ex_poly(all_polyhedra, cation_1)
        target_isite1 = find_in_coord_list_pbc(supercell.frac_coords, example_poly1.central_ion.frac_coords)[0]
        print supercell[target_isite1]
        print light_se._coordination_environments
        site_environment1 = light_se._coordination_environments[target_isite1]['ce_symbol']
        with open(path_to_jsons+"/%s.json" % site_environment1) as json_file:
            data1 = json.load(json_file)
        cn1 = data1['name']

        descriptions[cation_1] += cation_1 + " are " + cn1 + ". \n"
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

                example_poly2 = get_ex_poly(all_polyhedra, cation_1)
                target_isite2 = find_in_coord_list_pbc(supercell.frac_coords, example_poly2.central_ion.frac_coords)[0]
                print supercell[target_isite2]
                print light_se._coordination_environments
                site_environment2 = light_se._coordination_environments[target_isite2]['ce_symbol']
                with open(path_to_jsons+"/%s.json" % site_environment2) as json_file:
                    data2 = json.load(json_file)
                cn2 = data2['name']
                descriptions[cation_1] += "to " + cn2 + " " + cation_2 + ". "
    return descriptions


def get_connectivity_description_2(structure, radius=2.6, peripheral_species=None, central_species=None):
    """
    Writes a verbal description of the connectivity between cations in a structure and does not differentiate between
    species on different sites in the unit cell. Uses the Chemenv module in Pymatgen to calculate the coordination
    environment for specific species.

    :param structure: (Structure) target Structure
    :param radius: (float) radius within which to determine whether a nearby atom is a peripheral ion
    :param peripheral_species: (List of Strings) list of species strings of what we consider the peripheral species of
    the polyhedra in the structure
    :param central_species: (List of Strings) list of species strings of what we consider the central species of the
    polyhedra in the structure
    :return: (dict) dictionary of strings containing verbal descriptions of connectivites between cations in the given
    structure; keys = cation tag (eg: Li1, Na2, where the letters represent the species and the number distinguishes
    the species from the other instances of that same species in the unit cell); dict values = String descriptions of
    connectivity
    """

    if peripheral_species is None:
        peripheral_species = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S', 'N', 'N3-']
    if central_species is None:
        central_species = []

    connectivity_matrix, all_polyhedra, supercell = \
        get_connectivity_matrix(structure, True, radius, peripheral_species, central_species)

    lgf = polyfinder.LocalGeometryFinder()
    lgf.setup_parameters(centering_type='standard', structure_refinement='none')
    lgf.setup_structure(supercell)
    environments = lgf.compute_structure_environments_detailed_voronoi(maximum_distance_factor=1.5)
    light_se = se.LightStructureEnvironments(strategies.SimplestChemenvStrategy(), environments)
    path_to_jsons = os.path.dirname(cg_files.__file__)

    descriptions = {}
    for cation_1 in connectivity_matrix.keys():
        descriptions[cation_1] = ""
    for cation_1 in connectivity_matrix.keys():
        site_ces1 = find_species_ce_from_light_se(light_se, cation_1)
        cn1 = []
        for ce in site_ces1:
            with open(path_to_jsons+"/%s.json" % ce) as json_file:
                data1 = json.load(json_file)
            cn1.append(data1['name'])

        descriptions[cation_1] += cation_1 + " are " + cn1 + ". \n"
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

                site_ces2 = find_species_ce_from_light_se(light_se, cation_1)
                cn2 = []
                for ce in site_ces2:
                    with open(path_to_jsons+"/%s.json" % ce) as json_file:
                        data2 = json.load(json_file)
                    cn2.append(data2['name'])
                descriptions[cation_1] += "to " + cn2 + " " + cation_2 + ". "
    return descriptions


def get_surrounding_connectivity(structure, polyhedra, radius=2.6, peripheral_species=None, central_species=None):
    """
    Gives the surrounding connectivity of a specific polyhedra in a given structure

    :param structure: (Structure) target structure
    :param polyhedra: (Polyhedra) target polyhedra
    :param radius: (float) radius out to which we find surrounding neighbors
    :param peripheral_species: (list of Strings) List of strings with species names of the ions that can be peripheral
    ions
    :param central_species: (list of Strings) List of strings with species names of the ions that can be central ions;
    if there are no central species specified, then we consider all ions not in the list of peripheral ions to be
    central species
    :return: (list) list of polyhedra and their connections to the target polyhedra in the target structure, elements
    of the list are of the form: [[polyhedra, connection], ...]
    """
    if peripheral_species is None:
        peripheral_species = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S', 'N', 'N3-']
    if central_species is None:
        central_species = []

    target_site = polyhedra.central_ion

    supercell = structure.copy()
    n = get_supercell_size(structure, radius, peripheral_species, central_species)
    supercell.make_supercell([n, n, n])

    polyhedra_list = get_polyhedra_in_structure(supercell, radius, peripheral_species, central_species)

    target_polyhedra = polyhedra

    # identify equivalent polyhedra in one of center unit cells
    for poly in polyhedra_list:
        site = poly.central_ion
        if target_site.species_string == site.species_string:
            if abs(site.frac_coords[0] - 0.5-target_site.frac_coords[0]) < 0.01:
                if abs(site.frac_coords[1] - 0.5-target_site.frac_coords[1]) < 0.01:
                    if abs(site.frac_coords[2] - 0.5-target_site.frac_coords[2]) < 0.01:
                        target_polyhedra = poly

    connected_polyhedra = []
    for otherPolyhedra in polyhedra_list:
        if otherPolyhedra.central_ion != target_polyhedra.central_ion:
            connection = target_polyhedra.get_num_connections(otherPolyhedra)
            if connection > 0:
                connected_polyhedra.append([otherPolyhedra, connection])

    return connected_polyhedra
