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
        # TODO: document. What is the type of cation, e.g. is a Site object? Is peripheral ions a list of Sites?
        self.central_ion = cation
        self.central_ion_name = cation.species_string
        self.peripheral_ions = peripheralIons

        self.composition = self.central_ion._species
        for site in self.peripheral_ions:
            self.composition += site._species

    def get_central_name(self):  # TODO: why is a name needed? Why not just use self.central_ion.species_string?
        return self.central_ion_name

    def get_peripheral_sites(self): #TODO: a getter is NOT needed. The var peripheral_ions is already public
        return self.peripheral_ions

    def get_central_site(self): #TODO: a getter is NOT needed. The var peripheral_ions is already public
        return self.central_ion

    def set_central_name(self, new_name): #TODO: Why would someone wnat to rename the central_ion? Especially *after* constructing the object? Mutability is usually evil.
        self.central_ion_name = new_name

    def get_connection(self, other):  # TODO: probably rename this to get_num_connections or something
        """
        Gives the connectivity between the given polyhedra and another polyhedra by counting the number of atoms shared
        between the two polyhedra

        :param other: (Polyhedra) target Polyhedra against which we are checking for connectivity with the current
        polyhedra
        :return: (int) Integer giving the number of peripheral ions shared by the current polyhedra and the target
        polyhedra; returns -1 if the same polyhedra is being compared to itself
        """
        assert isinstance(other, Polyhedra)

        if self.central_ion == other.central_ion:
            return -1  # TODO: better to raise a ValueError than return some nonsensical result. The user should decide how to handle errors, not you. Please see "Writing Better Code", Tip #9 (newly added) for more details

        return len(self.get_connections(other))  # TODO: this is the only line needed in this function. The other lines above are not needed, they are repeated in the function call.

    def get_connections(self, other):
        """
        Gives the shared sites between the current Polyhedra and the given Polyhedra

        :param polyhedra2: (Polyhedra) target Polyhedra against which we are checking for connectivity with the current
        polyhedra
        :return: (list) list of sites shared between the current Polyhedra and the given Polyhedra
        """

        assert isinstance(other, Polyhedra)  # TODO: this is difficult to explain, but typically you don't need this. Python also encourages duck typing rather than assertions of type

        if self.central_ion == other.central_ion:  # TODO: it is not clear that this check will return True for two identical but different Polyhedra objects. You would need to first implement an __eq__ method.
            print "Checking connections between exact same polyhedra"

        shared_sites = []

        for site in other.peripheral_ions:
            for c_site in self.peripheral_ions:
                if c_site == site:
                    shared_sites.append(site)
        return shared_sites

    def is_connected(self, other):  # TODO: I don't think this is needed at all. You can just check the length of get_connection or even bool it. e.g., bool(get_connection)
        """
        Checks whether the current Polyhedra is connected to the given Polyhedra
        :param polyhedra2: (Polyhedra) target Polyhedra against which we are checking for connectivity with the current
        polyhedra
        :return: (boolean) Whether the two Polyhedra have any shared peripheralIons
        """

        assert isinstance(other, Polyhedra)

        if self.central_ion == other.central_ion:
            print "Checking connections between exact same polyhedra"

        for site in other.peripheralIons:
            for c_site in self.peripheralIons:
                if c_site == site:
                    return True
                    break

        return False

    def get_peripheral_distances(self):  # TODO: not clear this needs to be a separate function, unless it is used in multiple place. Use as many methods/abstractions/functions as necessary, but don't use more than that.
        """
        Gives the distances between the central site and the peripheral ion sites (for determining bond weights)

        :return: (list) list of distances between the central ion site and the peripheral ion sites
        """
        peripheral_distances = []
        for peripheral in self.peripheral_ions:
            peripheral_distances.append(peripheral.distance(self.central_ion))

        return peripheral_distances

    def __str__(self):
        return self.composition.formula


def get_connection(polyhedra1, polyhedra2):  # TODO: this function is completely pointless
    """
    Determines whether two polyhedra are connected (i.e. if they share any peripheral ions)

    :param polyhedra1: (Polyhedra) target polyhedra
    :param polyhedra2: (Polyhedra) target polyhedra we're getting connectivity with polyhedra1
    :return: (int) Integer giving the number of peripheral ions shared by the current polyhedra and the target
    polyhedra; returns -1 if the same polyhedra is being compared to itself
    """

    return polyhedra1.get_connectivity(polyhedra2)

    return connections  # TODO: you can't return after a return statement anyway


def get_surrounding_connectivity(structure, polyhedra, radius):
    """
    Gives the surrounding connectivity of a specific polyhedra in a given structure

    :param structure: (Structure) target structure
    :param polyhedra: (Polyhedra) target polyhedra
    :return: (list) list of polyhedra and their connections to the target polyhedra in the target structure, elements
    of the list are of the form: [[polyhedra, connection], ...]
    """

    polyhedraList = get_supercell_polyhedra(structure, radius)

    connectedPolyhedra = []
    for otherPolyhedra in polyhedraList:
        connection = polyhedra.get_connection(otherPolyhedra)
        if connection > 0:
            connectedPolyhedra.append([otherPolyhedra, connection])

    return connectedPolyhedra

def get_polyhedra(structure, site, radius=3.2):
    """
    Gives a polyhedra in a structure based on the site, where nearby atoms are considered peripheral ions to the site
    if their bond weights (as defined by Hoppe, 1979) are greater than given value (i.e. the atoms contribute to the
    ECoN)

    :param structure: (Structure) target structure
    :param site: (PeriodicSite) target site
    :param radius: (float) largest possible distance within which peripheral ions can be obtained
    :return: (Polyhedra) polyhedra object representing the polyhedra around the target site in the target structure
    """
    anions = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S']  # this should be a tuneable parameter

    anionSites = []
    bondlengths = []
    peripheralSites = []

    for entry in structure.get_neighbors(site, radius):
        if entry[0].species_string in anions and entry[1] < radius:
            anionSites.append(entry[0])
            bondlengths.append(entry[1])
    for potentialPeripheralIons in anionSites:
        #do not count nearby anions that do not contribute
        if effective_coordination.calculate_bond_weight(potentialPeripheralIons[1], bondlengths) > 1e-2:
            peripheralSites.append(potentialPeripheralIons)

    return Polyhedra(site, peripheralSites)



def count_connectivity(structure, radius = 3.2):

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
        site_coords = polyhedra.get_central_site().frac_coords
        if site_coords[0] >= (1.0/3)-1e-5 and site_coords[0] < (2.0/3):
            if site_coords[1] >= (1.0/3)-1e-5 and site_coords[1] < (2.0/3):
                if site_coords[2] >= (1.0/3)-1e-5 and site_coords[2] < (2.0/3):
                    center_cell_polyhedra.append(polyhedra)


    # Count the number of connections between the
    connections = {}
    for polyhedra1 in center_cell_polyhedra:
        cation = polyhedra1.get_central_name()
        connected = False #keep track of whether a polyhedra is not connected to any other polyhedra
        if not cation in connections.keys():
            connections[cation] = {"none": 0, "point": 0,"edge": 0, "face": 0}
        for polyhedra2 in polyhedra_list:
            connection = polyhedra1.get_connection(polyhedra2)
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

def get_supercell_polyhedra(structure, radius = 3.2):
    """
    Obtains all cation polyhedra needed to describe connectivity by creating a supercell and converting sites in the
    structure to cation polyhedra (peripheral ions of the polyhedra given by the

    :param structure: (Structure) target structure
    :param radius: (float) radius within which to determine whether a nearby atom is a peripheral ion
    :return: (list) list of Polyhedra objects in supercell large enough to determine connectivities between polyhedra
    """

    anions = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S']

    structure.make_supercell((3, 3, 3))

    polyhedra = []
    for site_index in range(structure.num_sites):
        if structure[site_index].species_string not in anions:
            polyhedra.append(get_polyhedra(structure, structure[site_index]))

    return polyhedra


def get_connectivity_matrix(structure, radius = 3.2):
    """
    Creates a connectivity matrix to describe the connectivity between cations in a structure; connections between
    cation polyhedra and reflections of itself (and reflections of other cation polyhedra) are also counted; different
    sites of the same cation species are NOT differentiated (that is, all instances of each cation are counted all
    together)

    :param structure: (Structure) target structure
    :param radius: (float) radius within which to determine whether a nearby atom is a peripheral ion
    :return: (dict) dictionary of dictionaries, with the first set of keys being the specified (unreflected) Polyhedra
    sites, the second (inner) set of keys being the specified (unreflected) Polyhedra sites which we are comparing
    connectivity to, and the values being a list [x, y, z] which give the numbers of [point-sharing, edge-sharing, and
    face-sharing] instances between the first Polyhedra and all of the reflections of the second Polyhedra
    """

    polyhedra_list = get_supercell_polyhedra(structure, radius)

    # get the polyhedra in the central cell of the supercell from which to calculate connectivities
    center_cell_polyhedra = []
    cation_names = []
    for polyhedra in polyhedra_list:
        site_coords = polyhedra.get_central_site().frac_coords
        if site_coords[0] >= (1.0/3)-1e-5 and site_coords[0] < (2.0/3)-1e-5:
            if site_coords[1] >= (1.0/3)-1e-5 and site_coords[1] < (2.0/3)-1e-5:
                if site_coords[2] >= (1.0/3)-1e-5 and site_coords[2] < (2.0/3)-1e-5:
                    center_cell_polyhedra.append(polyhedra)
                    cation_names.append(polyhedra.get_central_name())

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
        inner_cation = inner_polyhedra.get_central_name()
        for outer_polyhedra in polyhedra_list:
            outer_cation = outer_polyhedra.get_central_name()
            connection = inner_polyhedra.get_connection(outer_polyhedra)
            if connection < 0:
                continue
            if connection == 1:
                connections[inner_cation][outer_cation]["point"] += 1
            if connection == 2:
                connections[inner_cation][outer_cation]["edge"] += 1
            if connection >= 3:
                connections[inner_cation][outer_cation]["face"] += 1
        break

    return connections

def get_connectivity_matrix_2(structure, radius = 3.2):

    """
    Creates a connectivity matrix to describe the connectivity between cations in a structure; connections between
    cation polyhedra and reflections of itself (and reflections of other cation polyhedra) are also counted; different
    sites of the same species ARE differentiated (that is, connectivity for cations with different sites, not including
    reflections, are counted separately)

    :param structure: (Structure) target structure
    :param radius: (float) radius within which to determine whether a nearby atom is a peripheral ion
    :return: (dict) dictionary of dictionaries, with the first set of keys being the specified (unreflected) Polyhedra
    sites, the second (inner) set of keys being the specified (unreflected) Polyhedra sites which we are comparing
    connectivity to, and the values being a list [x, y, z] which give the numbers of [point-sharing, edge-sharing, and
    face-sharing] instances between the first Polyhedra and all of the reflections of the second Polyhedra
    """

    polyhedra_list = get_supercell_polyhedra(structure, radius)

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
        site_coords = polyhedra.get_central_site().frac_coords
        if site_coords[0] >= (1.0/3)-1e-5 and site_coords[0] < (2.0/3)-1e-5:
            if site_coords[1] >= (1.0/3)-1e-5 and site_coords[1] < (2.0/3)-1e-5:
                if site_coords[2] >= (1.0/3)-1e-5 and site_coords[2] < (2.0/3)-1e-5  :
                    if polyhedra.central_ion_name not in cation_iterators.keys():
                        cation_iterators[polyhedra.central_ion_name] = 1
                    else:
                        cation_iterators[polyhedra.central_ion_name] += 1
                    polyhedra.set_central_name(polyhedra.central_ion_name +
                                               str(cation_iterators[polyhedra.central_ion_name]))
                    center_cell_polyhedra.append(polyhedra)

    # classify supercell polyhedra into center-cell polyhedra by checking whether they are images
    for outer_polyhedra in polyhedra_list:
        for inner_polyhedra in center_cell_polyhedra:
            if check_image_in_supercell(outer_polyhedra.get_central_site(), inner_polyhedra.get_central_site()):
                outer_polyhedra.set_central_name(inner_polyhedra.get_central_name())
                break
        print outer_polyhedra.get_central_name()

    # instantiate nested dictionaries for connections
    connections = {}
    for cation_1 in center_cell_polyhedra:
        connections[cation_1.central_ion_name] = {}
        for cation_2 in center_cell_polyhedra:
            connections[cation_1.central_ion_name][cation_2.central_ion_name] = {"point": 0, "edge": 0, "face": 0}

    # count connections between specific polyhedra to create adjacency matrix
    # matrix should be symmetric, i.e. connections between cation1 and cation2 should be the same as connections
    # between cation2 and cation1
    for inner_polyhedra in center_cell_polyhedra:
        inner_cation = inner_polyhedra.get_central_name()
        for outer_polyhedra in polyhedra_list:
            outer_cation = outer_polyhedra.get_central_name()
            connection = inner_polyhedra.get_connection(outer_polyhedra)
            if connection < 0:
                continue
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
    s = Structure.from_file('Li5CoO4.cif', True, False)  # TODO: don't choose an 80-atom structure that kills my computer as the test case. Let the user run a simple BCC Fe or something.
    print s
    print get_connectivity_matrix_2(s, 3.0)  # TODO: this is completely eating my CPU. It prints nonsense like "Co8" repeated multiple times. The final output shows zero connections for almost all sites?

