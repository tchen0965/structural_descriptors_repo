import numpy as np
from pymatgen import Structure


__author__ = 'Tina'


# returns if there is a numMatches number of matches between list1 and list2
def check_matches(list1, list2, num_matches):
    """
    Check whether two lists of sites contain at least a certain number of matches between the two lists
    Args:
    :param list1: (List) list of sites
    :param list2: (List) list of sites
    :param num_matches: (int) number of matches we are looking for between list1 and list2
    :return: (boolean) whether there are numMatches of matches between list1 and list2
    """
    intersection = set(list1).intersection(set(list2))
    if len(intersection) >= num_matches:
        return True
    else:
        return False


def check_match(list1, list2):
    """
    Check whether two lists of sites contain the same sites between the two lists
    Args:
    :param list1: (List) list of sites
    :param list2: (List) list of sites
    :return: (boolean) whether list1 and list2 have the same sites within them
    """
    if set(list1) == set(list2):
        return True
    else:
        return False


def get_anion_neighbors(site, structure, radius, anions, get_distance=False):
    """
    Gets neighboring anions of sites
    Args:
    :param site: (Site) target site to find anion neighbors of
    :param structure: (Structure) structure that contains target site
    :param radius: (float) radius to which anion neighbors should be looked for
    :param anions: (List of Strings) list of species considered anions
    :param get_distance: (boolean) whether or not to get distance between cation and anion
    :return: (List of either Sites or [Site, float]) list of either anion neighboring sites or [Site, float] if
    getDistance is True
    """
    anion_neighbors = []
    neighbors = structure.get_neighbors(site, radius)
    for neighbor in neighbors:
        if neighbor[0].species_string in anions and neighbor[1] < radius:
            if get_distance:
                anion_neighbors.append(neighbor)
            else:
                anion_neighbors.append(neighbor[0])
    return anion_neighbors


def split_sites(structure, voronoi_point):
    """
    Split sites in a structure into sites that are Voronoi points and sites containing real species

    Args:
    :param structure: (Structure) target structure
    :return: (List of Sites, List of Sites) list of sites with Voronoi points, list of sites with real species
    """
    all_voronois = []
    other_sites = []
    for site in structure.sites:
        if site.species_string == voronoi_point:
            all_voronois.append(site)
        else:
            other_sites.append(site)
    return all_voronois, other_sites


def get_best_voronoi_from_dict(indices, structure, sites, radius, anions):
    """

    Args:
    :param indices: (list) list of indices
    :param structure: (Structure) target structure
    :param sites: (dict) dictionary of sites with indices as the key
    :return: index from indices associated with the site (in sites) with the most even bonds
    """
    total_errors = {}
    for i in indices:
        site = sites[i]
        neighbors = get_anion_neighbors(site, structure, radius, anions, get_distance=True)
        bond_lengths = []
        errors = []
        for neighbor in neighbors:
            bond_lengths.append(neighbor[1])
        for bond in bond_lengths:
            errors.append(np.average(bond_lengths) - bond)
        total_errors[site] = sum(errors)
    return min(indices, key=lambda x: total_errors[sites[x]])


def structure_from_sites(voronoi_points, other_sites, structure):
    """
    Creates a new structure with the sites from both the list of voronoi points and the sites from other_sites
    Args:
    :param voronoi_points: (List of Sites) list of sites with Voronoi points
    :param other_sites: (List of Sites) list of sites with real species on the sites
    :param structure: (Structure) structure with target lattice
    :return: (Structure) structure with both Voronoi sites and sites with real species
    """
    species = []
    frac_coords = []
    for site in voronoi_points:
        species.append(site.species_string)
        frac_coords.append(site.frac_coords)
    for site in other_sites:
        species.append(site.species_string)
        frac_coords.append(site.frac_coords)

    return Structure(structure.lattice, species, frac_coords)