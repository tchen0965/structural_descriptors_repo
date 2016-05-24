from __future__ import division
import math
import numpy


__author__ = 'Tina_Chen'
__contributor__ = 'Anubhav_Jain'


class EffectiveCoordFinder(object):

    """

    Finds the average effective coordination number for each cation in a given structure. It
    finds all cation-centered polyhedral in the structure, calculates the bond weight for each peripheral ion in the
    polyhedral, and sums up the bond weights to obtain the effective coordination number for each polyhedral. It then
    averages the effective coordination of all polyhedral with the same cation at the central site.

    We use the definition from Hoppe (1979) to calculate the effective coordination number of the polyhedrals:

    Hoppe, R. (1979). Effective coordination numbers (ECoN) and mean Active fictive ionic radii (MEFIR).
    Z. Kristallogr. , 150, 23-52.
    ECoN = sum(exp(1-(l_i/l_av)^6)), where l_av = sum(l_i*exp(1-(1_i/l_min)))/sum(exp(1-(1_i/l_min)))

    """

    def __init__(self, structure):
        self._structure = structure

    def get_site_cn(self, site, radius=2.6, min_weight=10e-5):

        neighboring_sites = []
        anion_sites = []  # list of all neighboring anions
        bond_lengths = []  # list of bond lengths of anions
        bond_weights = []  # list of bond weights of anions
        for entry in self._structure.get_neighbors(site, radius):  # entry = (site, distance)
            anion_sites.append(entry)
            bond_lengths.append(entry[1])

        for bond in anion_sites:
            # do not count nearby anions that do not contribute
            if calculate_bond_weight(bond[1], bond_lengths) > min_weight:
                neighboring_sites.append(bond)
                bond_weights.append(calculate_bond_weight(bond[1], bond_lengths))

        return sum(bond_weights)

    def get_isite_cn(self, isite, radius=2.6, min_weight=10e-5):
        return self.get_site_cn(self._structure[isite], radius, min_weight)

    def get_site_neighbors(self, site, radius=2.6, min_weight=10e-5):

        neighboring_sites = []
        anion_sites = []  # list of all neighboring anions
        bond_lengths = []  # list of bond lengths of anions
        bond_weights = []  # list of bond weights of anions
        for entry in self._structure.get_neighbors(site, radius):  # entry = (site, distance)
            anion_sites.append(entry)
            bond_lengths.append(entry[1])

        for bond in anion_sites:
            # do not count nearby anions that do not contribute
            if calculate_bond_weight(bond[1], bond_lengths) > min_weight:
                neighboring_sites.append(bond)
                bond_weights.append(calculate_bond_weight(bond[1], bond_lengths))

        return neighboring_sites

    def get_avg_cn(self, radius=2.6, anions=None):
        """
        Get the average coordination for all cations in structure.

        :param radius: (float) distance in Angstroms for bond cutoff
        :param anions: (List of Strings) list of species which we consider anions in the structure
        :return: (dict) A dictionary with keys corresponding to different cations and the values to the cation's ECoN
            coordination number averaged over all polyhedra with the same cation center in the structure
        """

        if anions is None:
            anions = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S']

        cation_cns = self.get_cation_cn(radius, anions)

        avg_cns = {}
        for cation in cation_cns.keys():
            avg_cns[cation] = numpy.mean(cation_cns[cation])

        return avg_cns

    def get_cation_cn(self, radius=2.6, min_weight=10e-5, anions=None):
        """
        Get all cation-centered polyhedra for a structure

        :param radius: (float) distance in Angstroms for bond cutoff
        :param anions: (List of Strings) list of species which we consider anions in the structure
        :return: (dict) A dictionary with keys corresponding to different cations and the values to the cation's
            ECoN coordination numbers
        """

        if anions is None:
            anions = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S']

        # pick out all sites with cations (not just "cation") as the species at the site
        cation_sites = []
        for site in self._structure.sites:
            if site.species_string not in anions:  # TODO: how do you handle mixed occupancy sites?
                cation_sites.append(site)
        # TODO: PolyhedralList needs to be a more organized/encapsulated data structure. Talk to Anubhav
        cn_list = {}  # list of polyhedrals with cations at the center
        for site in cation_sites:

            site_cn = self.get_site_cn(site, radius, min_weight)

            if site.species_string not in cn_list.keys():
                cn_list[site.species_string] = [site_cn]
            else:
                cn_list[site.species_string].append(site_cn)

        return cn_list


def get_effective_cn(polyhedra):
    """
    Get the effective coordination number for a Polyhedra object

    :param polyhedra: target polyhedra
    :return: (float) Effective coordination number of the polyhedra as given by Hoppe, 1979
    """
    bond_weights = []
    bond_lengths = []
    for peripheral in polyhedra.peripheral_ions:
        bond_lengths.append(peripheral.distance(polyhedra.central_ion))
    for peripheral in polyhedra.peripheral_ions:
        bond_weight = calculate_bond_weight(peripheral.distance(polyhedra.central_ion), bond_lengths)
        bond_weights.append(bond_weight)

    return sum(bond_weights)


def calculate_bond_weight(bond_length, bonds):
    """
    Get the weight of a given bond (defined by the effective coordination number formula in Hoppe, 1979) in a polyhedra

    :param bond_length: (float) target bond length
    :param bonds: (list) list of all bond lengths between the central ion and the peripheral ions in a polyhedra
    (including target bond)
    :return: (float) the bond's weight
    """

    weight = math.exp(1-(bond_length/calculate_weighted_avg(bonds))**6)
    return weight


def calculate_weighted_avg(bonds):
    """
    Get the weighted average bond length given by the effective coordination number formula in Hoppe (1979)

    :param bonds: (list) list of floats that are the bond distances between a cation and its peripheral ions
    :return: (float) exponential weighted average
    """

    minimum_bond = min(bonds)
    weighted_sum = 0.0
    total_sum = 0.0
    for entry in bonds:
        weighted_sum += entry*math.exp(1 - (entry/minimum_bond)**6)
        total_sum += math.exp(1-(entry/minimum_bond)**6)
    return weighted_sum/total_sum
