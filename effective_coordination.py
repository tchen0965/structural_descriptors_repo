from __future__ import division

__author__ = 'Tina_Chen'
__contributor__ = 'Anubhav_Jain'

import math
import numpy
from connectivity_from_structure import Polyhedra


class EffectiveCoordFinder(object):

    """

    Finds the average effective coordination number for each cation in a given structure. It
    finds all cation-centered polyhedral in the structure, calculates the bond weight for each peripheral ion in the
    polyhedral, and sums up the bond weights to obtain the effective coordination number for each polyhedral. It then
    averages the effective coordination of all polyhedral with the same cation at the central site.

    We use the definition from Hoppe (1979) to calculate the effective coordination number of the polyhedrals:

    Hoppe, R. (1979). Effective coordination numbers (ECoN) and mean Active fictive ionic radii (MEFIR). Z. Kristallogr. , 150, 23-52.
    ECoN = sum(exp(1-(l_i/l_av)^6)), where l_av = sum(l_i*exp(1-(1_i/l_min)))/sum(exp(1-(1_i/l_min)))

    """

    def __init__(self, structure):
        self._structure = structure
        #if target is None:
        #    self._target = structure.composition.elements
        #else:
        #    self._target = target

    def get_avg_CN(self, radius=3.0, anions = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S']):
        """
        Get the average coordination for all cations in structure.

        :param structure: (Structure) target structure
        :param radius: (float) distance in Angstroms for bond cutoff
        :return: (dict) A dictionary with keys corresponding to different cations and the values to the cation's ECoN
            coordination number averaged over all polyhedra with the same cation center in the structure
        """
        cationCNs = self.get_cation_CN(radius, anions)

        avgCNs = {}
        for cation in cationCNs.keys():
            avgCNs[cation] = numpy.mean(cationCNs[cation])

        return avgCNs


    def get_cation_CN(self, radius = 3.0, anions = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S']):
        """
        Get all cation-centered polyhedra for a structure

        :param structure: (Structure) target structure
        :param radius: (float) distance in Angstroms for bond cutoff
        :return: (dict) A dictionary with keys corresponding to different cations and the values to the cation's
            ECoN coordination numbers
        """

        # pick out all sites with cations (not just "cation") as the species at the site
        cationSites = []
        for site in self._structure.sites:
            if not site.species_string in anions:  # TODO: how do you handle mixed occupancy sites?
                cationSites.append(site)
        CNList = {} #list of polyhedrals with cations at the center  # TODO: PolyhedralList needs to be a more organized/encapsulated data structure. Talk to Anubhav
        for site in cationSites:
            sites = {}
            sites[site.species_string] = []
            anionSites = [] # list of all neighboring anions
            bondlengths = [] # list of bond lengths of anions
            bondweights = [] # list of bond weights of anions
            for entry in self._structure.get_neighbors(site, radius): #entry = (site, distance)
                if entry[0].species_string in anions and entry[1] < radius:
                    anionSites.append(entry)
                    bondlengths.append(entry[1])

            for bond in anionSites:
                if calculate_bond_weight(bond[1], bondlengths) > 10e-5: #do not count nearby anions that do not contribute
                    sites[site.species_string].append(bond)
                    bondweights.append(calculate_bond_weight(bond[1], bondlengths))

            if not site.species_string in CNList.keys():
                CNList[site.species_string] = [sum(bondweights)]
            else:
                CNList[site.species_string].append(sum(bondweights))

        return CNList

def get_effective_CN(polyhedra):
    """
    Get the effective coordination number for a Polyhedra object

    :param polyhedra: target polyhedra
    :return: (float) Effective coordination number of the polyhedra as given by Hoppe, 1979
    """
    bondweights = []
    bondlengths = []
    for peripheral in polyhedra.peripheral_ions:
        bondlengths.append(peripheral.distance(polyhedra.central_ion))
    for peripheral in polyhedra.peripheral_ions:
        bondweight = calculate_bond_weight(peripheral.distance(polyhedra.central_ion), bondlengths)
        bondweights.append(bondweight)

    return sum(bondweights)

def calculate_bond_weight(bondlength, bonds):
    """
    Get the weight of a given bond (defined by the effective coordination number formula in Hoppe, 1979) in a polyhedra

    :param bond: (float) target bond length
    :param bonds: (list) list of all bond lengths between the central ion and the peripheral ions in a polyhedra
    (including target bond)
    :return: (float) the bond's weight
    """

    weight = math.exp(1-(bondlength/calculate_weighted_avg(bonds))**6)
    return weight


def calculate_weighted_avg(bonds):
    """
    Get the weighted average bond length given by the effective coordination number formula in Hoppe (1979)

    :param bonds: (list) list of floats that are the bond distances between a cation and its peripheral ions
    :return: (float) exponential weighted average
    """

    minimumBond = min(bonds)
    weightedSum = 0.0
    totalSum = 0.0
    for entry in bonds:
        weightedSum += entry*math.exp(1 - (entry/minimumBond)**6)
        totalSum += math.exp(1-(entry/minimumBond)**6)
    return weightedSum/totalSum