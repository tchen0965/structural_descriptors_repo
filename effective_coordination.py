from __future__ import division

__author__ = 'Tina'
__contributor__ = 'Anubhav Jain'

import math
import numpy


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

    def __init__(self, structure, target=None):
        self._structure = structure
        if target is None:
            self._target = structure.composition.elements
        else:
            self._target = target

    def getAvgCN(self, radius):
        """

        Get the average coordination for all cations in structure.

        :param structure: (Structure) target structure
        :param radius: (float) distance in Angstroms for bond cutoff
        :return: (dict) A dictionary with keys corresponding to different cations and the values to the cation's ECoN
            coordination number averaged over all polyhedra with the same cation center in the structure
        """
        cationCNs = self.get_all_cation_polyhedral(radius)

        avgCNs = {}
        for cation in cationCNs.keys():
            avgCNs[cation] = numpy.mean(cationCNs[cation])

        return avgCNs


    def get_all_cation_polyhedral(self, radius = 3.0):
        """

        Get all cation-centered polyhedra for a structure

        :param structure: (Structure) target structure
        :param radius: (float) distance in Angstroms for bond cutoff
        :return: (dict) A dictionary with keys corresponding to different cations and the values to the cation's
            ECoN coordination numbers
        """

        anions = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S']

        # pick out all sites with cations (not just "cation") as the species at the site
        cationSites = []
        for site in self._structure.sites:
            if not site.species_string in anions:  # TODO: how do you handle mixed occupancy sites?
                cationSites.append(site)
        #print structure.formula
        polyhedralList = {} #list of polyhedrals with cations at the center  # TODO: PolyhedralList needs to be a more organized/encapsulated data structure. Talk to Anubhav
        for site in cationSites:
            sites = {}
            sites[site.species_string] = []
            anionSites = [] # list of all neighboring anions
            bondlengths = [] # list of bond lengths of anions
            bondweights = [] # list of bond weights of anions
            for entry in self._structure.get_neighbors(site, radius):
                if entry[0].species_string in anions and entry[1] < radius:
                    anionSites.append(entry)
                    bondlengths.append(entry[1])

            for bond in anionSites:
                if self.calculate_bond_weight(bond[1], bondlengths) > 10.0**-5: #do not count nearby anions that do not contribute
                    sites[site.species_string].append(bond)
                    bondweights.append(self.calculate_bond_weight(bond[1], bondlengths))

            if not site.species_string in polyhedralList.keys():
                polyhedralList[site.species_string] = [round(sum(bondweights))]
            else:
                polyhedralList[site.species_string].append(round(sum(bondweights)))

        return polyhedralList


    def calculate_bond_weight(self, bond, bonds):
        """

        Get the weight of a given bond (defined by the effective coordination number formula in Hoppe, 1979) in a polyhedra

        :param bond: (float) target bond
        :param bonds: (list) list of all bonds to peripheral ions in a polyhedra (including target bond)
        :return: (float) the bond's weight
        """

        weight = math.exp(1-(bond/self.calculate_weighted_avg(bonds))**6)
        return weight

    # calculates the weighted average bond length
    # takes a list of bond lengths as the input
    def calculate_weighted_avg(self, bonds):
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