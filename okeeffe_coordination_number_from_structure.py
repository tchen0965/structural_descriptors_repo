from __future__ import division

__author__ = 'Tina'

from substructure import substructures_from_structure
import numpy

"""

Finds the average O'Keeffe (Voronoi) coordination number for each cation in a given structure. Finds all
cation-centered polyhedral in the structure, calculates the bond weight for each peripheral ion in the
polyhedral, and sums up the bond weights to obtain the effective coordination number for each polyhedral. It then
averages the effective coordination of all polyhedral with the same cation at the central site.

TO-DO: Change method of getting Voronoi structures from substructures_from_structure to Pymatgen's
VoronoiCoordFinder (need to figure out how to import from Pymatgen for Github)

"""


def getAvgCN(structure):
    """

    Get average O'Keeffe coordination number for all cations in structure

    :param structure:(Structure) target structure
    :return: (dict) A dictionary with keys corresponding to different cations and the values to the cation's O'Keeffe
        coordination number averaged over all polyhedra with the same cation center in the structure
    """
    substructList = substructures_from_structure(structure)
    cationsAlreadySeen = {}
    listCN = {}
    for substruct in substructList:
        cation = substruct.central_ion
        if not cation in cationsAlreadySeen.keys():
            cationsAlreadySeen[cation] = 0
            listCN[cation] = [substruct.weight_sum()]
        else:
            listCN[cation].append(substruct.weight_sum())
    avgCNs = {}
    for cation in listCN.keys():
        avgCNs[cation] = numpy.mean(listCN[cation])
    return avgCNs