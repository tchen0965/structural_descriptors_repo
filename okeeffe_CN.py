from __future__ import division
from substructure import substructures_from_structure
import numpy

__author__ = 'Tina'


"""

Finds the average O'Keeffe (Voronoi) coordination number for each cation in a given structure. Finds all
cation-centered polyhedral in the structure, calculates the bond weight for each peripheral ion in the
polyhedral, and sums up the bond weights to obtain the effective coordination number for each polyhedral. It then
averages the effective coordination of all polyhedral with the same cation at the central site.

"""


def get_avg_cn(structure):
    """

    Get average O'Keeffe coordination number for all cations in structure

    :param structure:(Structure) target structure
    :return: (dict) A dictionary with keys corresponding to different cations and the values to the cation's O'Keeffe
        coordination number averaged over all polyhedra with the same cation center in the structure
    """
    substruct_list = substructures_from_structure(structure)
    cations_already_seen = {}
    list_cn = {}
    for substruct in substruct_list:
        cation = substruct.central_subspecies.specie.symbol
        if not cation in cations_already_seen.keys():
            cations_already_seen[cation] = 0
            list_cn[cation] = [substruct.weight_sum()]
        else:
            list_cn[cation].append(substruct.weight_sum())
    avgCNs = {}
    for cation in list_cn.keys():
        avgCNs[cation] = numpy.mean(list_cn[cation])
    return avgCNs