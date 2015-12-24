__author__ = 'Tina'
__contributor__ = 'Anubhav Jain'

from pymatgen.core.structure import Structure
import math
import numpy

#TODO: the file name "effective_coordination_number_from_structure.py" is way too long. you can do "effective_coordination.py"

"""
TODO: provide docs that give an overview of this file. In particular,
you might want to give a reference to a paper or web site regarding what is the effective coordination number
"""


def getAvgCN(structure, radius):
    """
    # TODO: Tina, please use this documentation example for how to write docs (delete this TODO when read)
    # TODO: Return a dictionary, not an awkward list of lists

    Get the average coordination for all cations in structure.

    :param structure: (Structure) target structure
    :param radius: (float) distance in Angstroms for bond cutoff
    :return: ([[(str), (float)]]) A list of entries, each entry being an array containing cation and coordination #
    """
    polyhedraList = getAllCationPolyhedral(structure, radius)
    cationList = []
    listCN = {}
    for polyhedra in polyhedraList:
        cation = polyhedra[0][0].species_string  # TODO: the [0][0] is confusing. Underlying data structure should be fixed
        if not cation in cationList:  # TODO: again, replace with dict or defaultdict
            cationList.append(cation)
            listCN[cation] = [polyhedra[1]]
        else:
            listCN[cation].append(polyhedra[1])
    avgCNs = []
    for key in listCN.keys():
        avgCNs.append([key, numpy.mean(listCN[key])])
    return avgCNs


def getAllCationPolyhedral(structure, radius = 3.0):
    # TODO: add docs (see example in previous method)

    anions = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S']

    # pick out all sites with cations (not just "cation") as the species at the site
    cationSites = []
    for site in structure.sites:
        if not site.species_string in anions:  # TODO: how do you handle mixed occupancy sites?
            cationSites.append(site)
    #print structure.formula
    polyhedralList = [] #list of polyhedrals with cations at the center  # TODO: PolyhedralList needs to be a more organized/encapsulated data structure. Talk to Anubhav
    for site in cationSites:
        sites = [site] #list of sites that will be in polyhedral
        anionSites = [] # list of all neighboring anions
        bondlengths = [] # list of bond lengths of anions
        bondweights = [] # list of bond weights of anions
        for entry in structure.get_neighbors(site, radius):
            if entry[0].species_string in anions and entry[1] < radius:
                anionSites.append(entry)
                bondlengths.append(entry[1])

        for bond in anionSites:
            if site.species_string == 'Rn':
                sites.append(bond)
                bondweights.append(1.0)
            else:
                if calculateBondWeight(bond[1], bondlengths) > 0: #do not count nearby anions that do not contribute
                    sites.append(bond)
                    bondweights.append(calculateBondWeight(bond[1], bondlengths))
        polyhedralList.append([sites, round(sum(bondweights))])

    return polyhedralList


def calculateBondWeight(bond, bonds):  # TODO: use PEP naming, e.g. calculate_bond_weight
    # TODO: add docstrings

    weight = math.exp(1-(bond/calculateWeightedAvg(bonds))**6)
    return weight

# calculates the weighted average bond length
# takes a list of bond lengths as the input
def calculateWeightedAvg(bonds):
    # TODO: move docs inside the function, see docstring example

    minimumBond = min(bonds)
    weightedSum = 0.0
    totalSum = 0.0
    for entry in bonds:
        weightedSum += entry*math.exp(1 - (entry/minimumBond)**6)
        totalSum += math.exp(1-(entry/minimumBond)**6)
    return weightedSum/totalSum  # TODO: although you are using floats and this code is fine, you should add "from __future__ import division" to top of code.

if __name__ == '__main__':
    # TODO: move this code to an "examples" module. Include the file LiCoO2.cif in the repo in that module, otherwise no one else can run your code.
    s = Structure.from_file('LiCoO2.cif', True, True)
    print s
    print getAvgCN(s, 3.2)