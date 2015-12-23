__author__ = 'Tina'
__contributor__ = 'Anubhav Jain'

from pymatgen.core.structure import Structure
import math
import numpy

def getAvgCN(structure, radius):
    polyhedraList = getAllCationPolyhedral(structure, radius)
    cationList = []
    listCN = {}
    for polyhedra in polyhedraList:
        cation = polyhedra[0][0].species_string
        if not cation in cationList:
            cationList.append(cation)
            listCN[cation] = [polyhedra[1]]
        else:
            listCN[cation].append(polyhedra[1])
    avgCNs = []
    for key in listCN.keys():
        avgCNs.append([key, numpy.mean(listCN[key])])
    return avgCNs


def getAllCationPolyhedral(structure, radius = 3.0):

    anions = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S']

    # pick out all sites with cations (not just "cation") as the species at the site
    cationSites = []
    for site in structure.sites:
        if not site.species_string in anions:
            cationSites.append(site)
    #print structure.formula
    polyhedralList = [] #list of polyhedrals with cations at the center
    for site in cationSites:
        sites = [site] #list of sites that will be in polyhedral
        anionSites = [] # list of all neighboring anions
        bondlengths = [] # list of bond lengths of anions
        bondweights = [] # list of bond weights of anions
        for entry in structure.get_neighbors(site, radius):
            if entry[0].species_string in anions and entry[1] < radius:
                # sites.append(entry)
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
                    #print bond
        #print "next central site"
        polyhedralList.append([sites, round(sum(bondweights))])
        #polyhedralList.append(sites)
            # list of polyhedral (-> list of list of sites) + coordnum

    return polyhedralList


def calculateBondWeight(bond, bonds):
    weight = math.exp(1-(bond/calculateWeightedAvg(bonds))**6)
    return weight

# calculates the weighted average bond length
# takes a list of bond lengths as the input
def calculateWeightedAvg(bonds):
    minimumBond = min(bonds)
    weightedSum = 0.0
    totalSum = 0.0
    for entry in bonds:
        weightedSum += entry*math.exp(1 - (entry/minimumBond)**6)
        totalSum += math.exp(1-(entry/minimumBond)**6)
    return weightedSum/totalSum

if __name__ == '__main__':
    s = Structure.from_file('LiCoO2.cif', True, True)
    print getAvgCN(s, 3.0)