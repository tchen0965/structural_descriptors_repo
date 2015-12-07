__author__ = 'Tina'

from pymatgen.core.structure import Structure
from pymatgen import PeriodicSite
from substructure import substructures_from_structure
import math
import re
import numpy


def getAvgOKeeffeCN(structure, radius):
    substructList = substructures_from_structure(structure)
    cationList = []
    listCN = {}
    for substruct in substructList:
        cation = substruct.central_ion
        if not cation in cationList:
            cationList.append(cation)
            listCN[cation] = [substruct.weight_sum()]
        else:
            listCN[cation].append(substruct.weight_sum())
    avgCNs = []
    for key in listCN.keys():
        avgCNs.append([key, numpy.mean(listCN[key])])
    return avgCNs

if __name__ == '__main__':
    s = Structure.from_file('LiCoO2.cif', True, True)
    print getAvgOKeeffeCN(s, 3.0)
