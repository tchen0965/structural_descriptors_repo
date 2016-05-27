'''
Created on Mar 12, 2015

@author: orvanano
'''
'''
Created on May 7, 2014
outputs the structure with the voronoi points
'''
import random

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from voronoi import VoronoiInsertionTransformation
from pymatgen import Structure

#structure1 = Structure, name = string
def add_voronoi(structure):

    """
    Args:
    :param structure1: (Structure) target structure
    :return structure3: (Structure) structure with all Voronoi points added to the structure
    """

    #choose an element that is not present in any of the structures to use to identify the voronoi sites
    transformer = VoronoiInsertionTransformation("Rn", midpoints=True)
    structure_with_voronoi = transformer.apply_transformation(structure)

    findsim=SpacegroupAnalyzer(structure_with_voronoi, symprec=1e-1)
    symmetrized_structure_with_voronoi = findsim.get_symmetrized_structure()

    return symmetrized_structure_with_voronoi

def add_voronoi_to_file(structure, name):

    """
    Args:
    :param structure1: (Structure) target structure
    :return structure3: (Structure) structure with all Voronoi points added to the structure
    """

    #choose an element that is not present in any of the structures to use to identify the voronoi sites
    structure_with_voronoi = add_voronoi(structure)

    structure_with_voronoi.to(filename="%s.cif"%(name))


if __name__ == '__main__':
    cation = 'Li'
    anion = 'O'

    #filename="structures_for_empty_polyhedra/%s.cif"%(task_id)
    s1 = Structure.from_file("LiCoO2.cif", True, True)
    add_voronoi_to_file(s1, "LiCoO2_with_voronoi")

    """

    print "getting data from database"
    rawdata = getEntries(cation, anion)

    print "picking random data 500 structures out of 10,000"
    data = {}
    x = 0
    while x < 500:
        task_id = random.choice(rawdata.keys())
        if not task_id in data.keys():
            data[task_id] = rawdata[task_id]
            x += 1
            print "Obtaining data point:", x

    print "outputting structures to folder structures_for_empty_polyhedra"
    x = 0
    for task_id in data.keys():
        data[task_id][0].to(filename="structures_for_empty_polyhedra/%s.cif"%(task_id))
        print "Printing original structure", x, ":", task_id
        x += 1

    x = 0
    print "adding voronoi points to structures"
    for task_id in data.keys():
        addVoronoi(data[task_id][0], str(task_id))
        print "Adding voronoi", x, ":", task_id
        x += 1

    """