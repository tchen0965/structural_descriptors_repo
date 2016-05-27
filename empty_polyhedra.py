from pymatgen import PeriodicSite
from Voronoi_sites import add_voronoi
from empty_polyhedra_util import *

__author__ = 'Tina'

"""

Created March 17, 2015
Takes code with all Voronoi points inserted into structure (as cif file from Voronoi_sites) as input
These points are indicated by Rn, an element not present in any of the structures
    -Outputted as struct_with_voronoi.cif
Check to see if sites correspond to the center of an empty polyhedra (and removes Voronoi points/Rn elements from
    sites which are not empty polyhedra)
Remaining instances of Rn are those corresponding to an empty polyhedron (not necessarily at center)

"""


# remove duplicates of voronois
# remove voronois that have fewer than 4 neighboring oxygen atoms
def step1(all_voronois, structure, radius, anions, min_neighbors):
    """
    Remove sites with Voronoi points that have fewer than 4 neighboring anions

    :param all_voronois: (List of Sites) list of all sites in structure1 containing Voronoi points
    :param structure: (Structure) target structure
    :param radius: (float) maximum radius to look out to for neighbors
    :param min_neighbors: (int) minimum number of neighbors a point must have to still be considered
    :return: (List of Sites) list of sites with Voronoi points to keep (Voronoi sites with fewer than 4 neighboring
    anions removed)
    """
    new_voronois = []
    for site in all_voronois:

        # take out duplicates of voronois
        if site in new_voronois:
            continue

        # include only those sites surrounded by at least 3 (or some number) of peripheral ions
        neighbors = []
        cation_neighbors = []
        for entry in structure.get_neighbors(site, radius):
            if entry[0].species_string in anions and entry[1] < radius:
                neighbors.append(entry)
            else:
                if entry[1] < radius:
                    cation_neighbors.append(entry)
        if len(neighbors) < min_neighbors:
            continue

        new_voronois.append(site)

    return new_voronois


def step3(voronoi_points, structure, radius, anions, num_matches):
    """
    Removes Voronoi points from the list voronoi_points from the given the structure which are basically the same as
    another Voronoi

    Args:
    :param voronoi_points: (List of Sites) list of sites in structure containing candidate Voronoi points
    (empty polyhedra centers)
    :param structure: (Structure) target Structure
    :param radius: (float) distance to which we check for anions
    :param anions: (List of Strings) list of species which we consider anions
    :param num_matches: (int) number of matches between neighboring anions to consider a Voronoi site the same as
    another
    :return: remaining Voronoi points

    """
    key = {}
    similar_voronoi = {}  # sets of voronoi points that are the same
    keys_viewed = []
    to_remove = []
    for x in range(0, len(voronoi_points)):
        key[x] = voronoi_points[x]
    for num1 in key.keys():
        print "voronoi point being analyzed:", num1
        similar_voronoi[num1] = []
        if num1 in to_remove:
            continue
        for num2 in key.keys():
            if num1 == num2:
                continue
            if num2 in keys_viewed or num2 in to_remove:
                continue
            neighbors1 = get_anion_neighbors(key[num1], structure, radius, anions)
            neighbors2 = get_anion_neighbors(key[num2], structure, radius, anions)
            if check_match(neighbors1, neighbors2) or check_matches(neighbors1, neighbors2, num_matches):
                similar_voronoi[num1].append(num2)
                to_remove.append(num2)
            else:
                # PeriodicSite(atoms_n_occu, coords, lattice)
                x = 0
                y = 0
                z = 0
                reflections = []
                to_reflect = key[num1]
                if 0.0 <= to_reflect.frac_coords[0] <= 0.1:
                    #reflections.append(PeriodicSite(toReflect.species_string, toReflect.frac_coords + [1.0, 0, 0],
                    #                                toReflect.lattice))
                    x = 1
                if 0.9 <= to_reflect.frac_coords[0] <= 1.0:
                    #reflections.append(PeriodicSite(toReflect.species_string, toReflect.frac_coords + [-1.0, 0, 0],
                    #                                toReflect.lattice))
                    x = -1
                if 0.0 <= to_reflect.frac_coords[1] <= 0.1:
                    #reflections.append(PeriodicSite(toReflect.species_string, toReflect.frac_coords + [0, 1.0, 0],
                    #                                toReflect.lattice))
                    y = 1
                if 0.9 <= to_reflect.frac_coords[1] <= 1.0:
                    #reflections.append(PeriodicSite(toReflect.species_string, toReflect.frac_coords + [0, -1.0, 0],
                    #                                toReflect.lattice))
                    y = -1
                if 0.0 <= to_reflect.frac_coords[2] <= 0.1:
                    #reflections.append(PeriodicSite(toReflect.species_string, toReflect.frac_coords + [0, 0, 1.0],
                    #                                toReflect.lattice))
                    z = 1
                if 0.9 <= to_reflect.frac_coords[2] <= 1.0:
                    #reflections.append(PeriodicSite(toReflect.species_string, toReflect.frac_coords + [0, 0, -1.0],
                    #                                toReflect.lattice))
                    z = -1

                if not x == 0:
                    reflections.append(PeriodicSite(to_reflect.species_string, to_reflect.frac_coords + [x*1.0, 0, 0],
                                                    to_reflect.lattice))
                    if not y == 0:
                        reflections.append(PeriodicSite(to_reflect.species_string,
                                                        to_reflect.frac_coords + [x*1.0, y*1.0, 0], to_reflect.lattice))
                    if not z == 0:
                        reflections.append(PeriodicSite(to_reflect.species_string,
                                                        to_reflect.frac_coords + [x*1.0, 0, z*1.0], to_reflect.lattice))
                if not y == 0:
                    reflections.append(PeriodicSite(to_reflect.species_string, to_reflect.frac_coords + [0, y*1.0, 0],
                                                    to_reflect.lattice))
                    if not z == 0:
                        reflections.append(PeriodicSite(to_reflect.species_string,
                                                        to_reflect.frac_coords + [0, y*1.0, z*1.0], to_reflect.lattice))
                if not z == 0:
                    reflections.append(PeriodicSite(to_reflect.species_string, to_reflect.frac_coords + [0, 0, z*1.0],
                                                    to_reflect.lattice))
                    if not x == 0 and not y == 0:
                        reflections.append(PeriodicSite(to_reflect.species_string,
                                                        to_reflect.frac_coords + [x*1.0, y*1.0, z*1.0],
                                                        to_reflect.lattice))

                if len(reflections) > 0:
                    temp_species = []
                    temp_frac_coords = []
                    for site in structure.sites:
                        temp_species.append(site.species_string)
                        temp_frac_coords.append(site.frac_coords)
                    for reflection in reflections:
                        temp_species.append(reflection.species_string)
                        temp_frac_coords.append(reflection.frac_coords)
                    temp_struct = Structure(structure.lattice, temp_species, temp_frac_coords)
                    for reflection in reflections:
                        temp_neighbors1 = get_anion_neighbors(reflection, temp_struct, radius, anions)
                        if check_match(temp_neighbors1, neighbors2) or \
                                check_matches(temp_neighbors1, neighbors2, num_matches):
                            similar_voronoi[num1].append(num2)
                            to_remove.append(num2)

            keys_viewed.append(num1)

    voronois_copy = voronoi_points
    for item in to_remove:
        if key[item] in voronois_copy:
            voronois_copy.remove(key[item])
    return voronois_copy


def step2(voronoi_points, structure, radius, anions, voronoi_point, min_neighbors):
    """
    Remove sites with Voronoi points that sit on the bonds between anions

    :param voronoi_points: (List of Sites) list of sites with candidate Voronoi points
    :param structure: (Structure) target structure
    :param radius: (float) distance to which we check for anions
    :param anions: (List of Strings) list of species which we consider anions

    """
    big_cations = ['K', 'Ba', 'Sr']
    to_remove = []
    counter = 0
    for point in voronoi_points:
        counter += 1
        #print counter
        neighbors = structure.get_neighbors(point, radius)
        for neighbor in neighbors:
            if not neighbor[0].species_string in anions and not neighbor[0].species_string == voronoi_point \
                    and neighbor[1] < radius:
                cv_distance = neighbor[1]
                co_distances = []

                if neighbor[0].species_string in big_cations:
                    for neighbor2 in structure.get_neighbors(neighbor[0], 3.5):
                        if neighbor2[0].species_string in anions and neighbor2[1] < radius:
                            co_distances.append(neighbor2[1])
                    for co_distance in co_distances:
                        # I can see a small area within which this might not work, but I don't think that area happens
                        # in actual materials
                        if cv_distance <= co_distance:
                            to_remove.append(point)
                else:
                    for neighbor2 in structure.get_neighbors(neighbor[0], radius):
                        if neighbor2[0].species_string in anions and neighbor2[1] < radius:
                            co_distances.append(neighbor2[1])
                    for co_distance in co_distances:
                        # I can see a small area within which this might not work, but I don't think that area happens
                        # in actual materials
                        if cv_distance <= co_distance*(1.1*(2**0.5)*0.5):
                            to_remove.append(point)

    for point in voronoi_points:
        neighbors = get_anion_neighbors(point, structure, radius, anions)
        if len(neighbors) < min_neighbors:
            to_remove.append(point)
    voronois_copy = voronoi_points
    # go through similarVoronoi to reduce to only one voronoi with a given set of peripheral anions
    for item in to_remove:
        if item in voronois_copy:
            voronois_copy.remove(item)
    return voronois_copy

if __name__ == '__main__':
    """
    voronoipoint = 'Rn'
    radius = 3.0
    minneighbors = 4
    anions = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S', 'N', 'N3-']
    numMatches = 4

    print "accessing structure with voronois points"
    task_ids = []
    toexclude = []
    allstructures = {}
    taskidsList = open('structures_with_voronoi/taskid.txt', 'r')
    for line in taskidsList:
        splits = line.split("\t")
        for task_id in splits:
            task_ids.append(int(task_id.strip().replace(".cif", "")))
    print task_ids
    print len(task_ids)
    print "Loading structure files"
    counter = 1
    for id in task_ids:
        print "Loading file", counter
        allstructures[id] = Structure.from_file("structures_with_voronoi/%s.cif"%(id))
        counter += 1
    structuresWithVoronoi = {}
    counter = 1
    for id in allstructures.keys():
        print counter, ":", id
        structure = allstructures[id]

        # find all voronoi sites in the structure
        print "splitting sites into voronoi and non-voronoi"
        allvoronois, othersites = splitSites(structure)

        # process potential empty polyhedral center sites
        print "processing structure, step 1"
        voronois = step1(allvoronois, structure)

        # create intermediate structure2
        print "creating structure2"
        structure2 = createStructure(voronois, othersites, structure)
        structure2.to(filename="structures_completed/%s_2.cif" % id)

        if len(voronois) <= 1000:
            #counter += 1
            continue

        print "processing structure, step 2 (may take a few minutes)"
        # find polyhedral of voronoi, organizing them by pairing those with the same set of peripheral anions
        # match each voronoi point in voronois to a number - dictionary
        # create a method like getCationPolyhedral that finds the polyhedral for the specific site
        # check polyhedral in the same way
        # use the number from 2nd step to choose which sites to remove
        voronois = step3(voronois, structure2, radius, anions)

        print "creating structure3"
        structure3 = createStructure(voronois, othersites, structure2)
        structure3.to(filename="structures_completed/%s_3.cif"%(id))

        print "processing structure, step 3"
        # add criterion based on C-V and C-O distances
        # if C-O < 2.5, if C-V distance > C-O distance, then remove
        # for cation closest to V with all nearby oxygen
        voronois = step2(voronois, structure3, radius, anions)

        print "creating structure4"
        structure4 = createStructure(voronois, othersites, structure3)
        structure4.to(filename="structures_completed/%s.cif"%(id))

        structuresWithVoronoi[id] = structure4
        counter += 1

    """

"""
structure1 = Structure.from_file("LiCoO2.cif")

structure2 = add_voronoi(structure1)

voronoi_species = 'Rn'
max_radius = 2.5
min_num_neighbors = 4
anions_list = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S', 'N', 'N3-']
num_matches = 4

# find all voronoi sites in the structure
print "splitting sites into voronoi and non-voronoi"
allvoronois, othersites = split_sites(structure2, voronoi_species)

# process potential empty polyhedral center sites
print "processing structure, step 1"
voronois = step1(allvoronois, structure2, max_radius, anions_list, min_num_neighbors)


# create intermediate structure2
print "creating structure2"
structure3 = structure_from_sites(voronois, othersites, structure2)


print "processing structure, step 3"
# add criterion based on C-V and C-O distances
# if C-O < 2.5, if C-V distance < C-O distance, then remove
# for cation closest to V with all nearby oxygen
voronois = step3(voronois, structure3, max_radius, anions_list, voronoi_species)


print "creating structure3"
structure4 = structure_from_sites(voronois, othersites, structure3)

print "processing structure, step 2 (may take a few minutes)"
# find polyhedral of voronoi, organizing them by pairing those with the same set of peripheral anions
# match each voronoi point in voronois to a number - dictionary
# create a method like getCationPolyhedral that finds the polyhedral for the specific site
# check polyhedral in the same way
# use the number from 2nd step to choose which sites to remove
voronois = step2(voronois, structure4, max_radius, anions_list, num_matches)

print "creating structure4"
structure5 = structure_from_sites(voronois, othersites, structure4)

voronois = step3(voronois, structure5, max_radius, anions_list, voronoi_species)

structure6 = structure_from_sites(voronois, othersites, structure5)

structure6.to(filename="LiCoO2_test.cif")
"""
