__author__ = 'Tina'


from pymatgen.core.structure import Structure
from pymatgen import PeriodicSite
from effective_coordination_number_from_structure import calculateBondWeight, getAllCationPolyhedral
import math
import numpy


def getConnection(polyhedral1, polyhedral2):
    if polyhedral1[0] == polyhedral2[0]:
        #print "Same central cation"
        return -1
    peripheral1 = polyhedral1[1:len(polyhedral1)]
    peripheral2 = polyhedral2[1:len(polyhedral2)]

    connections = 0
    for site1 in peripheral1:
        for site2 in peripheral2:
            if site1[0] == site2[0]: #check whether the sites are equal
                connections += 1
                break

    return connections

def calculateConnectivity(structure):

    polyhedraList = getAllCationPolyhedralWReflections(structure)
    print len(polyhedraList)

    connections = {}
    cationList = []
    for polyhedra1 in polyhedraList:
        cation = polyhedra1[0][0].species_string
        connected = False
        if not cation in cationList:
            cationList.append(cation)
            connections[cation] = [0,0,0,0]
        for polyhedra2 in polyhedraList:
            connection = getConnection(polyhedra1[0], polyhedra2[0])
            if connection < 0:
                continue
            if connection == 1:
                connections[cation][1] += 1
                connected = True
            if connection == 2:
                connections[cation][2] += 1
                connected = True
            if connection >= 3:
                connections[cation][3] += 1
                connected = True
        if connected == False:
            connections[cation][0] += 1


    return connections

def findReflections(polyhedra):
    return null

def getAllCationPolyhedralWReflections(structure, margin = 0.05, radius = 3.0):

    anions = ['O2-', 'O', 'F-', 'F', 'Cl-', 'Cl', 'I-', 'I', 'Br-', 'Br', 'S2-', 'S']

    # pick out all sites with cations (not just "cation") as the species at the site
    cationSites = []
    for site in structure.sites:
        if not site.species_string in anions:
            print site.frac_coords
            cationSites.append(site)
    reflections = []
    for site in cationSites:
        reflections = reflections + checkReflections(site, margin)
    useTemp = False
    if len(reflections) > 0:
        for site in reflections:
            cationSites.append(site)
        useTemp = True
        tempspecies = []
        tempfraccoords = []
        for site in structure.sites:
            tempspecies.append(site.species_string)
            tempfraccoords.append(site.frac_coords)
        for reflection in reflections:
            tempspecies.append(reflection.species_string)
            tempfraccoords.append(reflection.frac_coords)
        tempstruct = Structure(structure.lattice, tempspecies, tempfraccoords)
    if useTemp:
        structure = tempstruct
    print structure.formula
    polyhedralList = [] #list of polyhedrals with cations at the center
    for site in cationSites:
        sites = [site] #list of sites that will be in polyhedral
        anionSites = [] # list of all neighboring anions
        bondlengths = [] #list of bondslengths of anions
        bondweights = []
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

def checkReflections(site, margin):
    x = 0
    y = 0
    z = 0
    reflections = []
    toReflect = site
    if toReflect.frac_coords[0] >= 0.0 and toReflect.frac_coords[0] <= margin:
        x = 1
    if toReflect.frac_coords[0] <= 1.0 and toReflect.frac_coords[0] >= 1.0 - margin:
        x = -1
    if toReflect.frac_coords[1] >= 0.0 and toReflect.frac_coords[1] <= margin:
        y = 1
    if toReflect.frac_coords[1] <= 1.0 and toReflect.frac_coords[1] >= 1.0 - margin:
        y = -1
    if toReflect.frac_coords[2] >= 0.0 and toReflect.frac_coords[2] <= margin:
        z = 1
    if toReflect.frac_coords[2] <= 1.0 and toReflect.frac_coords[2] >= 1.0 - margin:
        z = -1

    if x != 0:
        reflections.append(PeriodicSite(site.species_string, site.frac_coords + [x*1.0, 0, 0],
                                            site.lattice))
        if y != 0:
            reflections.append(PeriodicSite(site.species_string, site.frac_coords + [x*1.0, y*1.0, 0],
                                                site.lattice))
            if z != 0:
                reflections.append(PeriodicSite(site.species_string, site.frac_coords + [x*1.0, y*1.0, z*1.0],
                                                site.lattice))
        if z != 0:
            reflections.append(PeriodicSite(site.species_string, site.frac_coords + [x*1.0, 0, z*1.0],
                                            site.lattice))

    if y != 0:
        reflections.append(PeriodicSite(site.species_string, site.frac_coords + [0, y*1.0, 0],
                                            site.lattice))
        if z != 0:
            reflections.append(PeriodicSite(site.species_string, site.frac_coords + [0, y*1.0, z*1.0],
                                                site.lattice))

    if z != 0:
        reflections.append(PeriodicSite(site.species_string, site.frac_coords + [0, 0, z*1.0],
                                            site.lattice))

    return reflections

if __name__ == '__main__':
    s = Structure.from_file('LiCoO2.cif', True, False)
    print s
    print calculateConnectivity(s)