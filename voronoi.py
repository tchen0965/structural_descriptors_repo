from pyhull.voronoi import VoronoiTess
from pyhull.delaunay import DelaunayTri
from pymatgen.transformations.transformation_abc import AbstractTransformation
from pymatgen.core.sites import PeriodicSite
from pymatgen.util.coord_utils import barycentric_coords
from pymatgen.analysis.structure_analyzer import solid_angle
from pymatgen.core.periodic_table import get_el_sp, DummySpecie
from pymatgen.core.structure import IStructure


import itertools
import logging
import numpy as np


class VoronoiInsertionTransformation(AbstractTransformation):
    """
    This transformation inserts sites at the voronoi points of
    the structure. Can also insert at the center of the voronoi
    faces. It can also enforce bond distances (by scaling these
    points back towards other atoms.
    """

    def __init__(self, sp, num = None, bond_lengths = None, 
                 interior_only = False,
                 midpoints = True,
                 min_solid_angle = 0,
                 site_distance_tol = 0.1):
        """
        Args:
            sp:
                specie to add
            num:
                number to add (will end up with a disordered structure)
            bond_lengths:
                dictionary of bond lengths. For example, if {'O' : 1} is
                used, all points will be scaled from their voronoi locations
                to be exactly 1 angstrom from the oxygen (if the oxygen is
                at the center of the voronoi region). In the case of a point
                neighboring 2 oxygen atoms, 2 points will be added - one for
                each oxygen
            interior_only:
                only uses voronoi points that are inside the tetrahedra
                formed by their 4 nearest neighbors. This is useful because
                when this is not true, there is another point accessible
                without going through a bottleneck that is strictly larger.
                Currently this doesn't work well with bond lengths because 
                finding the correct association between sites and voronoi 
                points is difficult. Doesn't affect the delaunay midpoints
            midpoints:
                whether to add the points at the center of the voronoi
                faces (midpoints of the delaunay triangulation)
            min_solid_angle:
                Threshold for creating points at the delaunay midpoints. The
                solid angle is a metric for how much two sites are neighbors of
                each other (the solid angle swept out by the voronoi face 
                corresponding to the two sites). Typical values are 0 to Pi, with
                Pi corresponding to a tetrahedrally coordinated site.
            site_distance_tol:
                Tolerance passed to to_unit_cell(). When there are multiple
                sites within this tolerance, only one is used
        """
        self._sp = get_el_sp(sp)
        self._n = num
        self._interior = interior_only
        self._midpoints = midpoints
        self._msa = min_solid_angle
        self._sdt = site_distance_tol
        
        #transform keys to species
        if bond_lengths:
            self._bl = dict(zip(map(lambda x: get_el_sp(x)
                                    , bond_lengths.keys())
                                , bond_lengths.values()))
        else:
            self._bl = {}

    def apply_transformation(self, structure):
        #get the new sites
        new_sites = self._voronoi_points(structure)
        if self._midpoints:
            dms = self._delaunay_midpoints(structure)
            for i in range(len(new_sites)):
                new_sites[i].extend(dms[i])
        
        #insert the new sites with dummy species
        se = structure
        dummy = DummySpecie()
        for i, sites_list in enumerate(new_sites):
            for ns in sites_list:
                if self._bl == {}:
                    se.append(self._sp, ns, coords_are_cartesian = True
                                   , validate_proximity = False)
                elif structure[i].specie in self._bl.keys():
                    bv = ns - structure[i].coords
                    bv *= self._bl[structure[i].specie] / np.linalg.norm(bv)
                    coords = structure[i].coords + bv
                    se.append_site(self._sp, coords, coords_are_cartesian = True
                                   , validate_proximity = False)
        logging.debug('Starting moving to unit cell')
        se=IStructure.get_reduced_structure(se)
        logging.debug('Finished moving to unit cell')
        
        structure = se#.modified_structure
        
        #delete sites that are closer than the bond distance
        to_delete = []
        for sp, d in self._bl.items():
            neighbors = structure.get_all_neighbors(d*0.999)
            for i in range(len(structure)):
                if structure[i].specie == dummy:
                    for n in neighbors[i]:
                        if n[0].specie == sp:
                            to_delete.append(i)
                            break
        se.remove_sites(to_delete)
        
        #replace the dummy species to get the correct composition
        if self._n:
            #amt = self._n / se.modified_structure.composition[dummy]
            amt = self._n / se.composition[dummy]
            se.replace_species({dummy : {self._sp: amt}})
        
        return se#.modified_structure
                    
    def _fake_periodic(self, structure):
        """
        Makes a 5x5x5 supercell to fake a periodic boundary for pyhull
        """
        coords = np.array(structure.cart_coords)
        lattice = structure.lattice.matrix
        all_coords = []
        for x in itertools.product([0,1,-1,2,-2], repeat = 3):
            offset = np.dot(x, lattice)
            offset_points = coords + offset
            all_coords.extend(offset_points.tolist())
        return all_coords
    
    def _delaunay_midpoints(self, structure):
        """
        returns the midpoints of the delaunay triangulation for a structure. 
        returns in the format of a list of lists, where the first index 
        corresponds to the index of the site in the original structure that 
        neighbors the point. I.e. dms[0] is the list of delaunay midpoints 
        around the site structure[0]
        
        We can't simply use the ridges from the voronoi tesselation because
        the ridges don't include nearest neighbors with a solid angle of 0
        (whose voronoi volumes are edge sharing)
        """
        logging.debug('Starting calculating Delaunay points')
        all_coords = self._fake_periodic(structure)
        vt = VoronoiTess(all_coords)
        dt = DelaunayTri(all_coords)
        
        v_array = np.array(vt.vertices)
        p_array = np.array(vt.points)
        dms = [ [] for i in range(len(structure))]
        
        for simplex in dt.vertices:
            if min(simplex) >= len(structure):
                continue
            for pt1, pt2 in itertools.permutations(simplex, 2):
                if pt1 >= len(structure):
                    continue
                center = structure[pt1].coords
                
                #get points that make up the shared plane
                if pt2>pt1:
                    vertices = vt.ridges.get((pt1, pt2), None)
                else:
                    vertices = vt.ridges.get((pt2, pt1), None)
                
                #get the solid angle, if pt1 and pt2 share a plane
                if vertices is not None:
                    sa = solid_angle(center, v_array[vertices])
                else:
                    sa = 0
                    
                new_coord = (p_array[pt1] + p_array[pt2]) / 2
                if sa >= self._msa:
                    dms[pt1].append(new_coord)
                
        logging.debug('Finished calculating Delaunay points')
        
        return dms
        
    
    def _voronoi_points(self, structure):
        """
        returns the voronoi points for a structure. returns in the 
        format of a list of lists, where the first index corresponds
        to the index of the site in the original structure that neighbors
        the point. I.e. simplices[0] is the list of points around the site
        structure[0]
        """
        all_coords = self._fake_periodic(structure)
        logging.debug('Starting calculating Voronoi points')
        vt = VoronoiTess(all_coords)
        simplices = []
        #only return the coordinates associated with the
        #points in the original structure
        for i, r in enumerate(vt.regions[:len(structure)]):
            points = []
            for v in r:
                if self._interior:
                    new_site = PeriodicSite('X', vt.vertices[v], structure.lattice,
                                            coords_are_cartesian = True)
                    d = structure[i].distance(new_site)
                    s_points = structure.get_neighbors_in_shell(new_site.coords, 
                                                              d, 0.01)
                    s_points = map(lambda x: x[0].coords, s_points)
                    #brute force check that point is in its polyhedron 
                    #(len(s_points) should usually be 4, except in high 
                    #symmetry cases)
                    for simplex in itertools.combinations(s_points, 4):
                        try:
                            bcs = barycentric_coords(vt.vertices[v], 
                                                     np.array(simplex))
                            if np.min(bcs) > 0:
                                points.append(vt.vertices[v])
                                break
                        except np.linalg.LinAlgError as err:
                            if 'Singular matrix' not in err.message:
                                raise
                            
                else:
                    points.append(vt.vertices[v])
            simplices.append(points)
        logging.debug('Finished calculating Voronoi points')
        return simplices

    def __str__(self):
        return "Voronoi Insertion Transformation : " + \
            "Species to insert = {}".format(str(self._sp))

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return False
    
    @property
    def use_multiprocessing(self):
        return True

    @staticmethod
    def as_dict(self):
        bl = dict(zip(map(lambda x: str(x)
                          , self._bl.keys())
                          , self._bl.values()))
        return {"name": self.__class__.__name__, "version": None,
                "init_args": {"sp" : str(self._sp),
                              "num" : self._n,
                              "bond_lengths" : bl,
                              "midpoints" : self._midpoints},
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}
