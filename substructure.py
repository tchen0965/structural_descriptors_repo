#!/usr/bin/env python

"""
This module provides classes to define substructural components of periodic structures
"""

from __future__ import division

__author__ = "Stephen_Dacek"
__credits__ = 'Lusann Wren Yang'
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Stephen Dacek"
__email__ = "sdacek@mit.edu"

from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder, VoronoiConnectivity
from pymatgen.serializers.json_coders import MSONable
from substructure_specie import SubStructureSpecie, SubStructureSite
from pymatgen.core.ion import Ion
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.transformations.site_transformations import TranslateSitesTransformation
from collections import Counter
from copy import deepcopy
import numpy as np


class SubStructure(MSONable):
    """
    A collection of SubStructureSpecies defining a substructure.
    """

    def __init__(self, peripheral_sub_species, central_sub_species=None, round_weight_to=0.05):

        """
        Create a Substructure Object

        Args:
            peripheral_sub_species (iterable): iterable of SubstructureSpecie objects
                coordinating the central point in the substructure

            central_sub_species (SubstructureSpecie or None): specie at the center of
                the substructure if it exists

            round_weight_to (float): round the weight of all substructural
                components to the specified precision
        """

        if all(isinstance(subspecies, SubStructureSpecie) for subspecies in peripheral_sub_species):
            self.peripheral_subspecies = list(peripheral_sub_species)
        else:
            raise ValueError('Every object in iterable peripheral_sub_species'
                             ' must be an instance of SubStructureSpecie ')

        if central_sub_species and isinstance(central_sub_species, SubStructureSpecie):
            self.central_subspecies = central_sub_species

        else:
            raise ValueError('central_sub_species must be an instance of SubStructureSpecie'
                             ' or None ')

        if round_weight_to:
            self._round_weight(round_weight_to)

        self.round_weight_to = round_weight_to

    def self_product(self):

        return 1.0 + self.weight_sum() if self.central_subspecies \
            else self.weight_sum()

    def _round_weight(self, precision):

        """
        Round the weights of all component species
        to the specified precision
        """

        self.central_subspecies.weight = round(self.central_subspecies.weight / precision) * precision

        for peripheralSpecies in self.peripheral_subspecies:
            peripheralSpecies.weight = round(peripheralSpecies.weight / precision) * precision

    def __len__(self):
        return 1 + len(self.peripheral_subspecies) if self.central_subspecies \
            else len(self.peripheral_subspecies)

    def __eq__(self, other):
        if isinstance(other, SubStructure):
            return (self.central_subspecies == other.central_subspecies) and \
                   Counter(self.peripheral_subspecies) == Counter(other.peripheral_species)
        else:
            return False

    def __str__(self):
        return 'Composition: {}, Product: {}'.format(self.to_Ion(),
                                                     self.self_product())

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash(self.__repr__())

    @classmethod
    def from_dict(cls, d):

        p_subspecies = [SubStructureSpecie.from_dict(x) for x in d['peripheral_subspecies']]

        central_subspecies = SubStructureSpecie.from_dict(d.get('central_subspecies')) if d.get('central_subspecies', None) else None

        #return SubStructure(p_subspecies, central_subspecies,
        #                    d.get('round_weight_to', None))

        return SubStructure(p_subspecies, central_subspecies,
                            0.05)

    @property
    def as_dict(self):
        """""
        Json-serialization dict representation of the substructure.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "central_subspecies": self.central_subspecies.as_dict,
                "peripheral_subspecies": [i.to_dict for i in self.peripheral_subspecies],
                "round_weight_to": self.round_weight_to}

    def weight_sum(self):
        """
        Sum of all peripheral ion weights
        """
        return sum([subspecies.weight for subspecies in self.peripheral_subspecies])

    def to_Ion(self, weight=True):
        composition = self.to_composition(weight)
        charge = sum([k.oxi_state * v for k, v in composition.iteritems()])
        return Ion(composition, charge)

    def to_composition(self, weight=True):
        comp_dict = {}
        for subspecies in self.peripheral_subspecies + [self.central_subspecies]:
            if comp_dict.get(subspecies.specie):
                comp_dict[subspecies.specie] += subspecies.weight if weight else 1
            else:
                comp_dict[subspecies.specie] = subspecies.weight if weight else 1
        return Composition(comp_dict)


class ExtendedSubStructure(MSONable):
    """
    A collection of SubStructureSites defining an ExtendedSubStructure.
    """

    def __init__(self, peripheral_sub_sites, central_sub_site=None, round_weight_to=0.05):

        """
        Create an ExtendedSubStructure Object

        Args:
            peripheral_sub_sites (iterable): iterable of SubstructureSite objects
                coordinating the central point in the substructure

            central_sub_site (SubstructureSite or None): specie at the center of
                the substructure if it exists

            round_weight_to (float): round the weight of all substructural
                components to the specified precision

            translate_coords (bool): Translate central ion to 0, 0, 0. If central Ion
             is None, translates the average coordinate of the substructure to 0,0,0

        """

        if all(isinstance(subsite, SubStructureSite) for subsite in peripheral_sub_sites):
            self.peripheral_subsites = list(peripheral_sub_sites)
        else:
            raise ValueError('Every object in iterable peripheral_sub_sites'
                             ' must be an instance of SubStructureSite')

        if central_sub_site and isinstance(central_sub_site, SubStructureSite):
            self.central_subsite = central_sub_site
        else:
            raise ValueError('central_sub_site must be an instance of SubStructureSite'
                             ' or None ')

        if round_weight_to:
            self._round_weight(round_weight_to)

        self.round_weight_to = round_weight_to

    def self_product(self):

        return 1.0 + self.weight_sum() if self.central_subsite \
            else self.weight_sum()

    def _round_weight(self, precision):

        """
        Round the weights of all component species
        to the specified precision
        """

        self.central_subsite.weight = round(self.central_subsite.weight / precision) * precision

        for peripheralSubSite in self.peripheral_subsites:
            peripheralSubSite.weight = round(peripheralSubSite.weight / precision) * precision

    def __len__(self):
        return 1 + len(self.peripheral_subsites) if self.central_subsite \
            else len(self.peripheral_subsites)

    def __eq__(self, other):
        if isinstance(other, SubStructure):
            return (self.central_subsite == other.central_subsite) and \
                   Counter(self.peripheral_subsites) == Counter(other.peripheral_subsites)
        else:
            return False

    def __str__(self):
        return 'Composition: {}, Product: {}'.format(self.to_Ion(),
                                                     self.self_product())

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash(self.__repr__())

    @classmethod
    def from_dict(cls, d):

        p_subsites = [SubStructureSite.from_dict(x) for x in d['peripheral_subsites']]

        central_subsite = SubStructureSite.from_dict(d.get('central_subsite')) if d.get('central_subsite', None) else None

        return SubStructure(p_subsites, central_subsite,
                            0.05)

    def translate_sites(self, lattice):
        s = Structure(lattice, [self.central_subsite.specie] + [site.specie for site in self.peripheral_subsites], [self.central_subsite.site.coords] + [site.site.coords for site in self.peripheral_subsites], coords_are_cartesian=True)
        trans = TranslateSitesTransformation(range(len(s)), -(lattice.get_fractional_coords(s[0].coords)))
        new_s = trans.apply_transformation(s)
        trans2 = TranslateSitesTransformation(range(len(s)), (-0.5, -0.5, -0.5))
        new_s = trans2.apply_transformation(Structure.from_sites(new_s.sites, to_unit_cell=True))
        self.peripheral_subsites = [SubStructureSite.from_coords_and_specie(site.coords, site.specie) for site in new_s.sites[1:]]
        self.central_subsite = SubStructureSite.from_coords_and_specie(new_s[0].coords, new_s[0].specie)

    def to(self, fmt, filename, weight=True, lattice=None):
        composition = self.to_composition(weight)
        charge = sum([k.oxi_state * v for k, v in composition.iteritems()])
        s = Structure(lattice,
                      [self.central_subsite.specie] + [site.specie for site in self.peripheral_subsites],
                      [self.central_subsite.site.coords] + [site.site.coords for site in self.peripheral_subsites], coords_are_cartesian=True)
        s.to(fmt, filename)


    @property
    def as_dict(self):
        """""
        Json-serialization dict representation of the substructure.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "central_subsite": self.central_subsite.as_dict,
                "peripheral_subsites": [i.to_dict for i in self.peripheral_subsites],
                "round_weight_to": self.round_weight_to}

    def weight_sum(self):
        """
        Sum of all peripheral ion weights
        """
        return sum([subsite.weight for subsite in self.peripheral_subsites])

    def to_Ion(self, weight=True):
        composition = self.to_composition(weight)
        charge = sum([k.oxi_state * v for k, v in composition.iteritems()])
        return Ion(composition, charge)

    def to_composition(self, weight=True):
        comp_dict = {}
        for subsite in self.peripheral_subsites + [self.central_subsite]:
            if comp_dict.get(subsite.specie):
                comp_dict[subsite.specie] += subsite.weight if weight else 1
            else:
                comp_dict[subsite.specie] = subsite.weight if weight else 1
        return Composition(comp_dict)


def find_cutoff_weight(weights):
    cut_list = [i + 1 for i, dw in enumerate(np.abs(np.diff(weights))) if dw > 0.5]

    return cut_list[0] if cut_list else None


def substructures_from_structure(structure, weight_cutoff=1e-2):
    """
    Helper method to calculate substructural components from a
    pymatgen structure object.

    Args:
        structure (Structure): Input structure
        weight_cutoff (float): minimum solid angle weight for
                               inclusion in a substructure

    Returns:
            A list of Substructure objects for the given structure
    """

    vcf = VoronoiCoordFinder(structure)
    substructures = []

    for i, site in enumerate(structure):

        charge = sum([getattr(specie, "oxi_state", 0) * amt
                      for specie, amt in site.species_and_occu.items()])

        central_subspecies = SubStructureSpecie(site.specie.symbol,
                                         oxidation_state=charge,
                                         properties=site.properties,
                                         weight=1.0)

        peripheral_subspecies = []
        for peripheral_site, weight in vcf.get_voronoi_polyhedra(i).iteritems():

            if weight > weight_cutoff:
                charge = sum([getattr(specie, "oxi_state", 0) * amt
                              for specie, amt in
                              peripheral_site.species_and_occu.items()])

                peripheral_subspecies.append(
                    SubStructureSpecie(peripheral_site.specie.symbol,
                                       oxidation_state=charge,
                                       properties=peripheral_site.properties,
                                       weight=weight))

        substructures.append(SubStructure(peripheral_subspecies, central_sub_species=central_subspecies))

    return substructures


def substructures_from_structure2(structure, weight_cutoff=0.01, use_cutoff=False, include_coords=True):
    """
    Helper method to calculate substructural components from a
    pymatgen structure object, including coords.

    Args:
        structure (Structure): Input structure
        weight_cutoff (float): minimum solid angle weight for
                               inclusion in a substructure

    Returns:
            A list of Substructure objects for the given structure
    """

    substructures = []
    vcf = VoronoiConnectivity(structure)
    connectivity_array = vcf.connectivity_array
    valid_sites = np.where(connectivity_array > 0)

    central_ions = np.empty((structure.num_sites), dtype=SubStructureSpecie)
    #create list of sites in structure
    for i, site in enumerate(structure):
        charge = sum([getattr(specie, "oxi_state", 0) * amt
                      for specie, amt in site.species_and_occu.items()])
        if include_coords:
            central_ion = SubStructureSite(site.coords, site.specie.symbol,
                                           oxidation_state=charge,
                                           properties=site.properties,
                                           weight=1.0)
        else:
            central_ion = SubStructureSpecie(site.specie.symbol,
                                             oxidation_state=charge,
                                             properties=site.properties,
                                             weight=1.0)
        central_ions[i] = central_ion

    for i, site in enumerate(structure):

        indicies = np.where((valid_sites[0] == i))[0]

        weights = connectivity_array[valid_sites][indicies] / np.max(connectivity_array[valid_sites][indicies])

        sort_inds = np.argsort(weights)[::-1]

        meets_weight_criteria = np.where(weights[sort_inds] > weight_cutoff)[0]

        if use_cutoff:
            cutoff = find_cutoff_weight(weights[sort_inds][meets_weight_criteria])
        else:
            cutoff = None
        coordinating_sites = valid_sites[1][indicies][sort_inds][meets_weight_criteria][: cutoff]

        substructure_sites = deepcopy(central_ions[coordinating_sites])

        [setattr(sites, 'weight', weight) for weight, sites in
         zip(weights[sort_inds][meets_weight_criteria][:cutoff], substructure_sites)]

        if include_coords:
            substructures.append(ExtendedSubStructure(substructure_sites, central_sub_site=central_ions[i]))
        else:
            substructures.append(SubStructure(substructure_sites, central_sub_species=central_ions[i]))

    return substructures


if __name__ == '__main__':
    s = Structure.from_file('LiCoO2.cif', True, True)

    for substruct in substructures_from_structure(s):
        print substruct.central_subspecies, substruct.weight_sum()
    for substruct in substructures_from_structure2(s):
        print substruct.central_subsite, substruct.weight_sum()