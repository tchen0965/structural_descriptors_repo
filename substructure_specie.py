#!/usr/bin/env python

"""
This module provides classes to define species components of a substructure
"""

from __future__ import division

__author__ = "Stephen_Dacek"
__credits__ = 'Lusann Wren Yang'
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Stephen Dacek"
__email__ = "sdacek@mit.edu"

from pymatgen.serializers.json_coders import MSONable
from pymatgen.core.periodic_table import Specie
from pymatgen.core.sites import Site
import numpy as np


class SubStructureSpecie(MSONable):

    """
    An extension of species to incorporate solid angle weights
    No-coordinates are necesarry
    """

    def __init__(self, symbol, oxidation_state=0.0, weight=1.0, properties=None):

        """
        Create a Substructure Object

        Args:
            symbol (str): Elemental symbol

            oxidation_state (float): Oxidation state of the species

            weight (float): Weight of the species in a substructure

            properties (dict): Dictionary of supported properties for species

        """

        self.specie = Specie(symbol, oxidation_state, properties)
        self.weight = weight

    @staticmethod
    def from_specie(specie, weight=0.0):
        """
        Helper method to create SubStructureSpecie from a species and weight
        """
        return SubStructureSpecie(symbol=specie.symbol, oxidation_state=specie.oxi_state, properties=specie._properties, weight=weight)

    @classmethod
    def from_dict(cls, d):
        return SubStructureSpecie.from_specie(Specie.from_dict(d['specie']), d.get('weight', 0.0))

    @property
    def as_dict(self):
        """
        JSON serialization of a SubStructureSpecie
        """
        sp = Specie(self.specie.symbol, self.specie.oxi_state, self.specie._properties)
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "specie": sp.as_dict(),
                "weight": self.weight}

    def __str__(self):
        return 'Specie: {}, Weight: {}'.format(self.specie, self.weight)

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, other, tol=1e-2):
        if isinstance(other, SubStructureSpecie):
            return self.specie == other.specie and np.abs(self.weight - other.weight) < tol
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)


class SubStructureSite(MSONable):

    """
    An extension of species to incorporate solid angle weights
    This version extends substructures to incorporate coordinates
    """

    def __init__(self, coordinates, symbol, oxidation_state=0.0, weight=1.0, properties=None):

        """
        Create a Substructure Object

        Args:
            coordinates (iter): Coordinates of substructure site

            symbol (str): Elemental symbol

            oxidation_state (float): Oxidation state of the species

            weight (float): Weight of the species in a substructure

            properties (dict): Dictionary of supported properties for species

        """

        self.specie = Specie(symbol, oxidation_state, properties)
        self.site = Site(self.specie, coordinates)
        self.weight = weight

    @staticmethod
    def from_coords_and_specie(coordinates, specie, weight=0.0):
        """
        Helper method to create SubStructureSite from coordinates species and weight
        """
        return SubStructureSite(coordinates, symbol=specie.symbol, oxidation_state=specie.oxi_state,
                                properties=specie._properties, weight=weight)


    @staticmethod
    def from_site(site, weight=0.0):
        """
        Helper method to create SubStructureSite from coordinates species and weight
        """
        return SubStructureSite.from_coords_and_specie(site.coords, site.specie, weight=weight)

    @classmethod
    def from_dict(cls, d):
        return SubStructureSite.from_coords_and_specie(d['coords'],
                                                       Specie.from_dict(d['specie']), d.get('weight', 0.0))

    @property
    def as_dict(self):
        """
        JSON serialization of a SubStructureSpecie
        """
        sp = Specie(self.specie.symbol, self.specie.oxi_state, self.specie._properties)
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "coords": self.site.coords,
                "specie": sp.as_dict(),
                "weight": self.weight}

    def __str__(self):
        return 'Coords: {}, Specie: {}, Weight: {}'.format(self.site.coords, self.specie, self.weight)

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, other, tol=1e-2):
        if isinstance(other, SubStructureSite):
            return self.specie == other.specie and np.abs(self.weight - other.weight) < tol
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)