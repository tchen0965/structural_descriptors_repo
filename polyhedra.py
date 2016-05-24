__author__ = 'Tina_Chen'


class Polyhedra(object):
    """
    Object representing a polyhedra in a structure with central ion site "cation" and surrounding ions site
    "peripheralIons"

    """

    def __init__(self, cation, peripheral_ions):
        """

        :param cation: (Site) Site object representing the central cation
        :param peripheral_ions: (list of Sites) list of Site objects representing the central cation's surrounding ions
        :param isite: (integer) index of site in original structure (not supercell)
        """

        self.central_ion = cation
        self.central_ion_name = cation.species_string
        self.peripheral_ions = peripheral_ions
        # Number to identify unique cation site in center cell from other sites with same species
        self.cation_num = 1

    def set_site_number(self, cation_number):
        """
        :param cation_number: number differentiating sites with the same species sitting on the site
        """
        self.cation_num = cation_number

    def get_num_connections(self, other):
        """
        Gives the connectivity between the given polyhedra and another polyhedra by counting the number of atoms shared
        between the two polyhedra

        :param other: (Polyhedra) target Polyhedra against which we are checking for connectivity with the current
        polyhedra
        :return: (int) Integer giving the number of peripheral ions shared by the current polyhedra and the target
        polyhedra; returns -1 if the same polyhedra is being compared to itself
        """
        return len(self.get_connections(other))

    def get_connections(self, other):
        """
        Gives the shared sites between the current Polyhedra and the given Polyhedra

        :param other: (Polyhedra) target Polyhedra against which we are checking for connectivity with the current
        polyhedra
        :return: (list) list of sites shared between the current Polyhedra and the given Polyhedra
        """

        if self == other:
            raise ValueError('Checking connections between exact same polyhedra')

        shared_sites = []

        for site in other.peripheral_ions:
            for c_site in self.peripheral_ions:
                if c_site == site:
                    shared_sites.append(site)
        return shared_sites

    def __eq__(self, other):
        """

        :param other: another Polyhedra
        :return: true if the central ion site of the current polyhedra is the same as the other polyhedra
        """
        if self.central_ion == other.central_ion:
            return True
        else:
            return False

    def __str__(self):
        """
        List of properties of the polyhedra: central species name, central species site, list of peripheral ions of
        the polyhedra
        """

        return str([self.central_ion_name,
                    [self.central_ion.frac_coords[0], self.central_ion.frac_coords[1], self.central_ion.frac_coords[2]],
                    self.peripheral_ions])
