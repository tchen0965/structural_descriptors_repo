from pymatgen.analysis.chemenv.coordination_environments import coordination_geometry_finder as polyfinder
from pymatgen.analysis.chemenv.coordination_environments import chemenv_strategies as strategies
from pymatgen.analysis.chemenv.coordination_environments import structure_environments as se
from pymatgen.util.coord_utils import find_in_coord_list_pbc


__author__ = 'Tina'

"""
Using the coordination_geometry_finder from pymatgen's chemenv to find
"""


def find_site_ce(structure, target_site):
    """
    Returns mp_symbol of site environment of given site
    :param structure: (Structure) Structure containing target Site
    :param target_site: (Site) target Site
    :return: (String) mp_symbol of site coordination environment of given site
    """

    # Find index of target_site
    #target_isite = find_in_coord_list_pbc(structure.frac_coords, target_site.frac_coords)[0]

    # Set up LocalGeometryFinder
    s1_finder = polyfinder.LocalGeometryFinder()
    s1_finder.setup_parameters(centering_type='standard', structure_refinement='none')
    s1_finder.setup_structure(structure)

    # Find site environment from LocalGeometryFinder
    environments = s1_finder.compute_structure_environments_detailed_voronoi(only_indices=[target_site],
                                                                             maximum_distance_factor=1.5)

    # Calculate actual site polyhedra using 'strategy'
    strategy = strategies.SimplestChemenvStrategy()
    strategy.set_structure_environments(environments)
    site_environment = strategy.get_site_coordination_environment(site=None, isite=target_site)

    return site_environment[0]


def find_species_string_ce(structure, species_string, min_fraction=0.0):
    """
    Returns mp_symbol's of all site coordination environments for the given species_string
    :param structure: (Structure) Structure containing target species
    :param species_string: (String) String representing the species we are looking for
    :param min_fraction: (float) minimum fraction that coordination environment of a site is calculated to be to be
    added to the list of site coordination environments
    :return: (List of Strings) List of mp_symbols of all site coordination environments for the given species_string
    """
    s1_finder = polyfinder.LocalGeometryFinder()
    s1_finder.setup_parameters(centering_type='standard', structure_refinement='none')
    s1_finder.setup_structure(structure)
    environments = s1_finder.compute_structure_environments_detailed_voronoi(maximum_distance_factor=1.5)

    light_se = se.LightStructureEnvironments(strategies.SimplestChemenvStrategy(), environments)

    all_ces = {}
    element = species_string
    for isite, site in enumerate(light_se.structure):
        if element in [sp.symbol for sp in site.species_and_occu]:
            for ce_dict in light_se._coordination_environments[isite]:
                if ce_dict['fraction'] < min_fraction:
                    continue
                if ce_dict['ce_symbol'] not in all_ces:
                    all_ces[ce_dict['ce_symbol']] = {'isites': [], 'fractions': [], 'csms': []}
                all_ces[ce_dict['ce_symbol']]['isites'].append(isite)
                all_ces[ce_dict['ce_symbol']]['fractions'].append(ce_dict['fraction'])
                all_ces[ce_dict['ce_symbol']]['csms'].append(ce_dict['csm'])
    return all_ces.keys()


def find_species_ce_from_light_se(light_structure, species_string, min_fraction=0.0):
    """
    Returns mp_symbol's of all site coordination environments for the given species_string
    :param structure: (Structure) Structure containing target species
    :param species_string: (String) String representing the species we are looking for
    :param min_fraction: (float) minimum fraction that coordination environment of a site is calculated to be to be
    added to the list of site coordination environments
    :return: (List of Strings) List of mp_symbols of all site coordination environments for the given species_string
    """

    all_ces = {}
    element = species_string
    for isite, site in enumerate(light_structure.structure):
        if element in [sp.symbol for sp in site.species_and_occu]:
            for ce_dict in light_structure._coordination_environments[isite]:
                if ce_dict['fraction'] < min_fraction:
                    continue
                if ce_dict['ce_symbol'] not in all_ces:
                    all_ces[ce_dict['ce_symbol']] = {'isites': [], 'fractions': [], 'csms': []}
                all_ces[ce_dict['ce_symbol']]['isites'].append(isite)
                all_ces[ce_dict['ce_symbol']]['fractions'].append(ce_dict['fraction'])
                all_ces[ce_dict['ce_symbol']]['csms'].append(ce_dict['csm'])
    return all_ces.keys()