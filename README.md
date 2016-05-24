# Structural Descriptors for Local Topology

The goal of this project is to develop structural descriptors for the local topology of given structures. Specifically,
we have descriptors describing the local coordination environment of cations and describing the connectivity between
different coordination polyhedra in the structure. The intention of these descriptors is firstly to aid people in
understanding the individual structures they're working with and secondly to act as descriptors for data mining
properties that may depend on local topology.


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing
purposes. See deployment for notes on how to deploy the project on a live system.

The structural_descriptors_repo repository contains the code for the both the local coordination number descriptor
and the local connectivity descriptor. It should be run with Python (preferably 2.7) and relies heavily on code from
Pymatgen. It also uses Numpy.

### Prerequisities

Updated version of Pymatgen including the Chemenv module (and all its dependencies, including Python and Numpy)

```
http://pymatgen.org/installation.html
```

### Installing

Once the code is downloaded from source and prerequisites are installed, the code should be immediately usable.

Download code from source

```
$ git clone https://github.com/tchen0965/structural_descriptors_repo.git
```

### Starting

The code can then be run to obtain connectivities of structures by editing and running run_connectivity.py.
Specifically, 'path' should be changed to the path to the structure being analyzed. The central and peripheral species
can be specified, or else the default is that only atoms (the anions O2-, F-, Cl-, I-, Br-, S2-, and N3- and atoms O, F,
Cl, I, Br, S, and N) will be treated as peripheral atoms for the polyhedra whose connectivities we are calculating. A
radius can also be satisfied if the target (central) atoms we're interested in are very large (i.e. Ba, Sr).

Ex: for LiCoO2 (which is already in the structural_descriptors_repo folder), if we want to find the connectivities
between all cations, we could leave the parameters as the following (as the default will consider O the anion and Li
and Co as the cations).

```
path = 'LiCoO2.cif'
central_species = None
peripheral_species = None
radius = 2.6
```
Running run_connectivity.py with these lines then gives us structure information and connectivity information:

```
Full Formula (Li1 Co1 O2)
Reduced Formula: LiCoO2
abc   :   2.844940   2.844940   4.969293
angles:  73.366304  73.366304  60.000000
Sites (4)
  #  SP           a         b         c
---  ----  --------  --------  --------
  0  Li    0.5       0.5       0.5
  1  Co    0         0         0
  2  O     0.739749  0.739749  0.780752
  3  O     0.260251  0.260251  0.219248

Co1 are 6-fold coordinated.
They are edge-connected to 6-fold coordinated Co1. They are edge-connected and point-connected to 6-fold coordinated Li1.

Li1 are 6-fold coordinated.
They are edge-connected and point-connected to 6-fold coordinated Co1. They are edge-connected to 6-fold coordinated Li1.
```

Ex: for LiCoO2, if we want to find the connectivities between only Li, we could specify the central_species to be only
Li.

```
path = 'LiCoO2.cif'
central_species = ['Li']
peripheral_species = None
radius = 2.6
```

Running run_connectivity.py with these lines then gives us the same structure information but different connectivity
information:

```
Li1 are 6-fold coordinated.
They are edge-connected to 6-fold coordinated Li1.
```

This can give less information, but the information can also be more easy to interpret for larger or more complex
structures.

## Running the tests

The automated tests for this project are given in the 'examples' folder. These tests test the code to calculate
coordination number (coordination_test.py) on some simple structures, the code to calculate connectivity
(connectivity_test.py) on some simple structures, and the code for both coordination number and connectivity on a variety
 of common structures (structures_test.py).

### Break down into end to end tests

The coordination number tests (coordination_test.py) test the accuracy of the O'Keeffe and effective coordination number
methods from effective_coordination.py and okeeffe_CN.py.

Ex: Finding the effective coordination number for a calcium site in CaF2
```
from pymatgen.core.structure import Structure
import effective_coordination as econ

caf2_structure = Structure.from_file('CaF2.cif', True, True)

# create an EffectiveCoordFinder object given the structure of interest
caf2_econ_finder = econ.EffectiveCoordFinder(caf2_structure)

for isite, site in enumerate(caf2_structure._sites):
    if site.species_string == 'Ca':
        first_site = site
        ifirst_site = isite

# calculate the coordination number of site using either the site itself or the index of the site for the
# structure
self.assertAlmostEqual(caf2_econ_finder.get_site_cn(first_site), 8, 1,
                      "Ca should be, on average, 8-fold coordinated in CaF2")
self.assertEqual(caf2_econ_finder.get_site_cn(first_site), caf2_finder.get_isite_cn(ifirst_site),
                "CN from site and CN from index of site should be the same")
```

The connectivity tests (connectivity_test.py) test the accuracy of the connectivity matrices built by the connectivity
code.

Ex: Finding the connectivity matrix of LiCoO2 (how Li polyhedra are connected to Co polyhedra and vice versa)

```
from pymatgen.core.structure import Structure
import connectivity_from_structure as connectivity

# the get_connectivity_matrix method gives several outputs, the first of which is the actual connectivity matrix, the
# second of which are the polyhedra in the structure, the third of which is the supercell structure created to calculate
# the connectivity matrix
licoo2_matrix, licoo2_polyhedra, supercell = \
            connectivity.get_connectivity_matrix(self.licoo2_structure, False)

# the keys of the matrix are the species_string of the central sites
self.assertIn('Li', licoo2_matrix.keys(), "Li not found in LiCoO2 matrix")
self.assertIn('Co', licoo2_matrix.keys(), "Co not found in LiCoO2 matrix")

# the matrix is given as nested dictionaries, with the first key being the first type of species, the second key being
the second type of species, and the third key being the type of connectivity-sharing between the two species
self.assertEqual(licoo2_matrix['Li']['Li']['point'], 0, "Li should not be point-sharing with other Li polyhedra")
self.assertEqual(licoo2_matrix['Li']['Li']['edge'], 6, "Li should be edge-sharing with other Li polyhedra 6 times")
self.assertEqual(licoo2_matrix['Li']['Li']['face'], 0, "Li should not be face-sharing with each other")
self.assertEqual(licoo2_matrix['Co']['Co']['point'], 0, "Co should not be point-sharing with each other")
self.assertEqual(licoo2_matrix['Co']['Co']['edge'], 6, "Co should be edge-sharing with 6 other Co polyhedra")
self.assertEqual(licoo2_matrix['Co']['Co']['face'], 0, "Co should not be face-sharing with each other")
self.assertEqual(licoo2_matrix['Li']['Co']['point'], 6, "Li should be point-sharing with 6 Co polyhedra")
self.assertEqual(licoo2_matrix['Li']['Co']['edge'], 6, "Li should be edge-sharing with 6 Co polyhedra")
self.assertEqual(licoo2_matrix['Li']['Co']['face'], 0, "Li and Co polyhedra should not be face-sharing")
self.assertEqual(licoo2_matrix['Li']['Co'], licoo2_matrix['Co']['Li'],
                         "Connectivity matrix should be symmetric")
```

Ex: Finding the connectivity description of LiCoO2 (obtaining a verbal description of the connectivities between
polyhedra in LiCoO2)

```
from pymatgen.core.structure import Structure
import connectivity_from_structure as connectivity

# find the connectivity matrix of the structure
licoo2_matrix, licoo2_polyhedra, supercell = \
   connectivity.get_connectivity_matrix(self.licoo2_structure, True)

# give the connectivity matrix to the get_connectivity_description method
licoo2_descriptions = connectivity.get_connectivity_description(licoo2_matrix, licoo2_polyhedra, self.licoo2_structure, False)

# the output is a dictionary of descriptions, the keys of which are the cation species of interest and the values of
# which are the connectivity description for the polyhedra of that cation species
for cation in licoo2_descriptions.keys():
   in_matrix_keys = False
   for keys in licoo2_matrix.keys():
       if cation in keys:
           in_matrix_keys = True
   self.assertTrue(in_matrix_keys, "Cation being described not in connectivity matrix")
   self.assertTrue(isinstance(licoo2_descriptions[cation], str)
                  or isinstance(licoo2_descriptions[cation], unicode),
                  "Descriptions are not type str or unicode")
```

The structures test (structures_test.py) test the accuracy of the connectivity matrix method on a variety of common
structure types, including the barite structure, the fluorite structure, the pervoskite structure, the spinel structure,
and various other structures. The tests are almost identical to the connectivity matrix tests (in connectivity_test.py).



## Built With

* PyCharm


## Contributing

Please email Tina Chen at tchen0965@gmail.com for information on how to contribute.


## Authors

* Tina Chen - *Initial work* - [tchen0965](https://github.com/tchen0965)


## License

This project is licensed under the MIT License (see License.md for details).


## Acknowledgments

* Materials Project
* Pymatgen
* Lawrence Berkeley National Lab
* UC Berkeley
* Anubhav Jain
