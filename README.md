#PHYS4061 - MgO & Ne crystal modeling

This project can simulate crystal structure and apply periodic boundary conditions to better simulate
the behaviour of individual atoms in a lattice structure. This allows for relaxation through conjugate gradient and steepest descent methods.


Dependencies:

* numpy
* matplotlib


This repository holds python code for modelling crystals. The crystal structure is displayed using [VMD](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD) by importing the xyz files. It is recommended to change the graphics options to CPK to properly view the structure.

Lab 4 code comprises the finsihed project which can take a lattice of atoms, randomly perturb them and then relax them into their equilibrium position using either conjugate gradient or steepest descent algorithms. Example structures for before and after relaxation can be fond in the lab4 folder.


The parameters used are as shown below:

Lennard Jones Parameters

|Interaction        | sigma (Angstroms) |epsilon (eV)       |
| -----------       | -------------     |:-------------:    |
| Ne <-> Ne         | 2.782             | 3.084 x 10^(-3)   |


Source:
https://journals.aps.org/prb/abstract/10.1103/PhysRevB.65.014112
Isotope effects in structural and thermodynamic properties of solid neon
Herrero, Carlos P.
American Physical Society. p. 014112. 2001




Buckingham Parameters

|Interaction        | A (eV)        |rho (Angstroms) |C (Angstrom^(6)eV)|
| -----------       | ------------- |:-------------:|   ---------------|
| Mg 1.7+ <-> O 1.7-| 926.69        | 0.29909       | 0                |
| O 1.7- <-> O 1.7- | 4870.0        | 0.2670        | 77.0             |


Source:
https://journals.aps.org/prb/abstract/10.1103/PhysRevB.72.115437
MgO addimer diffusion on MgO(100): A comparison of ab initio and empirical models
Henkelman, Graeme and Uberuaga, Blas P. and Harris, Duncan J. and Harding, John H. and Allan, Neil L.
American Physical Society. p. 115437. 2005
