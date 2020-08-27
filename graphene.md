# Graphene 
A tutorial for the calculation of the band structure of graphene via the *density functional theory*  
In this tutorial you will learn how to do the following: 
- define the crystal structure of graphene 
- set up the essential first-principles parameters
- perform a self consistent field calculation 
- test the convergence of total energy with respect to the cut-off energy and the number of k points. 


## Set the path 
All the necessary files are available in `graphene` folder. First change the directory to `graphene`. 
```bash
cd graphene/
```
Then you need to set the path to the root folder of Quantum ESPRESSO on your machine. You can do this by editing the `qe-path` file: 
```bash
#!/bin/sh

QE_PATH=$HOME/code/espresso6.0 # replace with your quantum espresso path
```

## Input files
There are various input files that we use to feed data and parameters to the Quantum ESPRESSO executables. 
Here are the descriptions: 
- `run-graphene`: this is the main `bash` script that we use to run the executables and to set and change the parameters of the calculations. 
- `in/graphene.scf`: the input file related to the *self consistent field* calculations that is used to define the crystal structure and the calculations type and other parameters. The code in this file are mainly in `fortran` language which Quantum ESPRESSO is written in. Therefore, a little bit of familiarity with `fortran` could greatly help in resolving bugs and issues with the input files of this sort.
- `in/graphene.vc`: very similar input file to the `.scf` with the same definition of the crystal structure except for a different type of calculation which is called *variable cell* calculation. 
- `pseudo/`: this folder holds the *pseudo-potential* files which for now we consider them as a smoothed approximation to the real Coulomb potential of the atomic nuclei. We will discuss these in depth in the future. 

**Note**: for a detailed description of input files refer to the Quantum ESPRESSO documentations. You can find them in `PW/Doc/INPUT_PW.html` which is located in the root directory of Quantum ESPRESSO. 

## Define the crystal structure 
To define the crystal structure of graphene you need to know the crystal system and the relative atomic positions. In Quantum ESPRESSO you can easily choose your [*Bravais* lattice](https://en.wikipedia.org/wiki/Bravais_lattice) and other related parameters. The following lines in the `graphene.scf` file in the `&system` section defines the crystal structure. 
```fortran
ibrav = 4           ! index of the Bravais lattice (here a hexagonal lattice)
A = 2.46000000      ! lattice constant in angstrom 
C = 20.00000000     ! vacuum size along the c axis 
nat = 2             ! number of atoms in the unit cell (2 carbon atoms)
ntyp = 1            ! different types of atoms (only 1: carbon)
```
**Note**: since these lines of code are in `fortran` language, the exclamation mark `!` stands for commenting.  
**Note**: The reason that the value for `C` is much larger that the lattice constant is that all the structures in Quantum ESPRESSO are considered periodic. Therefore, to define graphene which is only 1-atom thick and is not periodic in the out of plane direction, we need to introduce a large periodicity. Usually `20` Angstrom is reasonable. 

After defining the crystal structure you will need to determine the `ATOMIC_SPECIES`. 
```fortran 
ATOMIC_SPECIES
 C 12.0106000  $pseudo_c
```
The first entry determines the element carbon. The second entry is the atomic weight which you can find from [NIST](https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses). The third entry is the name of the pseudo-potential you are using which is located in `pseudo/` folder. Here we use a bash placeholder `$pseudo_c` which is in fact defined in the `run-graphene` file as follows.
```bash
pseudo_c='C.pbe-n-kjpaw_psl.1.0.0.UPF'
```

Next, you need to enter the atomic positions.
```fortran
ATOMIC_POSITIONS {angstrom}
 C 0.00000000 -0.710140831  0.00000000
 C 0.00000000  0.710140831  0.00000000
```
The first entry is the atomic symbol. The next three are the coordinates of the atom in the unit cell in Angstrom units. Graphene has two atoms and each is defined in a separate line. 
**Note**: in the case of graphene the second coordinate `y=0.710140831` of atomic positions is dependent on the lattice constant. In fact, `y = A/(2*sqrt(3))`.


## Set up main parameters
The main parameters that you may want to change at different runs are listed in the beginning of the `run-graphene` files. 
```bash
#!/bin/sh

# parameters
prefix='graphene'                       # used to name output files
pseudo_c='C.pbe-n-kjpaw_psl.1.0.0.UPF'  # pseudo-potential file 
ecutwfc=40.0                            # kinetic energy cutoff
ecutrho=326.0                           # charge density energy cutoff
num_kx=8                                # number of k point along x
num_ky=8                                # number of k point along y
num_kz=1                                # number of k point along z   
num_bands=16                            # number of bands 
```
**Note**: since this file is written in `bash`, the hash symbol `#` is used for commenting.  
**Note**: when defining variables in `bash` there should be no space before and after the equal sign. Therefore, `a = 5` won't work and the correct command is `a=5` with no spaces.  
**Note**: The recommended values of kinetic energy cutoff and the charge density cutoff are usually mentioned in the beginning of the pseudo-potential files.  
**Note**: The three-dimensional k mesh will be of size `num_kx X num_ky X num_kz`


## Perform the self-consistent field calculation 
As you may recall from the *density functional theory*, the Schrodinger's equation depends on an effective potential that depends on the charge density which is the wavefunction squared. Therefore, the Schrodinger equation needs to be solved self consistently. That is what the `scf` calculations does. You can perform this calculation by executing the `run-graphene` script inside the `graphene/` folder.  
```bash
bash run-graphene
```
When the calculation is done you will see `done` in the terminal window. The output of the calculation is written in the following file `out/graphene.scf.out`. This file is usually verbose. Using the following commands you can get the important information about your system quickly. 
```bash
grep -e "\!" out/graphene.scf.out
```
This command gives you the total energy (in Ry units) of the system at the end of self consistent calculations.  

```bash
grep -e "Total force" out/graphene.scf.out
```
This command gives you the total force (in Ry/au units) of the system. 

## Convergence tests 
You can test the convergence of the calculations by redoing the calculations with a different set of parameters and plotting the following.  
1. Total energy vs. kinetic energy cutoff: make sure the total energy converges as the kinetic energy cutoff increases.  
2. Total energy vs. number of k points: make sure the total energy converges as the number of k points increases, e.g., from `4X4X1` to `5X5X1` and etc.  
3. Total energy vs. the lattice constant `A`: does the minimum energy happen for `A=2.46` Angstrom or some other number? Redo the calculations for different values of `A` and see how the energy behaves as a function of `A`.  
**Note**: make sure to update the `y` coordinate of the `ATOMIC_POSITIONS` as you change `A`. 