# APSOgenEAM
EAM Potential Generator with Adaptive Particle Swarm Optimization (APSO) using MPI Parallelization

This code generates DYNAMO style tabulated EAM files. These files can also be used by LAMMPS MD code.

Initially, this code is written to develop Cu-Ni alloy potential as it is detailed in 
[Berk Onat and Sondan Durukanoğlu 2014 J. Phys.: Condens. Matter 26 035404](http://iopscience.iop.org/article/10.1088/0953-8984/26/3/035404).

As it is stated in the LICENCE of this code, any publication that uses this code and/or the results of this code should cite the paper: [Berk Onat and Sondan Durukanoğlu 2014 J. Phys.: Condens. Matter 26 035404](http://iopscience.iop.org/article/10.1088/0953-8984/26/3/035404)

### How to Compile
mpif77 APSOgenEAM.f -o APSOgenEAM.x

#### Example command to compile the code:
mpiifort APSOgenEAM.f -o APSOgenEAM.x -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -g

### How to Run
./APSOgenEAM.x < input_parameters

### Library Dependencies
- Intel MKL
