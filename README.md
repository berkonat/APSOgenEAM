# APSOgenEAM
EAM Potential Generator with Adaptive Particle Swarm Optimization (APSO) using MPI Parallelization

### How to Compile
gfortran APSOgenEAM.f -o APSOgenEAM.x

#### Example command to compile the code:
mpif77 APSOgenEAM.f -o APSOgenEAM.x -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -g

### How to Run
./APSOgenEAM.x < input_parameters

### Library Dependencies
- Intel MKL
