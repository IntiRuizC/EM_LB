# Lattice-Boltzmann for Electrodynamics

This project implements the 3D lattice-Boltzmann model for the Maxwell equations in media by M. Mendoza and J. D. Muñoz (arXiv:0806.2678). The test case is an electromagnetic Gaussian pulse crossing the interface between vacuum and a dielectric with relative permittivity 2.5.

## Requirements

- **C++ compiler**
- **MPI** (`mpic++`, `mpirun`)
- **Make utility**

## Running

```
make run NP=4
```

The lattice is split in slabs along z, one per MPI rank; output goes to `pulse.dat`.
