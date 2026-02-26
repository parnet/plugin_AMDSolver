# AMDSolver Plugin #

**AMDSolver Plugin** is a **UG4 plugin** that provides a simple and rudimentary interface to use the **rocALUTION** solver library together with UG4.

The plugin bundles a variety of algebraic solvers, including both **direct** and **iterative** approaches:

## Included solver types

### Direct solvers
* LU decomposition
* Matrix inversion
* QR factorization

### Iterative solvers for sparse systems
* BiCGStab
* Conjugate gradients
* Algebraic multigrid (AMG)

These solvers enable efficient assembled systems within UG4 workflows.

For detailed information about solver algorithms, configuration options, and backend capabilities, please refer to the official [rocALUTION documentation](https://rocm.docs.amd.com/projects/rocALUTION/en/latest/).

---

## Project status

 **Experimental draft**

This plugin currently provides a first draft interface.  
APIs, naming, configuration mechanisms, and usage patterns are expected to evolve and may change significantly.

Limitation: The rocalution library requires a initialization and finish (compared to MPI). since there is no workflow to finish the usage of a plugin the 


## Installation
Please install or clone this repository directly into the plugin folder of a ug4 installation since it is not yet included into the package manager.
Before building just enable this with -DAMDSolver=ON 

To install the libraries needed use: 

on Fedora
```bash 
sudo dnf install rocalution rocalution-devel
```
if needed also
```bash
sudo dnf install rocm-hip hipblas hiprand rocsparse rocblas hipcub rocprim
sudo dnf install cmake gcc-c++ git openmp-devel
```

on Ubuntu refer to the official installation guide
https://rocm.docs.amd.com/projects/install-on-linux/en/latest/install/quick-start.html 