[![DOI](https://zenodo.org/badge/199072669.svg)](https://zenodo.org/badge/latestdoi/199072669)

# Hidden Fluid Mechanics
We present hidden fluid mechanics (HFM), a physics informed deep learning framework capable of encoding an important class of physical laws governing fluid motions, namely the Navier-Stokes equations. In particular, we seek to leverage the underlying conservation laws (i.e., for mass, momentum, and energy) to infer hidden quantities of interest such as velocity and pressure fields merely from spatio-temporal visualizations of a passive scaler (e.g., dye or smoke), transported in arbitrarily complex domains (e.g., in human arteries or brain aneurysms). Our approach towards solving the aforementioned data assimilation problem is unique as we design an algorithm that is agnostic to the geometry or the initial and boundary conditions. This makes HFM highly flexible in choosing the spatio-temporal domain of interest for data acquisition as well as subsequent training and predictions. Consequently, the predictions made by HFM are among those cases where a pure machine learning strategy or a mere scientific computing approach simply cannot reproduce. The proposed algorithm achieves accurate predictions of the pressure and velocity fields in both two and three dimensional flows for several benchmark problems motivated by real-world applications. Our results demonstrate that this relatively simple methodology can be used in physical and biomedical problems to extract valuable quantitative information (e.g., lift and drag forces or wall shear stresses in arteries) for which direct measurements may not be possible.

# Source Codes:
## Nektar-3D
3D spectral/hp element solver for unsteady incompressible Navier-Stokes equations coupled with advection-diffusion-reaction equations.

Please refer to [Installation Guide](https://github.com/alirezayazdani1/HFM/blob/master/SourceCode/InstallationGuide.md) for compilation and building Utility tools (the code has been successfully compiled with Intel Compilers 2018 and GCC 7.5).

## Nektar++
Nektar++ is a tensor product based finite element package designed to allow one to construct efficient classical low polynomial order h-type solvers (where h is the size of the finite element) as well as higher p-order piecewise polynomial order solvers.

Nektar++ is available both as pre-compiled binaries for some operating systems and as source-code tarballs at the following address: https://www.nektar.info/downloads/

# Benchmark Problems
The benchmark problems include the input files needed for the simulations to run. These input files contain the physical and solver parameters, grid, initial and boundary conditions. The input file for Nektar++ is in xml format; inputs for Nektar-3D are ASCII data files.
