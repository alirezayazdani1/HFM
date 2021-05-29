# Dependencies
1. C, C++ and Fortran compilers (preferably GCC compilers Release 7 and lower or Intel Compilers 2018 and lower)
2. MPI library built by the compiler suite of step 1 (Open-MPI library is recommended)
3. LAPACK and BLAS libraries built by the compiler suite of step 1

# Compilation Instruction
1. **Download and install Open-MPI.** Download Open-MPI source code from https://www.open-mpi.org/software/ompi/v4.1/ and unzip the files. In the unzipped folder, run the configure script as follows (assuming a Linux machine with GCC compilers):

```bash
sudo  ./configure --prefix=/usr/local/openmpi CC=gcc CXX=g++ F77=gfortran FC=gfortran
sudo make all install
```
Make sure the compilers installation directory are in your default `$PATH` environment variable.

2. **Download and install LAPACK and BLAS libraries.** LAPACK and BLAS can be cloned from its GIT repo https://github.com/Reference-LAPACK/lapack and installed using the same compilers. Make sure the built libraries can be found and linked by adding their location to `$LD_LIBRARY_PATH` environment variable.

Note that on Linux machines it is possible to download and install binaries, lib and include files of the above dependencies using `yum` or `apt-get`.

3. **Modify the flag file** in the `SourceCode/Flags` directory depending on your system's architecture. For a local Linux workstation use `Linux.inc` file. 
You need to set the correct path to the mpi compiler as well as the math libraries to be linked (i.e., LAPACK and BLAS).

4. **Compile each of Nektar dependency library** and create the `*.a` static libraries:

- **gs:** Verify the Makefile and run `make mopt`.
- **Veclib**: Same as above then run `make OPTM=1`.
- **metis-4.0.3:** Modify `Makefile.in` accordingly (if needed) and run `make`.
- **Hlib:** Depending on your machine go to the corresponding directory (`Linux` or `Darwin` for MacOS) in `Hlib`. Double check the Makefile and MakeHybrid files to see everything is fine or make any changes needed. Then run `make PARALLEL=1 mopt`.
  - Copy all the built static libraries (`libgs.a, libvec.a, libmetis.a`) to `Hlib/Linux`.
- **Fluid Solver:** In `Nektar3d` go to the right directory (`Linux` or `Darwin`). Double check the Makefile and MakeNek files to see everything is fine. Run `make PARALLEL=1 ADR=1 mopt`. `ADR` flag adds advection-diffusion-reaction module to the fluid solver, which is needed for the [BenchmarkProblems](https://github.com/alirezayazdani1/HFM/tree/master/BenchmarkProblems); don't use this flag if you want to compile the fluid solver only.

If everything goes well, then you will have the executable: `Nektar.adr`.

Note that all the above steps are also inclueded in the `Makefile` in `HFM/SourceCode/` directory for convenience. Just run `make` to build the static libraries and the final executable.

# Running the Simulation

Put the `.rea` and `.adr` files in a directory say `RUN` and the executable next to it. Then run the simulation using `mpirun -np N ./Nektar.adr -V -chk -S -n3 RUN/Cylinder2D.rea`, where N is the number of processors.

The flags could be different depending on the specific problems. The `-V` flag is for time-dependent and position-dependent boundary conditions; `-chk` is for dumping results every IO checkpoint (defined in the `.rea` file) and `-S` is for every processor dumping its own results; `-n` is for number of polynomial modes.

# Results and Visualization

1. You will need the `p2sfld_s` and `nek2tec` executables for creating restart `RST` files and `dat or plt` files for tecplot. Their source codes can be found in the `Utilities/src` directory.

2. Go to `Utilities/Linux` directory, and after checking the `Makefile and MakeUtils` for compatibility run `make p2sfldmopt PARALLEL=1`. 
This will give the `p2sfld` executable to convert the IO files to `RST` files readable by [ParaView](https://www.paraview.org/) using the [pvNektarReader-plugin](https://github.com/alirezayazdani1/pvNektarReader-Plugin). Also, to convert the `RST` files to tecplot readable data format you need the `nek2tec` executable. Run `make nek2tecmopt PARALLEL=1` to compile that file as well.
