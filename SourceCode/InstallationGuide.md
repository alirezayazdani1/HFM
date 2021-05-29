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
2. **Download and install LAPACK and BLAS libraries.** LAPACK and BLAS can be cloned from its GIT repo https://github.com/Reference-LAPACK/lapack and installed using the same compilers.

Note that on Linux machines it is possible to download and install binaries, lib and include files of the above dependencies using `yum` or `apt-get`.

3. **Modify the flag file** in the `Flags` directory depending on your system's architecture. For a local Linux workstation use `Linux.inc` file. 
You need to set the correct path to the mpi compiler as well as the math libraries to be linked (i.e., LAPACK and BLAS).

4. **Compile each of Nektar dependency library** and create the `*.a` static libraries:

- **gs:** Verify the Makefile and run `make mopt`.
- **Veclib**: Same as above then run `make OPTM=1`.
- **metis-4.0.3:** Modify `Makefile.in` accordingly (if needed) and run `make`.
- **Hlib:** Depending on your machine go to the corresponding directory (`Linux` or `Darwin` for MacOS) in `Hlib`. Double check the Makefile and MakeHybrid files to see everything is fine or make any changes needed. Then, run `make PARALLEL=1 mopt`.
- **Fluid Solver:** In `Nektar3d` go to the right directory (`Linux` or `Darwin`). Double check the Makefile and MakeNek files to see everything is fine. Run `make PARALLEL=1 ADR=1 mopt`.

If everything goes well, then you will have the executable: `Nektar.adr`

# Running the Simulation

1. Put the .rea file in a directory say "RUN". Then run the simulation by typing: "mpirun -np N ./Nektar -V -chk -S -n4 RUN/fluid.rea".

where N is the number of processors. The flags could be different based on different problems. The -V flag is for time-dependent and position-dependent boundary condition; 
-chk is for dumping results every IO checkpoint and -S is for every processor dumping its own results; -n is for number of modes.

Visualization:

1. You need the "p2sfld_s and nek2tec" executables for creating restart "RST" files and "dat/plt" files for tecplot. Source files can be found in the Utilities/src directory.

2. Go to "Utilities/Linux" directory, and after checking the "Makefile and MakeUtils" for compatibility issues type "make p2sfldmopt". 
This will give the "p2sfld" executable to convert the IO files to "RST" files readable by ParaView using the pvNektar-plugin. Also to convert the "RST" files to tecplot readable data format you need the "nek2tec" executable. Type "make nek2tecmopt" to compile that file as well.
