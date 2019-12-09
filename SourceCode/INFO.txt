The follwoing is a short instruction manual to compile Nektar and run a simulation.

Compilation instruction:

1. Modify the flag file in the Flags directory depending on your system's architecture. For a local workstation use "Linux.inc" file. 
You need to set the correct path to the mpi compiler as well as the math libraries to be linked (e.g., Lapack, Blas, ACML or MKL math libraries...).

2. Start compiling each of the libraries first and create the *.a static libraries; by executing "Make" in the main directory the following will be compiled:

--gs: verify the Makefile for the architecture variable ARCH (should be Linux if your flag file is Linux.inc) and run "make mopt".
--Veclib: Same as above then type "make OPTM=1"
--metis-4.0.3: Modify "Makefile.in" accordingly. Go to "Lib" directory and type "make".
--Hlib: Depending on your machine go to the corresponding directory (e.g. Linux). Double check the Makefile for the "defines" first. Then run "make PARALLEL=1 mopt".
--3D Fluid Solver: Go to "Nektar3d" directory. Double check the Makefile and MakeNek files to see everything is fine. Run "make PARALLEL=1 mopt".

If everything goes well, then you will have the executable: "Nektar".

Running a simulation:

5. Put the .rea file in a directory say "RUN". Then run the simulation by typing: "mpirun -np N ./Nektar -V -chk -S -n4 RUN/fluid.rea".

where N is the number of processors. The flags could be different based on different problems. The -V flag is for time-dependent and position-dependent boundary condition; 
-chk is for dumping results every IO checkpoint and -S is for every processor dumping its own results; -n is for number of modes.

Visualization:

1. You need the "p2sfld_s and nek2tec" executables for creating restart "RST" files and "dat/plt" files for tecplot. Source files can be found in the Utilities/src directory.

2. Go to "Utilities/Linux" directory, and after checking the "Makefile and MakeUtils" for compatibility issues type "make p2sfldmopt". 
This will give the "p2sfld" executable to convert the IO files to "RST" files readable by ParaView using the pvNektar-plugin. Also to convert the "RST" files to tecplot readable data format you need the "nek2tec" executable. Type "make nek2tecmopt" to compile that file as well.
