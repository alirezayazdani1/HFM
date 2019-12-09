COMPATH	= 
CXX	= g++
CC	= gcc 
FC	= gfortran

MPILIB	= -L/usr/local/packages/mpich2-gcc/lib
MPIINC	= -I/usr/local/packages/mpich2-gcc/include/

ifdef MPIPRG
# CFM default rule
message	+= CFM default uses MPICH, ch_p4, sun4u. 
CXX	= mpicxx
CC	= mpicc
FC	= gfortran
endif
message	+= You can specify MCPU=-mcpu=pentiumpro,ev5 etc or MCPU=-march=pentium. See the gcc info pages for more information. 
OPTXX	= -O3 $(MCPU) -ffast-math  -w
MOPTXX	= -O3 $(MCPU) -ffast-math -funroll-loops -fstrict-aliasing 
OPT	= -O3 $(MCPU) -ffast-math  -w
MOPT	= -O3 $(MCPU) -ffast-math -funroll-loops -fstrict-aliasing 
OPTF	= -O3 $(MCPU) -ffast-math  -w
MOPTF	= -O3 $(MCPU) -ffast-math -funroll-loops  
DEBUG	= -g -w 
DEBUGF	= -g -w
LD	= $(CXX)
DEBUGLD	= $(DEBUG)
MOPTLD	= $(MOPTXX)
OPTLD	= $(OPTXX)
ifdef FFTPRG
FFTLIB	= -lfftw
endif
LLAPACK = -L/usr/local/packages/lapack-3.3.1 -llapack 
LBLAS   = -L/usr/local/packages/BLAS -lblas
LPG   	= #-L/usr/lib -lpgftnrtl -lpgc   
LIBS	= $(FFTLIB) $(LLAPACK) $(LBLAS) 
