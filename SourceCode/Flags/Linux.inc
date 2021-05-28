CXX	= g++
CC	= gcc 
FC	= gfortran

COMPATH	= /usr/lib64/openmpi
MPILIB	= -L$(COMPATH)/lib
MPIINC	= -I$(COMPATH)/include

ifdef MPIPRG
message	+= CFM default uses OpenMPI
CXX	= $(COMPATH)/bin/mpic++
CC	= $(COMPATH)/bin/mpicc
FC	= $(COMPATH)/bin/mpifort
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
LLAPACK = -L/usr/lib64 -llapack
LBLAS   = -L/usr/lib64 -lblas
LIBS	  = $(FFTLIB) $(LLAPACK) $(LBLAS)
