CXX   = mpicxx
CC    = mpicc
FC    = mpif90

#MPILIB  = #-L/opt/apps/intel13/modulefiles/mvapich2/1.9a2/lib
#MPIINC  = #-I/opt/apps/intel13/modulefiles/mvapich2/1.9a2/include

OPTXX	= -O3 -g
MOPTXX	= -O3 -g

OPT	= -O3 -g
MOPT	= -O3 -g

OPTF	= -O3 -g
MOPTF	= -O3 -g

DEBUG	= -g 
DEBUGF	= -g 

LD	= $(CXX)
DEBUGLD	= $(DEBUG)
MOPTLD	= $(MOPTXX)
OPTLD	= $(OPTXX)

ifdef FFTPRG
FFTLIB	= -lfftw
endif

LIBS = -L$(MKLROOT)/lib/intel64 -lmkl_core -lmkl_sequential -lmkl_intel_lp64
