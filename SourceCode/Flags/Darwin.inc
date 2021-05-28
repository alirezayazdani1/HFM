COMPATH	= /opt/openmpi/
CXX	= g++-7
CC	= gcc-7 
ifeq ($(MACHTYPE),i386)
FC	= /sw/bin/g95
FFLAGS = -DNO_DLAMCH
CPPFLAGS := $(CPPFLAGS) -D__APPLE_INTEL__
else
FC 	= gfortran-7
endif

MPILIB	= -L$(COMPATH)/lib
MPIINC	= -I$(COMPATH)/include

ifdef MPIPRG
# CFM default rule
message	+= CFM default uses MPICH, ch_p4, sun4u. 
CXX	= $(COMPATH)/bin/mpic++
CC	= $(COMPATH)/bin/mpicc
FC	= $(COMPATH)/bin/mpifort
endif
message	+= You can specify MCPU=-mcpu=pentiumpro,ev5 etc or MCPU=-march=pentium. See the gcc info pages for more information. 
OPTXX	= -O3 $(MCPU) -ffast-math 
MOPTXX	= -O3 $(MCPU) -ffast-math -funroll-loops -fstrict-aliasing 
OPT	= -O3 $(MCPU) -ffast-math 
MOPT	= -O3 $(MCPU) -ffast-math -funroll-loops -fstrict-aliasing 
OPTF	= -O3 $(MCPU) -ffast-math 
MOPTF	= -O3 $(MCPU) -ffast-math -funroll-loops  
DEBUG	= -g  -gstabs+
DEBUGF	= -g $(FPPFLAGS)
LD	= $(CXX) -bind_at_load
DEBUGLD	= $(DEBUG)
MOPTLD	= $(MOPTXX)
OPTLD	= $(OPTXX)
ifdef FFTPRG
FFTLIB	= -lfftw
endif
LLAPACK = #-L/usr/local/opt/lapack/lib -llapack
LBLAS   = #-L/usr/local/opt/openblas/lib -lblas
LIBS    = $(FFTLIB) $(LLAPACK) $(LBLAS) -framework Accelerate
