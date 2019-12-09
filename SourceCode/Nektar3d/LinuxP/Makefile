ARCH = $(shell uname)
MACH    := $(shell uname -m)
PROG    = nektar
VPATH   = ../src
LIBDIR  = ../../Hlib/$(ARCH)P
HDRDIR  = ../../include
TARFILE = $(PROG).tar
MPIINC  = /usr/local/packages/mpich2-gcc/include
MPIPRG  = yes

ifeq ($(ARCH),Linux)
  ifeq ($(MACH),ia64)
        ARCH =  Linux-ia64
  endif
endif

defines = METIS PARALLEL
CPPFLAGS = -I$(HDRDIR)
ifdef PARALLEL
MPIPRG  = yes
endif

ALL     = 	analyser    bwoptim     cfl  convective     comm   	\
		dgalerkin   drive       forces              io          \
		mlevel 	    prepost     pressure            rotational	\
		sections    stokes      subcycle     	    womersley   \
		velinterp   

include ../../Flags/$(ARCH).inc

ifeq ($(ARCH),cray-pvp)			# Cray PVP
include ../../Flags/CRAYPVP.inc
endif

ifeq ($(ARCH),cray-t3d)			# Cray T3D
include ../../Flags/T3D.inc
endif

ifeq ($(ARCH),unicosmk)			# Cray T3E
include ../../Flags/T3E.inc
endif

ifeq ($(ARCH),alpha)			# DEC Alpha AXP
include ../../Flags/DECAXP.inc
endif

ifeq ($(ARCH),OSF1)                     # DEC Alpha AXP
include ../../Flags/DECAXP.inc
endif

ifeq ($(ARCH),XT3)                      # CRAY  XT3
include ../../Flags/XT3.inc
endif


ifdef PARALLEL
ifeq ($(ARCH), AIX)
CXX = $(MPICXX)
endif
override CPPFLAGS :=  -DPARALLEL  $(CPPFLAGS) $(MPIINC)

ifneq ($(ARCH), XT3)
override LIBS     := $(LIBS) -lmetis -lgs $(MPILIB) -lmpi
endif

ifeq ($(ARCH), XT3)
override LIBS     := $(LIBS) -lmetis -lgs $(MPILIB)
endif


ifeq ($(ARCH),ap3000)
override LIBS     := $(LIBS) -L/opt/FSUNaprun/lib -lmpl -lemi -lthread 
endif
ifeq ($(ARCH),OSF1)
override LIBS     := $(LIBS) -lelan -lelan3
endif
endif

ifdef WOMERR
override ADDONS := WOMERR
PROG = nektar.wom
endif

ifdef VELINTERP
override ADDONS := VELINTERP
PROG = nektar.vel
endif

ifdef    defines
DEF     = $(foreach def, $(defines) $(ADDONS), -D$(def))
override CPPFLAGS :=  $(DEF) $(CPPFLAGS)
endif

export


help:
	@echo ''
	@echo ' Nektar Makefiles:'
	@echo ''
	@echo ' Run make with one of the following options:'	
	@echo '	a) dbx.... dbx version'
	@echo '	b) opt.... optimised version'
	@echo ''


dbx: 
	$(MAKE) LTYPE=g -f MakeNek
opt:
	$(MAKE) LTYPE=opt -f MakeNek
mopt:
	$(MAKE) LTYPE=mopt -f MakeNek

tar: Makefile
	@echo 'creating a tape archive in' $(TARFILE)
	@tar cvf $(TARFILE) Makefile MakeNek src/*.C
clean:
	$(MAKE) -f MakeNek clean


