ARCH = $(shell uname)
MACH    := $(shell uname -m)
PROG    = Nektar
VPATH   = ../src_adr
ifdef FCM
VPATH   = ../src_bio
endif
LIBDIR  = ../../Hlib/$(ARCH)
HDRDIR  = ../../include
TARFILE = $(PROG).tar

ifeq ($(ARCH),Linux)
  ifeq ($(MACH),ia64)
        ARCH =  Linux-ia64
  endif
endif

defines = METIS PBC_1D NEK2VTK
CPPFLAGS = -I$(HDRDIR)
ifdef PARALLEL
MPIPRG  = yes
override defines := PARALLEL $(defines)
endif

ALL     = 	analyser    bwoptim     cfl  convective     comm  	\
		dgalerkin   drive       forces              io          \
		mlevel 	    prepost     pressure            rotational	\
		sections    stokes      subcycle     				womersley   \
		flowrate 		pbc_1d 			small_vessel				out_vtkXML

ifdef VELINTERP
override SPECIAL := velinterp $(SPECIAL)
endif

ifdef ADR
adr			=		bcflux chemical
override SPECIAL := $(adr) $(SPECIAL)
endif

ifdef FCM
fcm 		=		fcm_bgmesh	fcm_boundary	fcm_domain	fcm_particle	fcm_bio
override	SPECIAL := $(fcm) $(SPECIAL)
endif

include ../../Flags/$(ARCH).inc

ifdef ADR
override ADDONS := $(ADDONS) ADR
PROG = Nektar.adr
endif

ifdef FCM
override ADDONS := $(ADDONS) FCM
PROG = Nektar.fcm
ifdef ADR
PROG = Nektar.fcm.adr
endif
endif

ifdef CCHF
override ADDONS := $(ADDONS) CCHF
endif

ifdef CCLF
override ADDONS := $(ADDONS) CCLF
endif

ifdef defines
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
