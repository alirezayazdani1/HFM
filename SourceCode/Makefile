ARCH := $(shell uname)

default:
	cd gs && make mopt
	cd Veclib && make OPTM=1
	cd metis-4.0.3 && make
	cd Hlib/$(ARCH) && make PARALLEL=1 mopt && \
	cp ../../gs/libgs.a . && \
	cp ../../Veclib/libvec.a . && \
	cp ../../metis-4.0.3/libmetis.a .
	cd Nektar3d/$(ARCH) && make ADR=1 PARALLEL=1 mopt
	cp Nektar3d/$(ARCH)/Nektar* .
	
clean:
	cd gs && rm -f *.a
	cd Veclib && rm -f *.a
	cd metis-4.0.3 && make clean && rm -f *.a
	cd Hlib/$(ARCH) && rm *.a *.tab.c
	cd Nektar3d/$(ARCH) && make clean && rm Nektar*
