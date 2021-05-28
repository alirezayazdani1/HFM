default:
	cd gs && make ARCH=Linux mopt
	cd Veclib && make ARCH=Linux OPTM=1
	cd metis-4.0.3 && make
	cd Hlib/Linux && make ARCH=Linux PARALLEL=1 mopt && \
	cp ../../gs/libgs.a . && \
	cp ../../Veclib/libvec.a . && \
	cp ../../metis-4.0.3/libmetis.a .
	cd Nektar3d/Linux && make ARCH=Linux ADR=1 PARALLEL=1 mopt
	cp Nektar3d/Linux/Nektar* .
	
clean:
	cd gs && rm -f *.a
	cd Veclib && rm -f *.a
	cd metis-4.0.3 && make clean && rm -f *.a
	cd Hlib/Linux && rm *.a *.tab.c
	cd Nektar3d/Linux && make clean && rm Nektar*
