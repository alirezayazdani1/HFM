default:
	cd gs && make ARCH=Darwin mopt
	cd Veclib && make ARCH=Darwin OPTM=1
	cd metis-4.0.3/Lib && make
	cd Hlib/Darwin && make ARCH=Darwin PARALLEL=1 mopt && \
	cp ../../gs/libgs.a . && \
	cp ../../Veclib/libvec.a . && \
	cp ../../metis-4.0.3/libmetis.a .
	cd Nektar3d/Darwin && make ARCH=Darwin ADR=1 PARALLEL=1 mopt
	cp Nektar3d/Darwin/Nektar* .
	
clean:
	cd gs && rm -f *.a
	cd Veclib && rm -f *.a
	cd metis-4.0.3 && make clean && rm -f *.a
	cd Hlib/Darwin && rm libgs.a libvec.a libmetis.a libhybridmopt.a *.tab.c
	cd Nektar3d/Darwin && make clean
