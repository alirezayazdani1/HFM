default:
	cd gs && make ARCH=STAMPEDE mopt
	cd Veclib && make ARCH=STAMPEDE OPTM=1
	cd metis-4.0.3 && make
	cd Hlib/STAMPEDE && make ARCH=STAMPEDE PARALLEL=1 mopt && \
	cp ../../gs/libgs.a . && \
	cp ../../Veclib/libvec.a . && \
	cp ../../metis-4.0.3/libmetis.a .
	cd Nektar3d/STAMPEDE && make ARCH=STAMPEDE ADR=1 mopt
	cp Nektar3d/STAMPEDE/Nektar* .
	
clean:
	cd gs && rm -f *.a
	cd Veclib && rm -f *.a
	cd metis-4.0.3 && make clean
	cd Hlib/STAMPEDE && rm libgs.a libvec.a libmetis.a libhybridmopt.a *.tab.c
	cd Nektar3d/STAMPEDE && make clean
