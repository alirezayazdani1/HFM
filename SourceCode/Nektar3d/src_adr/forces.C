/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $RCSfile:
 * $Revision:
 * $Author: 
 * $Date: 
 * $State:
 * ------------------------------------------------------------------------- */
#include "nektar.h"

static int init=1; /* externals */

static void set_SurGeofac(Bndry *Pbc, Bndry *Ubc);

/*-------------------------------------------------------------------*
 * This is a function to calculate the forces of the fluid acting on *
 * a body pressumed to be represented by the flag 'W'.               *
 *                                                                   *
 *   F_i = P n_i - RHO*KINVIS*T_ij n_j  (RHO=1)                      *
 *                                                                   *
 *-------------------------------------------------------------------*/

static int surf_number=0;

void forces(Domain *omega, int step, double time){
  register int i;
  int      eid, qa, qb, face;
  int  surfid, flag;
  Bndry   *B = omega->Ubc;
  double *za, *wa, *zb, *wb, *zc, *wc;
  double  fp[3],fv[3],**Fv;
  double  kinvis = dparam("KINVIS");
  double  *w  = dvector(0,QGmax*QGmax-1);
  double  *wk = dvector(0,3*QGmax*QGmax*QGmax-1);
  double  **D = dmatrix(0,10,0,QGmax*QGmax*QGmax-1);
  double  *ux = D[0], *uy = D[1], *uz = D[2];
  double  *vx = D[3], *vy = D[4], *vz = D[5];
  double  *wx = D[6], *wy = D[7], *wz = D[8], *p = D[9];
  Coord   X;
  Element_List *U  = omega->U, *V = omega->V, *W = omega->W, *P = omega->P;
	// to-do: set-up the wall flux calculations correctly
#ifdef ADR
//	Element_List *T = omega->T;
//	Element *eT;
//	double *t = D[10];
//	double  qy[0];
#endif
  static Element_List *Surf[4];
  FILE    *fout = omega->fce_file;
  Element *eU, *eV, *eW, *eP;
  fp[0] = fp[1] = fp[2] = 0.0;
  fv[0] = fv[1] = fv[2] = 0.0;

  static int srf=option("SurForce");

  Fv = dmatrix(0,2,0,QGmax*QGmax-1);
  
  /* print header */
  if(init){
    ROOTONLY{
      fprintf(fout,"# Force acting on body\n");
      fprintf(fout,"# \n");
#ifdef ADR
      fprintf(fout,"# Time  (Fx-press, Fx-visc) Fx  (Fy-press, Fy-visc)"
	      "Fy  (Fz-press, Fz-visc) Fz\n");
#else
      fprintf(fout,"# Time  (Fx-press, Fx-visc) Fx  (Fy-press, Fy-visc)"
        "Fy  (Fz-press, Fz-visc) Fz\n");
#endif
    }
    /* need to set up surface geometric factors in U from P */
    set_SurGeofac(omega->Pbc,omega->Ubc);

    if(srf&&!Surf[0]){ // setup surface structure
      Surf[0] = setup_surflist(omega,B,'W');

      if(Surf[0]){
	for(i = 0; i < Surf[0]->nel; ++i)
	  Surf[0]->flist[i]->type = 'p';

	Surf[1] = Surf[0]->gen_aux_field('s');
	Surf[2] = Surf[0]->gen_aux_field('r');
	Surf[3] = Surf[0]->gen_aux_field('t');
	
	//at this point Surf[] is still in physical space
      }
    }   
    init = 0;
  }

  surfid=0;
  
  for(B = omega->Ubc;B; B = B->next)
    if(B->type == 'W' || B->usrtype == 'm' || B->usrtype == 'M'){
      face = B->face;
      eid  = B->elmt->id;

     
      eU = U->flist[eid];
      eV = V->flist[eid];
      eW = W->flist[eid];
      eP = P->flist[eid];
#ifdef ADR
//      eT = T->flist[eid];
#endif

      // get zero weights and order of face quad 
      if(eU->identify() == Nek_Tet || eU->identify() == Nek_Hex){
        qa   = eU->qa;
        qb   = eU->qb;
        eU->GetZW(&za, &wa, &zb, &wb, &zc, &wc);
      }
      else if(eU->identify() == Nek_Prism){
        if(face == 1 || face == 3){
          qa = eU->qa;
          qb = eU->qc;
          eU->GetZW(&za, &wa, &zc , &wc, &zb, &wb);
        }
        else{
          qa   = eU->qa;
          qb   = eU->qb;
          eU->GetZW(&za, &wa, &zb, &wb, &zc, &wc);
        }
      }
      else{
	fprintf(stderr,"forces.c nead setting up for pyramids\n");
	exit(1);
      }

      X.x = dvector(0,qa*qb-1);
      X.y = dvector(0,qa*qb-1);
      X.z = dvector(0,qa*qb-1);
      eU->GetFaceCoord(B->face,&X);

      for(i = 0; i < qb; ++i)
	dcopy(qa,wa,1,w+i*qa,1);
      for(i = 0; i < qb; ++i)
	dsmul(qa,wb[i],w+i*qa,1,w+i*qa,1);
      
      if(eU->curvX)
        dvmul(qa*qb,B->sjac.p,1,w,1,w,1);
      else
        dsmul(qa*qb,B->sjac.d,w,1,w,1); 

      // put field into physical space
      eP->Trans(eP, J_to_Q); eP->state = 't';
      eU->Trans(eU, J_to_Q); eU->state = 't';
      eV->Trans(eV, J_to_Q); eV->state = 't';
      eW->Trans(eW, J_to_Q); eW->state = 't';
#ifdef ADR
//      eT->Trans(eT, J_to_Q); eT->state = 't';
#endif        
      
      eU->Grad_d(ux, uy, uz, 'a');
      eV->Grad_d(vx, vy, vz, 'a');
      eW->Grad_d(wx, wy, wz, 'a');
 

      /* get appropriate face values from D */
      for(i = 0; i < 9; ++i){
	eU->GetFace(D[i],face,wk);
	eU->InterpToFace1(face,wk,D[i]);
      }

      eP->GetFace(**eP->h_3d,face,wk);
      eP->InterpToFace1(face,wk,p);
#ifdef ADR
//      eT->GetFace(**eT->h_3d,face,wk);
//      eT->InterpToFace1(face,wk,t);
#endif
   
      flag = 1;
      for(i = 0; i < qa*qb; ++i){
	if ( X.x[i]<1.0 || abs(X.y[i]>1e-5) ){
	   flag = 0;
	   };

      };
      if ( flag){
      if(eU->curvX){
        for(i = 0; i < qa*qb; ++i){
          fp[0] += p[i]*B->nx.p[i]*w[i];
          fp[1] += p[i]*B->ny.p[i]*w[i];
          fp[2] += p[i]*B->nz.p[i]*w[i];
            
#ifdef ADR
//          qy[0]    += t[i]*B->ny.p[i]*w[i];
#endif

          Fv[0][i] = kinvis*(2.0*ux[i]*B->nx.p[i] + (uy[i] + vx[i])*B->ny.p[i]
                           + (uz[i] + wx[i])*B->nz.p[i]);
          Fv[1][i] = kinvis*(2.0*vy[i]*B->ny.p[i] + (uy[i] + vx[i])*B->nx.p[i]
                           + (vz[i] + wy[i])*B->nz.p[i]);
          Fv[2][i] = kinvis*(2.0*wz[i]*B->nz.p[i] + (uz[i] + wx[i])*B->nx.p[i]
                           + (vz[i] + wy[i])*B->ny.p[i]);
        }
      }
      else
        for(i = 0; i < qa*qb; ++i){
          fp[0] += p[i]*B->nx.d*w[i];
          fp[1] += p[i]*B->ny.d*w[i];
          fp[2] += p[i]*B->nz.d*w[i];
#ifdef ADR            
//          qy[0]    += t[i]*B->ny.d*w[i];
#endif

          Fv[0][i] = kinvis*(2.0*ux[i]*B->nx.d + (uy[i] + vx[i])*B->ny.d
                     + (uz[i] + wx[i])*B->nz.d);
          Fv[1][i] = kinvis*(2.0*vy[i]*B->ny.d + (uy[i] + vx[i])*B->nx.d
                     + (vz[i] + wy[i])*B->nz.d);
          Fv[2][i]= kinvis*(2.0*wz[i]*B->nz.d + (uz[i] + wx[i])*B->nx.d
			    + (vz[i] + wy[i])*B->ny.d);
        }

      if(srf){
        dcopy(qa*qb,p,1,Surf[0]->flist[surfid]->h[0],1);
        dcopy(qa*qb,Fv[0],1,Surf[1]->flist[surfid]->h[0],1);
        dcopy(qa*qb,Fv[1],1,Surf[2]->flist[surfid]->h[0],1);
        dcopy(qa*qb,Fv[2],1,Surf[3]->flist[surfid]->h[0],1);
      }
      
      dvmul(qa*qb,w,1,Fv[0],1,Fv[0],1);
      fv[0] -= dsum(qa*qb,Fv[0],1);
      dvmul(qa*qb,w,1,Fv[1],1,Fv[1],1);
      fv[1] -= dsum(qa*qb,Fv[1],1);
      dvmul(qa*qb,w,1,Fv[2],1,Fv[2],1);
      fv[2] -= dsum(qa*qb,Fv[2],1);

      free(X.x);free(X.y);free(X.z);
      surfid++;     
    }
};
  DO_PARALLEL{
    gdsum(fp,3,wk);
    gdsum(fv,3,wk);
#ifdef ADR
//    gdsum(qy,1,wk);	
#endif
  }

  ROOTONLY
#ifdef ADR
    fprintf(fout,"%lf ( %lf %lf ) %lf ( %lf %lf ) %lf ( %lf %lf ) %lf \n",
	    time,fp[0],fv[0],fp[0]+fv[0],fp[1],fv[1],fp[1]+fv[1],fp[2],
	    fv[2],fp[2]+fv[2]);
#else
    fprintf(fout,"%lf ( %lf %lf ) %lf ( %lf %lf ) %lf ( %lf %lf ) %lf \n",
      time,fp[0],fv[0],fp[0]+fv[0],fp[1],fv[1],fp[1]+fv[1],fp[2],
      fv[2],fp[2]+fv[2]);
#endif


  // Dump data at every hisstep. 
  if(srf && Surf[0]){/* dump surface data */
    
    if ((step % option("hisstep")) == 0){
      char buf[BUFSIZ];
      FILE *fp[2];
      
      for(i = 0; i < 4; ++i){
	Surf[i]->Trans(Surf[i], Q_to_J);
	Surf[i]->Set_state('t');
      }

      DO_PARALLEL{
        sprintf(buf,"%s_surf_%d.fld.hdr.%d",omega->name,surf_number,
		pllinfo.procid);
        fp[0] = fopen(buf,"w");
        sprintf(buf,"%s_surf_%d.fld.dat.%d",omega->name,surf_number,
		pllinfo.procid);
        fp[1] = fopen(buf,"w");

        option_set("Sur_dump",1);
        Writefld(fp, omega->name, step, time, 4, Surf);

        if(dparam("concat")){
#ifdef PARALLEL
	  gsync();
#endif
	  
	  ROOTONLY{
	    int nfiles;
	    nfiles = pllinfo.nprocs; 
	    char filelist[BUFSIZ],syscall[BUFSIZ];
	    
	    /*system call to cat the .dat files */
	    sprintf(filelist,"%s_surf_%d.fld.dat.0",omega->name,surf_number);
	    
	    for(i = 1; i < nfiles; ++i){
	      sprintf(filelist,"%s %s_surf_%d.fld.dat.%d ",filelist, 
		      omega->name, surf_number, i);
	    }

	    // concatinate files
	    sprintf(syscall,"/bin/cat %s > %s_surf_%d.fld.tot.dat",filelist,
		    omega->name, surf_number);
	    fprintf(stdout,"system call 1 : %d \n",system(syscall));
	    
	    // remove files
	    sprintf(syscall,"/bin/rm -rf %s ",filelist);
	    fprintf(stdout,"system call 2 : %d \n",system(syscall));
	  }
	  else if(pllinfo.procid == 1){
	    int nfiles;
	    nfiles = pllinfo.nprocs; 
	    char filelist[BUFSIZ],syscall[BUFSIZ];
	    
	    /*system call to cat the .dat files */
	    sprintf(filelist,"%s_surf_%d.fld.hdr.0",omega->name,surf_number);
	    
	    for(i = 1; i < nfiles; ++i)
	      sprintf(filelist,"%s %s_surf_%d.fld.hdr.%d ",filelist, 
		      omega->name, surf_number, i);
	    
	    // concatinate files
	    sprintf(syscall,"/bin/cat %s > %s_surf_%d.fld.tot.hdr",filelist,
		    omega->name, surf_number);
	    system(syscall);
	    
	    // remove files
	    sprintf(syscall,"/bin/rm -rf %s ",filelist);
	    system(syscall);
	  }
	}

        option_set("Sur_dump",0);

        fclose(fp[0]);
        fclose(fp[1]);
      }
      else{
        sprintf(buf,"%s_surf_%d.fld",omega->name,surf_number);
        fp[0] = fp[1] = fopen(buf,"w");

        option_set("Sur_dump",1);
        Writefld(fp, omega->name, step, time, 4, Surf);
        option_set("Sur_dump",0);

        fclose(fp[0]);
      }
    }
  }
  else
#ifdef PARALLEL
    gsync(); // place to syncronize the dump call
#endif
  
  ++surf_number;
  
  free(w); free(wk); free_dmatrix(D,0,0); free_dmatrix(Fv,0,0);
}


static void set_SurGeofac(Bndry *Pbc, Bndry *Ubc){
  Bndry *Ebc;
  
  for(;Pbc; Pbc = Pbc->next)
    if(Pbc->type == 'F')
      for(Ebc = Ubc; Ebc; Ebc = Ebc->next)
	if((Ebc->elmt->id == Pbc->elmt->id)&&(Ebc->face == Pbc->face)){
	  Ebc->sjac = Pbc->sjac;
	  Ebc->nx   = Pbc->nx;
	  Ebc->ny   = Pbc->ny;
	  Ebc->nz   = Pbc->nz;
	} 
}

