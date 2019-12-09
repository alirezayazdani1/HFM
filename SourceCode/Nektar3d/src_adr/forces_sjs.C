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

void forces(Domain *omega, double time){
  register int i;
  int      eid, qa, qb, face, lmax=0;
  Bndry   *B = omega->Ubc;

  double *za, *wa, *zb, *wb, *zc, *wc;
  
  double  fp[3],fv[3];
  double  kinvis = dparam("KINVIS");
  double  *w  = dvector(0,QGmax*QGmax-1);
  double  *wk = dvector(0,3*QGmax*QGmax*QGmax-1);
  double  **D = dmatrix(0,9,0,QGmax*QGmax*QGmax-1);
  double  *ux = D[0], *uy = D[1], *uz = D[2];
  double  *vx = D[3], *vy = D[4], *vz = D[5];
  double  *wx = D[6], *wy = D[7], *wz = D[8], *p = D[9];
  Element_List *U  = omega->U, *V = omega->V, *W = omega->W, *P = omega->P;
  FILE    *fout = omega->fce_file;
  Basis   *b;
  Element *eU, *eV, *eW, *eP;
  fp[0] = fp[1] = fp[2] = 0.0;
  fv[0] = fv[1] = fv[2] = 0.0;
  
  /* print header */
  if(init){
    ROOTONLY{
      fprintf(fout,"# Force acting on body\n");
      fprintf(fout,"# \n");
      fprintf(fout,"# Time  (Fx-press, Fx-visc) Fx  (Fy-press, Fy-visc)"
	      "Fy  (Fz-press, Fz-visc) Fz \n");
    }
    /* need to set up surface geometric factors in U from P */
    set_SurGeofac(omega->Pbc,omega->Ubc);
    
    init = 0;
  }

  for(B = omega->Ubc;B; B = B->next)
    if(B->type == 'W' || B->usrtype == 'm' || B->usrtype == 'M'){
      face = B->face;
      eid  = B->elmt->id;
      
      eU = U->flist[eid];
      eV = V->flist[eid];
      eW = W->flist[eid];
      eP = P->flist[eid];

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
      b    = eU->getbasis();
      lmax = eU->lmax;

      for(i = 0; i < qb; ++i)
	dcopy(qa,wa,1,w+i*qa,1);
      for(i = 0; i < qb; ++i)
	dsmul(qa,wb[i],w+i*qa,1,w+i*qa,1);
      
  
      eP->Trans(eP, J_to_Q); eP->state = 't';
      eU->Trans(eU, J_to_Q); eU->state = 't';
      eV->Trans(eV, J_to_Q); eV->state = 't';
      eW->Trans(eW, J_to_Q); eW->state = 't';
      
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

      if(eU->curvX)
	for(i = 0; i < qa*qb; ++i){
	  fp[0] += p[i]*B->nx.p[i]*B->sjac.p[i]*w[i];
	  fp[1] += p[i]*B->ny.p[i]*B->sjac.p[i]*w[i];
	  fp[2] += p[i]*B->nz.p[i]*B->sjac.p[i]*w[i];

	  fv[0] -= kinvis*(2.0*ux[i]*B->nx.p[i] + (uy[i] + vx[i])*B->ny.p[i]
			   + (uz[i] + wx[i])*B->nz.p[i])*B->sjac.p[i]*w[i];
	  fv[1] -= kinvis*(2.0*vy[i]*B->ny.p[i] + (uy[i] + vx[i])*B->nx.p[i]
			   + (vz[i] + wy[i])*B->nz.p[i])*B->sjac.p[i]*w[i];
	  fv[2] -= kinvis*(2.0*wz[i]*B->nz.p[i] + (uz[i] + wx[i])*B->nx.p[i]
		     + (vz[i] + wy[i])*B->ny.p[i])*B->sjac.p[i]*w[i];
	}
      else
	for(i = 0; i < qa*qb; ++i){
	  fp[0] += p[i]*B->nx.d*B->sjac.d*w[i];
	  fp[1] += p[i]*B->ny.d*B->sjac.d*w[i];
	  fp[2] += p[i]*B->nz.d*B->sjac.d*w[i];

	  fv[0] -= kinvis*(2.0*ux[i]*B->nx.d + (uy[i] + vx[i])*B->ny.d
		     + (uz[i] + wx[i])*B->nz.d)*B->sjac.d*w[i];
	  fv[1] -= kinvis*(2.0*vy[i]*B->ny.d + (uy[i] + vx[i])*B->nx.d
		     + (vz[i] + wy[i])*B->nz.d)*B->sjac.d*w[i];
	  fv[2] -= kinvis*(2.0*wz[i]*B->nz.d + (uz[i] + wx[i])*B->nx.d
		     + (vz[i] + wy[i])*B->ny.d)*B->sjac.d*w[i];
	}
    }

  DO_PARALLEL{
    gdsum(fp,3,wk);
    gdsum(fv,3,wk);
  }
  
  ROOTONLY
    fprintf(fout,"%lf ( %lf %lf ) %lf ( %lf %lf ) %lf ( %lf %lf ) %lf \n",
	    time,fp[0],fv[0],fp[0]+fv[0],fp[1],fv[1],fp[1]+fv[1],fp[2],
	    fv[2],fp[2]+fv[2]);

  free(w); free(wk); free_dmatrix(D,0,0);
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
