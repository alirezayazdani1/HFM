/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $Source: 
 * $Revision:
 * $Date: 
 * $Author:
 * $State: 
 *---------------------------------------------------------------------------*/

#include <string.h>
#include "nektar.h"

double Flux(Domain *omega, char *label);
double Flux(Domain *omega, int *ID_faces_per_outlet, int Nfaces);
double FluxInlet(Domain *omega, char *label);

double Flux(Domain *omega, char *labeli, double *momentum_flux);
double FluxInlet(Domain *omega, char *label, double *momentum_flux);

double FluxLinmod(Domain *omega, char *label);
double PressureMean(Domain *omega, char *label);
double PressureMeanInlet(Domain *omega, char *label);
double BoundaryArea(Domain *omega, char *label);
double BoundaryArea(Bndry *Bc, char *label);
double BoundaryAreaInlet(Domain *omega, char *label);
double BoundaryAreaInlet(Bndry *Bc, char *label);


/* Flux - computes flow rate through boundary labeled as "label" */
double Flux(Domain *omega, char *label){

  int i,j,k,face,qa,qb;
  double *utmp,*u,*v,*w,*za,*zb,*wa,*wb;
  double sum,sum_i,vel;
  double FlowRate = 0.0;

  Bndry *Bc,*BcV,*BcW;
  Element *F;
  utmp = dvector(0, QGmax*QGmax-1);
  u = dvector(0, QGmax*QGmax-1);
  v = dvector(0, QGmax*QGmax-1);
  w = dvector(0, QGmax*QGmax-1);

  for(Bc=omega->Ubc,BcV=omega->Vbc,BcW=omega->Wbc;Bc;Bc=Bc->next,BcV=BcV->next,BcW=BcW->next){
    if ( (Bc->type == 'o' || Bc->type == 'O') && (strcmp(Bc->blabel,label) == 0) ){
      sum = 0.0;

      F = Bc->elmt;
      face = Bc->face;

      qa=Bc->elmt->qa;
      qb=Bc->elmt->qb;

      F = Bc->elmt;
//      F->Surface_geofac(Bc); 
//      F->Trans(F,J_to_Q); 
//      F->state = 't';
      F->GetFace(F->h_3d[0][0],face,utmp);
      F->InterpToFace1(face,utmp,u);

      F = BcV->elmt;
//      F->Surface_geofac(BcV); 
//      F->Trans(F,J_to_Q);
//      F->state = 't';
      F->GetFace(F->h_3d[0][0],face,utmp);
      F->InterpToFace1(face,utmp,v);

      F = BcW->elmt;
//      F->Surface_geofac(BcW); 
//      F->Trans(F,J_to_Q);
//      F->state = 't';
      F->GetFace(F->h_3d[0][0],face,utmp);
      F->InterpToFace1(face,utmp,w);

      i = qa*qb;

      if(F->curvX){
         dvmul  (i, Bc->nx.p, 1, u, 1, utmp, 1);
         dvvtvp (i, Bc->ny.p, 1, v, 1, utmp, 1, utmp, 1);
         dvvtvp (i, Bc->nz.p, 1, w, 1, utmp, 1, utmp, 1);
      }
      else{
         dsmul  (i, Bc->nx.d, u, 1, utmp, 1);
         dsvtvp (i, Bc->ny.d, v, 1, utmp, 1, utmp, 1);
         dsvtvp (i, Bc->nz.d, w, 1, utmp, 1, utmp, 1);
      }
      /* now utmp = u*nx+v*ny+w*nz                       */
      /* that is utmp = velocity normal to the surface   */

      if(F->curvX)
         dvmul(qa*qb,Bc->sjac.p,1,utmp,1,utmp,1);
      else
         dscal(qa*qb,Bc->sjac.d,utmp,1);

      getzw(qa,&za,&wa,'a');
      getzw(qb,&zb,&wb,'b');
#if 0
      k = 0;
      for (j = 0; j < qb; ++j){
        sum_i = 0.0;
        for (i = 0; i < qa; ++i){
          sum_i += wa[i]*utmp[k];
          k++; 
        }
        sum += (sum_i*wb[j]);
      }
#else
      for(j = 0; j < qb; ++j){
         for (i = 0; i < qa; ++i){
           k = i+j*qa;
           vel = utmp[k];
           sum += wb[j]*wa[i]*vel;
	 }
      }
#endif
      FlowRate += sum; 

    }  //end of "if (B->bstring == label)"
  }    // end of for 

  free(utmp); free(u); free(v); free(w);

  return FlowRate;
}

double Flux(Domain *omega, int *ID_faces_per_outlet, int Nfaces){

  int i,j,k,face,qa,qb,face_index;
  double *utmp,*u,*v,*w,*za,*wa,*wb;
  double sum,sum_i,vel;
  double FlowRate = 0.0;

  Bndry *Bc,*BcV,*BcW;
  Element *F;
  i = QGmax*QGmax;
  utmp = dvector(0, i-1);
  u = dvector(0, i-1);
  v = dvector(0, i-1);
  w = dvector(0, i-1);

  Bndry *Ubc = omega->Ubc;
  Bndry *Vbc = omega->Vbc;
  Bndry *Wbc = omega->Wbc;

  for (face_index = 0; face_index < Nfaces; face_index++){
    
    i = ID_faces_per_outlet[face_index];
    Bc  = Ubc+i;
    BcV = Vbc+i;
    BcW = Wbc+i;

    sum = 0.0;

    F = Bc->elmt;
    face = Bc->face;
    qa=F->qa;
    qb=F->qb;
    i = qa*qb;

      if (face == 0){
	  F->GetFace(F->h_3d[0][0],face,u);
	  F = BcV->elmt;
	  F->GetFace(F->h_3d[0][0],face,v);
	  F = BcW->elmt;
	  F->GetFace(F->h_3d[0][0],face,w);

	  if(F->curvX){
	      dvmul  (i, Bc->nx.p, 1, u, 1, utmp, 1);
	      dvvtvp (i, Bc->ny.p, 1, v, 1, utmp, 1, utmp, 1);
	      dvvtvp (i, Bc->nz.p, 1, w, 1, utmp, 1, utmp, 1);
	      dvmul(i,Bc->sjac.p,1,utmp,1,utmp,1);
	  }
	  else{
	      dsmul  (i, Bc->nx.d, u, 1, utmp, 1);
	      dsvtvp (i, Bc->ny.d, v, 1, utmp, 1, utmp, 1);
	      dsvtvp (i, Bc->nz.d, w, 1, utmp, 1, utmp, 1);
	      dscal(i,Bc->sjac.d,utmp,1);
	  }
      }
      else{       
	  if (F->curvX){

	      F->GetFace(F->h_3d[0][0],face,utmp);
	      F->InterpToFace1(face,utmp,u);

	      F = BcV->elmt;
	      F->GetFace(F->h_3d[0][0],face,utmp);
	      F->InterpToFace1(face,utmp,v);

	      F = BcW->elmt;
	      F->GetFace(F->h_3d[0][0],face,utmp);
	      F->InterpToFace1(face,utmp,w);

	      dvmul  (i, Bc->nx.p, 1, u, 1, utmp, 1);
	      dvvtvp (i, Bc->ny.p, 1, v, 1, utmp, 1, utmp, 1);
	      dvvtvp (i, Bc->nz.p, 1, w, 1, utmp, 1, utmp, 1);
	      dvmul(i,Bc->sjac.p,1,utmp,1,utmp,1);
	  }
	  else{

	      F->GetFace(F->h_3d[0][0],face,u);
	      dsmul  (i, Bc->nx.d, u, 1, utmp, 1);

	      F = BcV->elmt;
	      F->GetFace(F->h_3d[0][0],face,v);
	      dsvtvp (i, Bc->ny.d, v, 1, utmp, 1, utmp, 1);

	      F = BcW->elmt;
	      F->GetFace(F->h_3d[0][0],face,w);
	      dsvtvp (i, Bc->nz.d, w, 1, utmp, 1, u, 1);

	      F->InterpToFace1(face,u,utmp);
              dscal(i,Bc->sjac.d,utmp,1);
	  }
      }
      getzw(qa,&za,&wa,'a');
      getzw(qb,&za,&wb,'b');

      k = 0;
      for (j = 0; j < qb; ++j){
        sum_i = 0.0;
        for (i = 0; i < qa; ++i){
          sum_i += wa[i]*utmp[k];
          k++; 
        }
        sum += (sum_i*wb[j]);
      }

      FlowRate += sum; 
  }    // end of for 

  free(utmp); free(u); free(v); free(w); 

  return FlowRate;
}

/* Flux - computes flow rate through inlet boundary  */
double FluxInlet(Domain *omega, char *label){

  int i,j,k,face,qa,qb;
  double *utmp,*u,*v,*w,*za,*zb,*wa,*wb;
  double sum,vel;
  double FlowRate = 0.0;

  Bndry *Bc,*BcV,*BcW;
  Element *F;
  utmp = dvector(0, QGmax*QGmax-1);
  u    = dvector(0, QGmax*QGmax-1);
  v    = dvector(0, QGmax*QGmax-1);
  w    = dvector(0, QGmax*QGmax-1);

  for(Bc=omega->Ubc,BcV=omega->Vbc,BcW=omega->Wbc;Bc;Bc=Bc->next,BcV=BcV->next,BcW=BcW->next){
    if ( (Bc->type == 'v' || Bc->type == 'V') && strcmp(Bc->blabel,label) == 0 ){
      sum = 0.0;

      F = Bc->elmt;
      face = Bc->face;

      qa=Bc->elmt->qa;
      qb=Bc->elmt->qb;

      F = Bc->elmt;
      //F->Trans(F,J_to_Q);
      //F->state = 't';
      F->GetFace(F->h_3d[0][0],face,utmp);
      F->InterpToFace1(face,utmp,u);

      F = BcV->elmt;
      //F->Trans(F,J_to_Q);
      //F->state = 't';
      F->GetFace(F->h_3d[0][0],face,utmp);
      F->InterpToFace1(face,utmp,v);

      F = BcW->elmt;
      //F->Trans(F,J_to_Q);
      //F->state = 't';
      F->GetFace(F->h_3d[0][0],face,utmp);
      F->InterpToFace1(face,utmp,w);

      i = qa*qb;

      if(F->curvX){
         dvmul  (i, Bc->nx.p, 1, u, 1, utmp, 1);
         dvvtvp (i, Bc->ny.p, 1, v, 1, utmp, 1, utmp, 1);
         dvvtvp (i, Bc->nz.p, 1, w, 1, utmp, 1, utmp, 1);
      }
      else{
         dsmul  (i, Bc->nx.d, u, 1, utmp, 1);
         dsvtvp (i, Bc->ny.d, v, 1, utmp, 1, utmp, 1);
         dsvtvp (i, Bc->nz.d, w, 1, utmp, 1, utmp, 1);
      }
      if(F->curvX)
         dvmul(qa*qb,Bc->sjac.p,1,utmp,1,utmp,1);
      else
         dscal(qa*qb,Bc->sjac.d,utmp,1);

      getzw(qa,&za,&wa,'a');
      getzw(qb,&zb,&wb,'b');

      for(j = 0; j < qb; ++j){
         for (i = 0; i < qa; ++i){
           k = i+j*qa;
           vel = utmp[k];
           sum += wb[j]*wa[i]*vel;
	 }
      }
      FlowRate += sum; 
    }  //end of "if (B->bstring == label)"
  }    // end of for 

  free(utmp); free(u); free(v); free(w);
  return FlowRate;
}


/* Flux - computes approximate flow-rate through boundary labeled as "label" */
/* flow-rate  is computed basing on linear modes only                           */
double FluxLinmod(Domain *omega, char *label){

  int i,j,vn,nvert,k,face,qa,qb,qc;
  double Area,vel;
  double FlowRate = 0.0;
  double *w1,*w2,*w3,*za,*zb,*zc;

  Bndry *Bc,*BcV,*BcW;
  Element *E;
  
  for(Bc=omega->Ubc,BcV=omega->Vbc,BcW=omega->Wbc;Bc;Bc=Bc->next,BcV=BcV->next,BcW=BcW->next){
    if ( (Bc->type == 'o' || Bc->type == 'O') && (strcmp(Bc->blabel,label) == 0) ){
      
      E = Bc->elmt;
      face = Bc->face;
      nvert = E->Nfverts(face);
      qa=E->qa;
      qb=E->qb;
      qc=E->qc;
      getzw(qa,&za,&w1,'a'); /* first  direction */
      getzw(qb,&zb,&w2,'b'); /* second direction */
      getzw(qc,&zc,&w3,'c'); /* third  direction */

      if (E->curvX){

	    Area = 0.0;

	    switch(face){
		case 0:
 		    for(j = 0; j < qb; ++j){
			for (i = 0; i < qa; ++i){
			    k = i+j*qa;
			    Area += w2[j]*w1[i]*(Bc->sjac.p[k]);
			}
		    }
		case 1:
		    for(j = 0; j < qc; ++j){
			for (i = 0; i < qa; ++i){
			    k = i+j*qa;
			    Area += w3[j]*w1[i]*(Bc->sjac.p[k])*2.0/(1.0-zc[j]);
			}
		    }
		case 2:
		case 3:
		    for(j = 0; j < qc; ++j){
			for (i = 0; i < qb; ++i){
			    k = i+j*qb;
			    Area += w3[j]*w2[i]*(Bc->sjac.p[k])/(1.0-zb[i])/(1.0-zc[j])*4.0;
			}
		    }
	    }

//assume: face is flat  - normal to the face at any point of the face grid is the same

      vel = 0.0;
      for(j = 0; j < nvert; ++j){      
      	vn = E->vnum(face,j);
      	vel  += Bc->elmt->vert[vn].hj[0]*Bc->nx.p[0]+BcV->elmt->vert[vn].hj[0]*BcV->ny.p[0]+BcW->elmt->vert[vn].hj[0]*BcW->nz.p[0];
	    } 
	    FlowRate += vel*(2.0/3.0)*Area*0.5;
	}
	else{
	    vel = 0.0;
      for(j = 0; j < nvert; ++j){      
      	vn = E->vnum(face,j);
        vel  += Bc->elmt->vert[vn].hj[0]*Bc->nx.d+BcV->elmt->vert[vn].hj[0]*BcV->ny.d+BcW->elmt->vert[vn].hj[0]*BcW->nz.d;
	    } 
	    FlowRate += vel*(2.0/3.0)*Bc->sjac.d;
	}

    }
  }

  return FlowRate;
}

/* Flux - computes flow rate and momentum flux through boundary labeled as "label" */
double Flux(Domain *omega, char *label, double *momentum_flux){

  int i,j,k,face,qa,qb;
  double *utmp,*u,*v,*w,*za,*zb,*wa,*wb;
  double sum,tmp,vel;
  double sum_flux_u,sum_flux_v,sum_flux_w;
  double FlowRate = 0.0;

  Bndry *Bc,*BcV,*BcW;
  Element *F;
  utmp = dvector(0, QGmax*QGmax-1);
  u = dvector(0, QGmax*QGmax-1);
  v = dvector(0, QGmax*QGmax-1);
  w = dvector(0, QGmax*QGmax-1);

  momentum_flux[0] = 0.0; momentum_flux[1] = 0.0; momentum_flux[2] = 0.0;
  
  for(Bc=omega->Ubc,BcV=omega->Vbc,BcW=omega->Wbc;Bc;Bc=Bc->next,BcV=BcV->next,BcW=BcW->next){
    if ( (Bc->type == 'o' || Bc->type == 'O') && (strcmp(Bc->blabel,label) == 0) ){
      sum = 0.0;
      sum_flux_u = sum_flux_v = sum_flux_w = 0.0;

      F = Bc->elmt;
      face = Bc->face;
      qa=Bc->elmt->qa;
      qb=Bc->elmt->qb;

      F = Bc->elmt;
      F->Surface_geofac(Bc);
      F->Trans(F,J_to_Q);
      F->state = 't';
      F->GetFace(F->h_3d[0][0],face,utmp);
      F->InterpToFace1(face,utmp,u);

      F = BcV->elmt;
      F->Surface_geofac(BcV);
      F->Trans(F,J_to_Q);
      F->state = 't';
      F->GetFace(F->h_3d[0][0],face,utmp);
      F->InterpToFace1(face,utmp,v);

      F = BcW->elmt;
      F->Surface_geofac(BcW);
      F->Trans(F,J_to_Q);
      F->state = 't';
      F->GetFace(F->h_3d[0][0],face,utmp);
      F->InterpToFace1(face,utmp,w);

      i = qa*qb;

      if(F->curvX){
         dvmul  (i, Bc->nx.p, 1, u, 1, utmp, 1);
         dvvtvp (i, Bc->ny.p, 1, v, 1, utmp, 1, utmp, 1);
         dvvtvp (i, Bc->nz.p, 1, w, 1, utmp, 1, utmp, 1);
      }
      else{
         dsmul  (i, Bc->nx.d, u, 1, utmp, 1);
         dsvtvp (i, Bc->ny.d, v, 1, utmp, 1, utmp, 1);
         dsvtvp (i, Bc->nz.d, w, 1, utmp, 1, utmp, 1);
      }
      /* now utmp = u*nx+v*ny+w*nz                       */
      /* that is utmp = velocity normal to the surface   */

      if(F->curvX)
         dvmul(qa*qb,Bc->sjac.p,1,utmp,1,utmp,1);
      else
         dscal(qa*qb,Bc->sjac.d,utmp,1);

      getzw(qa,&za,&wa,'a');
      getzw(qb,&zb,&wb,'b');

      for(j = 0; j < qb; ++j){
         for (i = 0; i < qa; ++i){
           k = i+j*qa;
           vel = utmp[k];
           tmp = wb[j]*wa[i]*vel;
           sum += tmp;
           sum_flux_u += tmp*u[k];
           sum_flux_v += tmp*v[k];
           sum_flux_w += tmp*w[k];
         }
      }

      FlowRate += sum;
      momentum_flux[0] += sum_flux_u;
      momentum_flux[1] += sum_flux_v;
      momentum_flux[2] += sum_flux_w;

    }  //end of "if (B->bstring == label)"

  }    // end of for

  free(utmp); free(u); free(v); free(w);

  return FlowRate;
}

/* Flux - computes flow rate and momentum flux through inlet boundary  */
double FluxInlet(Domain *omega, char *label, double *momentum_flux){

  int i,j,k,face,qa,qb;
  double *utmp,*u,*v,*w,*za,*zb,*wa,*wb;
  double sum,tmp,vel;
  double sum_flux_u,sum_flux_v,sum_flux_w;
  double FlowRate = 0.0;

  Bndry *Bc,*BcV,*BcW;
  Element *F;
  utmp = dvector(0, QGmax*QGmax-1);
  u    = dvector(0, QGmax*QGmax-1);
  v    = dvector(0, QGmax*QGmax-1);
  w    = dvector(0, QGmax*QGmax-1);

  momentum_flux[0] = 0.0; momentum_flux[1] = 0.0; momentum_flux[2] = 0.0;

  for(Bc=omega->Ubc,BcV=omega->Vbc,BcW=omega->Wbc;Bc;Bc=Bc->next,BcV=BcV->next,BcW=BcW->next){
    if ( (Bc->type == 'v' || Bc->type == 'V') && strcmp(Bc->blabel,label) == 0 ){
      sum = 0.0;
      sum_flux_u = sum_flux_v = sum_flux_w = 0.0;

      F = Bc->elmt;
      face = Bc->face;
      qa=Bc->elmt->qa;
      qb=Bc->elmt->qb;

      F = Bc->elmt;
      F->Surface_geofac(Bc);
      F->Trans(F,J_to_Q);
      F->state = 't';
      F->GetFace(F->h_3d[0][0],face,utmp);
      F->InterpToFace1(face,utmp,u);

      F = BcV->elmt;
      F->Surface_geofac(BcV);
      F->Trans(F,J_to_Q);
      F->state = 't';
      F->GetFace(F->h_3d[0][0],face,utmp);
      F->InterpToFace1(face,utmp,v);

      F = BcW->elmt;
      F->Surface_geofac(BcW);
      F->Trans(F,J_to_Q);
      F->state = 't';
      F->GetFace(F->h_3d[0][0],face,utmp);
      F->InterpToFace1(face,utmp,w);

      i = qa*qb;

      if(F->curvX){
         dvmul  (i, Bc->nx.p, 1, u, 1, utmp, 1);
         dvvtvp (i, Bc->ny.p, 1, v, 1, utmp, 1, utmp, 1);
         dvvtvp (i, Bc->nz.p, 1, w, 1, utmp, 1, utmp, 1);
      }
      else{
         dsmul  (i, Bc->nx.d, u, 1, utmp, 1);
         dsvtvp (i, Bc->ny.d, v, 1, utmp, 1, utmp, 1);
         dsvtvp (i, Bc->nz.d, w, 1, utmp, 1, utmp, 1);
      }
      if(F->curvX)
         dvmul(qa*qb,Bc->sjac.p,1,utmp,1,utmp,1);
      else
         dscal(qa*qb,Bc->sjac.d,utmp,1);

      getzw(qa,&za,&wa,'a');
      getzw(qb,&zb,&wb,'b');

      for(j = 0; j < qb; ++j){
         for (i = 0; i < qa; ++i){
           k = i+j*qa;
           vel = utmp[k];
           tmp = wb[j]*wa[i]*vel;
           sum += tmp;
           sum_flux_u += tmp*u[k];
           sum_flux_v += tmp*v[k];
           sum_flux_w += tmp*w[k];
	 }
      }
      FlowRate += sum; 
      momentum_flux[0] += sum_flux_u;
      momentum_flux[1] += sum_flux_v;
      momentum_flux[2] += sum_flux_w;

    }  //end of "if (B->bstring == label)"

  }    // end of for 

 
  free(utmp); free(u); free(v); free(w);

  return FlowRate;
}


/* PressureMean - integrates pressure over boundary with a type = 'o' or 'O' and Bc->blabel = label as    */
double PressureMean(Domain *omega, char *label){
                                                                                                                                                             
  int i,j,k,face,qa,qb;
  double *p,*ptmp,*za,*zb,*wa,*wb;
  double sum,vel;
  double Pressure_mean = 0.0;

  Bndry *Bc, *BcP;
  Element *F;
  p    = dvector(0, QGmax*QGmax-1);
  ptmp = dvector(0, QGmax*QGmax-1);                                                                                                                                          

  for(BcP=omega->Pbc;BcP;BcP=BcP->next){
    if ( (BcP->usrtype == 'O' || BcP->usrtype == 'o') && strcmp(BcP->blabel,label) == 0 ){
                                                                                     
      sum = 0.0;
      F = BcP->elmt;
      face = BcP->face;
      qa=F->qa;
      qb=F->qb;

      F->Surface_geofac(BcP);
      F->Trans(F,J_to_Q);
      F->state = 't';
      F->GetFace(F->h_3d[0][0],face,ptmp);
      F->InterpToFace1(face,ptmp,p);

      if(F->curvX)
         dvmul(qa*qb,BcP->sjac.p,1,p,1,p,1);
      else
         dscal(qa*qb,BcP->sjac.d,p,1);

      getzw(qa,&za,&wa,'a');
      getzw(qb,&zb,&wb,'b');

      for(j = 0; j < qb; ++j){
         for (i = 0; i < qa; ++i){
           k = i+j*qa;
           vel = p[k];
           sum += wb[j]*wa[i]*vel;
	 }
      }          

      Pressure_mean += sum;                                                                                                                                         
    }  //end of "if (B->type ==  'o')"
                                                                                                                                                             
  }    // end of for                                                                                                                                                                                                                                  
  free(ptmp); free(p);

  return Pressure_mean;                                                                                            
}

/* PressureMeanInlet - integrates pressure over boundary with a type = 'v' or 'V'  */
double PressureMeanInlet(Domain *omega,  char *label){
                                                                                           
  int i,j,k,face,qa,qb;
  double *ptmp,*p,*za,*zb,*wa,*wb;
  double sum,vel;
  double Pressure_mean = 0.0;

  Bndry *Bc, *BcP;
  Element *F;
  p    = dvector(0, QGmax*QGmax-1);
  ptmp = dvector(0, QGmax*QGmax-1);

  for(BcP=omega->Pbc;BcP;BcP=BcP->next){
    if ( (BcP->usrtype == 'v' || BcP->usrtype == 'V') && strcmp(BcP->blabel,label) == 0 ){

      sum = 0.0;
      F = BcP->elmt;
      face = BcP->face;
      qa=F->qa;
      qb=F->qb;

      F->Surface_geofac(BcP);
      F->Trans(F,J_to_Q);
      F->state = 't';
      F->GetFace(F->h_3d[0][0],face,ptmp);
      F->InterpToFace1(face,ptmp,p);

      if(F->curvX)
         dvmul(qa*qb,BcP->sjac.p,1,p,1,p,1);
      else
         dscal(qa*qb,BcP->sjac.d,p,1);

      getzw(qa,&za,&wa,'a');
      getzw(qb,&zb,&wb,'b');

      for(j = 0; j < qb; ++j){
         for (i = 0; i < qa; ++i){
           k = i+j*qa;
           vel = p[k];
           sum += wb[j]*wa[i]*vel;
	 }
      }                                 

      Pressure_mean += sum;                                                                                                                                         
    }  //end of "if (B->type == 'v')"
                                                                                                                                                             
  }    // end of for                                                                                                                                                                                                                                  
  free(ptmp); free(p);
                                                                                                                                                  
  return Pressure_mean;                                                                                            
}


/* BoundaryArea - computes the area of boundary with a type = 'o' or 'O' and Bc->blabel = label   */
double BoundaryArea(Domain *omega, char *label){
  return  BoundaryArea(omega->Ubc,label);
}


double BoundaryArea(Bndry *Ubc, char *label){
                                                                                                                                                          
  int i,j,k,qa,qb;
  double *za,*zb,*wa,*wb;
  double Area = 0.0;

  Bndry *Bc;
  Element *F;
                                                                 
  for(Bc=Ubc;Bc;Bc=Bc->next){
    if ( (Bc->type == 'O' || Bc->type == 'o') && strcmp(Bc->blabel,label) == 0 ){
                                                                                    
      F = Bc->elmt;  
      qa=F->qa;
      qb=F->qb;
                                                                
      F->Surface_geofac(Bc);
            
      getzw(qa,&za,&wa,'a');
      getzw(qb,&zb,&wb,'b');

      if(F->curvX){
	  for(j = 0; j < qb; ++j){
	      for (i = 0; i < qa; ++i){
		  k = i+j*qa;
		  Area += wb[j]*wa[i]*Bc->sjac.p[k];
	      }
	  }
      }
      else{
         for(j = 0; j < qb; ++j){
            for (i = 0; i < qa; ++i){
              k = i+j*qa;
              Area += wb[j]*wa[i]*Bc->sjac.d;
	    }
	 }         
      }
    }  //end of "if (B->bstring == label)"
                                                                                                                                                             
  }    // end of for                
                                                                                                                                                  
  return Area;  
}


/* BoundaryAreaInlet - computes the area of boundary with a type = 'v' or 'V'  */
double BoundaryAreaInlet(Domain *omega, char *label){
  return BoundaryAreaInlet(omega->Ubc,label);
}

double BoundaryAreaInlet(Bndry *Ubc, char *label){
                                                                                                                                                          
  int i,j,k,qa,qb;
  double *za,*zb,*wa,*wb;
  double  Area = 0.0;

  Bndry *Bc;
  Element *F;
  
  for(Bc=Ubc;Bc;Bc=Bc->next){
    if ( (Bc->type == 'v' || Bc->type == 'V') && strcmp(Bc->blabel,label) == 0 ){
      
      F = Bc->elmt;  
      qa=F->qa;
      qb=F->qb;
                                                                
      F->Surface_geofac(Bc);

      getzw(qa,&za,&wa,'a');
      getzw(qb,&zb,&wb,'b');

      if(F->curvX){
	  for(j = 0; j < qb; ++j){
	      for (i = 0; i < qa; ++i){
		  k = i+j*qa;
		  Area += wb[j]*wa[i]*Bc->sjac.p[k];
	      }
	  }
      }
      else{
         for(j = 0; j < qb; ++j){
            for (i = 0; i < qa; ++i){
              k = i+j*qa;
              Area += wb[j]*wa[i]*Bc->sjac.d;
	    }
	 }         
      }                                                                                        
    }  //end of "if (B->bstring == label)"
                                                                                                                                                             
  }    // end of for                
                                                                                                                                                  
  return Area;  
}
