/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $Source: 
 * $Revision:
 * $Date: 
 * $Author:
 * $State: 
 *---------------------------------------------------------------------------*/

#include <stdio.h>
#include <time.h>
#include "nektar.h"

static double *escales = (double*)0;

void   cfl_setup(Element_List *U, double dt){
  Element      *E;
  register int i;
  int      id0,id1;
  double   *w, *z;
  double   medge,elen;
  
  escales = dvector(0, U->nel-1);
  
  dzero(U->nel, escales, 1);
  
  medge = 1e10;
  for(E=U->fhead;E;E=E->next){
    for(i = 0; i < E->Nedges; ++i){
      id0 = E->edvnum(i,0);
      id1 = E->edvnum(i,1);
      elen = (E->vert[id1].x-E->vert[id0].x)*(E->vert[id1].x-E->vert[id0].x) + 
	(E->vert[id1].x-E->vert[id0].x)*(E->vert[id1].x-E->vert[id0].x) + 
	(E->vert[id1].y-E->vert[id0].y)*(E->vert[id1].y-E->vert[id0].y) + 
	(E->vert[id1].y-E->vert[id0].y)*(E->vert[id1].y-E->vert[id0].y) + 
	(E->vert[id1].z-E->vert[id0].z)*(E->vert[id1].z-E->vert[id0].z) + 
	(E->vert[id1].z-E->vert[id0].z)*(E->vert[id1].z-E->vert[id0].z);
      medge = min(elen,medge);
    }
    
    getzw(E->qa, &z, &w, 'a');
    
    escales[E->id] = dt/(0.5*(z[1]-z[0])*sqrt(medge));
  }
}

double cfl_checker(Domain *omega, double dt){

  Element   *U, *V, *W;
  double    cfl, umax;
  double    *tmp = dvector(0, QGmax*QGmax*QGmax-1);

  //if(!escales)
    cfl_setup(omega->U,dt);
  
  cfl = 0.;

  for(U=omega->U->fhead,V=omega->V->fhead,W=omega->W->fhead;U;
      U=U->next,V=V->next,W=W->next){
    
    dvmul (U->qtot, U->h_3d[0][0], 1, U->h_3d[0][0], 1, tmp, 1);
    dvvtvp(V->qtot, V->h_3d[0][0], 1, V->h_3d[0][0], 1, tmp, 1, tmp, 1);
    dvvtvp(W->qtot, W->h_3d[0][0], 1, W->h_3d[0][0], 1, tmp, 1, tmp, 1);

    umax = sqrt(tmp[idmax(U->qtot, tmp, 1)])*escales[U->id];
    
    cfl = max(cfl,umax);
  }
  free(tmp);
  free(escales);
  DO_PARALLEL
    gdmax(&cfl,1,&umax);
    
  return cfl;
}

    
double full_cfl_checker(Domain *omega, double dt, int *eid_max){
  Element   *U, *V, *W;
  double    cfl, umax;
  double    *uloc = dvector(0, 3*QGmax*QGmax*QGmax-1), *vloc, *wloc;
  int       qtot,id;
  Geom      *G;

  vloc = uloc + QGmax*QGmax*QGmax;
  wloc = vloc + QGmax*QGmax*QGmax;

  cfl = 0.;

  for(U=omega->U->fhead,V=omega->V->fhead,W=omega->W->fhead;U;
      U=U->next,V=V->next,W=W->next){

    qtot = U->qtot;
    G    = U->geom;

    // calculate local velocities
    if(U->curvX){
      dvmul  (qtot, G->rx.p, 1, U->h_3d[0][0], 1, uloc, 1);
      dvvtvp (qtot, G->ry.p, 1, V->h_3d[0][0], 1, uloc, 1, uloc, 1);
      dvvtvp (qtot, G->rz.p, 1, W->h_3d[0][0], 1, uloc, 1, uloc, 1);

      dvmul  (qtot, G->sx.p, 1, U->h_3d[0][0], 1, vloc, 1);
      dvvtvp (qtot, G->sy.p, 1, V->h_3d[0][0], 1, vloc, 1, vloc, 1);
      dvvtvp (qtot, G->sz.p, 1, W->h_3d[0][0], 1, vloc, 1, vloc, 1);

      dvmul  (qtot, G->tx.p, 1, U->h_3d[0][0], 1, wloc, 1);
      dvvtvp (qtot, G->ty.p, 1, V->h_3d[0][0], 1, wloc, 1, wloc, 1);
      dvvtvp (qtot, G->tz.p, 1, W->h_3d[0][0], 1, wloc, 1, wloc, 1);
    }
    else{
      dsmul  (qtot,G->rx.d, U->h_3d[0][0], 1, uloc, 1);
      dsvtvp (qtot,G->ry.d, V->h_3d[0][0], 1, uloc, 1, uloc, 1);
      dsvtvp (qtot,G->rz.d, W->h_3d[0][0], 1, uloc, 1, uloc, 1);

      dsmul  (qtot,G->sx.d, U->h_3d[0][0], 1, vloc, 1);
      dsvtvp (qtot,G->sy.d, V->h_3d[0][0], 1, vloc, 1, vloc, 1);
      dsvtvp (qtot,G->sz.d, W->h_3d[0][0], 1, vloc, 1, vloc, 1);

      dsmul  (qtot,G->tx.d, U->h_3d[0][0], 1, wloc, 1);
      dsvtvp (qtot,G->ty.d, V->h_3d[0][0], 1, wloc, 1, wloc, 1);
      dsvtvp (qtot,G->tz.d, W->h_3d[0][0], 1, wloc, 1, wloc, 1);
    }

    dvmul (U->qtot, uloc, 1, uloc, 1, uloc, 1);
    dvvtvp(V->qtot, vloc, 1, vloc, 1, uloc, 1, uloc, 1);
    dvvtvp(W->qtot, wloc, 1, wloc, 1, uloc, 1, uloc, 1);
    
    id = idmax(U->qtot,uloc,1);
    umax = sqrt(uloc[id])*0.2*(U->lmax-1)*(U->lmax-1)*dt;
    
    if(umax > cfl){
      eid_max[0] = U->id;
      cfl = umax;
    }
  }
  
  free(uloc);
  
  DO_PARALLEL{
    gdmax(&cfl,1,&umax);
    gimax(eid_max,1,&id);
  }

  return cfl;
}


