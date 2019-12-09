#include <stdio.h>
#include <math.h>
#include <veclib.h>
#include <hotel.h>
#include <time.h>
#include "nektar.h"

#define EulerStep(htot,Un,Uf,U,dt) \
 dsvtvp(htot,dt,Uf,1,Un,1,U,1); 

static void Upwind_edge_advection(Element_List *U, Element_List *V, Element_List *W, Element_List *UF, Element_List *VF, Element_List *WF);
static void Zero_Symm_BCs(Bndry *Ubc, Element_List *U, Element_List *V, 
			  Element_List *W,Element_List *Uf, Element_List *Vf, Element_List *Wf);
static void Fill_sides(Element_List *U, Element_List *Uf, Bndry *Ubc, 
		       Domain *Omega, double t);
static void Fill_BC_sides(Element_List *Uf, Bndry *Ubc, 
			  Domain *Omega, double t);
static void Project   (Element_List *U, Element_List *Uf, double *wk);
static void SubtractBC(Element_List *in, Element_List *out, Bndry *Bc, double dt);
static void ExtField(int htot, double t, double *un, double **up, double dt, 
		     int ord);

// RK arrays
static double BarrayC[2][2]    = {{1,0},{0.5,1}};
static double BarrayA[2][2][2] = {{{1,0},{0,0}},{{0.5,0},{-1,2}}};
static double BarrayB[2][3]    = {{0.5,0.5,0},{1/6.0,2/3.0,1.0/6.0}};

static void AverageVel(int N, double  **Uf,double *Barray, double  *Uout, 
		       int ntot);

#define TIMING

#ifdef TIMING
#define AddTime(s) \
ROOTONLY{ s += (clock()-st)/cps; \
 st = clock(); }
#else
#define AddTime(s) \
 /* Nothing */
#endif

static void Advance (int N, double dt, double *u, double *v, double *w,
		     Domain *Omega, int Je);

static double st, ts[7], cps = (double)CLOCKS_PER_SEC;

static void reshuffle(double **u,int N){
  int i;
  double *t = u[N-1];
  for(i = N-1; i; --i) u[i] = u[i-1];
  u[0] = t;
}

void SubCycle(Domain *Omega, int Order){
  register int i,j,N;
  double   dt  = Omega->dt;
  double   cfl;
  Element_List *U,*V,*W,*Uf,*Vf,*Wf;
  int      eidmax,Je = iparam("INTYPE");
  double   time = dparam("t");
  double   alpha[3];
  int      MinSteps = iparam("MinSteps");

  U   = Omega->U;
  V   = Omega->V;
  W   = Omega->W;
  Uf  = Omega->Uf;
  Vf  = Omega->Vf;
  Wf  = Omega->Wf;

  st = clock();
  dzero(7,ts,1);

#if 0
  cfl = cfl_checker(Omega,dt);
#else
  cfl = full_cfl_checker(Omega,dt,&eidmax);
#endif
  AddTime(ts[0]);  
  N = (int) (cfl/0.6) + 1;
  N = max(N,MinSteps);
  dt /= (double)N;

#if 0 
  ROOTONLY fprintf(stdout,"N: %d, cfl = %lf\n",(int) N, cfl/N);
#else
  ROOTONLY fprintf(stdout,"N: %d, cfl = %lf (orig cfl = %lf in %d)\n",
		   (int) N, cfl/N,cfl,eidmax);
#endif
  
  // put  u^n into Lagrangian vel ul
  dcopy(U->htot,Omega->u[0],1,Omega->ul[0],1);
  dcopy(V->htot,Omega->v[0],1,Omega->vl[0],1);
  dcopy(W->htot,Omega->w[0],1,Omega->wl[0],1);
  
  // advance Omega->u, Oemga->v by N*dt 
  // this is main part of subcycling 
  for(i = 0; i < Order; ++i)
    Advance(N,dt,Omega->ul[i],Omega->vl[i],Omega->wl[i],Omega,Order);
  
#ifdef TIMING
  ROOTONLY{
    fprintf(stdout,"\t cfl_checker... Took %lf second\n",ts[0]);
    fprintf(stdout,"\t Convective.... Took %lf second\n",ts[1]);
    fprintf(stdout,"\t Fill edges.... Took %lf second\n",ts[2]);
    fprintf(stdout,"\t Upwinding..... Took %lf second\n",ts[3]);
    fprintf(stdout,"\t Add Surf Cont. Took %lf second\n",ts[4]);
    fprintf(stdout,"\t Projection.... Took %lf second\n",ts[5]);
    fprintf(stdout,"\t Avg./Euler.... Took %lf second\n",ts[6]);
  }
#endif

  // calculate explicit part of backwards approximation of
  // Du/Dt*dt and store in [Je-1] for viscous solve
  getalpha(alpha);
  dsmul(U->htot, alpha[Order-1], Omega->ul[Order-1], 1, Omega->ul[Order-1], 1);
  for(i = 0; i < Order-1; ++i)
    daxpy(U->htot, alpha[i], Omega->ul[i], 1, Omega->ul[Order-1], 1);

  dsmul(V->htot, alpha[Order-1], Omega->vl[Order-1], 1, Omega->vl[Order-1], 1);
  for(i = 0; i < Order-1; ++i)
    daxpy(V->htot, alpha[i], Omega->vl[i], 1, Omega->vl[Order-1], 1);

  dsmul(W->htot, alpha[Order-1], Omega->wl[Order-1], 1, Omega->wl[Order-1], 1);
  for(i = 0; i < Order-1; ++i)
    daxpy(W->htot, alpha[i], Omega->wl[i], 1, Omega->wl[Order-1], 1);

  // Put ul into Uf for pressure BC's and viscous solve
  dcopy(U->htot,Omega->ul[Order-1],1,Uf->base_h,1);
  dcopy(V->htot,Omega->vl[Order-1],1,Vf->base_h,1);
  dcopy(W->htot,Omega->wl[Order-1],1,Wf->base_h,1);

  // put current BC into Uf for pressure BC's  
  Fill_BC_sides (Uf, Omega->Ubc,Omega,time+Omega->dt);
  Fill_BC_sides (Vf, Omega->Vbc,Omega,time+Omega->dt);
  Fill_BC_sides (Wf, Omega->Wbc,Omega,time+Omega->dt);

  // copy in u^n into u
  dcopy(U->htot,Omega->u[0],1,U->base_h,1);
  dcopy(V->htot,Omega->v[0],1,V->base_h,1);
  dcopy(W->htot,Omega->w[0],1,W->base_h,1);
  
  SetPBCs(Omega);

  dcopy(U->htot,Omega->ul[Order-1],1,U->base_h,1);
  dcopy(V->htot,Omega->vl[Order-1],1,V->base_h,1);
  dcopy(W->htot,Omega->wl[Order-1],1,W->base_h,1);

  reshuffle(Omega->ul,Je);
  reshuffle(Omega->vl,Je);
  reshuffle(Omega->wl,Je);

  reshuffle(Omega->u,Je+1);
  reshuffle(Omega->v,Je+1);
  reshuffle(Omega->w,Je+1);
}

static void Advance (int N, double dt, double *u, double *v, double *w,
		     Domain *Omega,  int Je){
  register int i,j;
  Element_List *U,*V,*W,*Uf,*Vf,*Wf,*Vcomm[6];
  double time = dparam("t");
  double Beta[3],*uf,*vf,*wf;
  int Iord = Je;

  U  = Omega->U;
  V  = Omega->V;
  W  = Omega->W;
  Uf = Omega->Uf;
  Vf = Omega->Vf;
  Wf = Omega->Wf;
  
  uf = Omega->uf[0];
  vf = Omega->vf[0];
  wf = Omega->wf[0];
  
#ifdef PARALLEL
  Vcomm[0] = U;  Vcomm[1] = V;  Vcomm[2] = W;
  Vcomm[3] = Uf; Vcomm[4] = Vf; Vcomm[5] = Wf;
#endif

#if 1 // if 0 then just use AB scheme
  //Do first substep using RK  
  dcopy(U->htot,u,1,U->base_h,1);
  dcopy(V->htot,v,1,V->base_h,1);
  dcopy(W->htot,w,1,W->base_h,1);
    
  // advance using  Extrapolated field
  // DG convection
  
  // extrapolate 
  ExtField(U->htot,0,Uf->base_h,Omega->u,Omega->dt,min(Omega->step,Je));
  ExtField(U->htot,0,Vf->base_h,Omega->v,Omega->dt,min(Omega->step,Je));
  ExtField(W->htot,0,Wf->base_h,Omega->w,Omega->dt,min(Omega->step,Je));
  Uf->Set_state('p');
  Vf->Set_state('p');
  Wf->Set_state('p');
  
  // fill edges with u
  Fill_sides            (U,U,Omega->Ubc,Omega,time);
  Fill_sides            (V,V,Omega->Vbc,Omega,time);
  Fill_sides            (W,W,Omega->Wbc,Omega,time);

  Fill_sides            (Uf,Uf,Omega->Ubc,Omega,time);
  Fill_sides            (Vf,Vf,Omega->Vbc,Omega,time);
  Fill_sides            (Wf,Wf,Omega->Wbc,Omega,time);
  AddTime(ts[2]);

#ifdef PARALLEL // Do parallel communication all at once
  DO_PARALLEL  exchange_sides(6,Vcomm);  
#endif  


  Upwind_edge_advection (U,V,W,Uf,Vf,Wf);
  Zero_Symm_BCs(Omega->Ubc,U,V,W,Uf,Vf,Wf);

  AddTime(ts[3]);

  CdgradV (Omega); 
  AddTime(ts[1]);
  
  // put upwinded flux minus local flux (stored in U) into Uf 
  Add_surface_contrib   (U,Uf);
  Add_surface_contrib   (V,Vf);
  Add_surface_contrib   (W,Wf);
  AddTime(ts[4]);
  
  InnerProduct_Orth     (Uf,Uf);
  InnerProduct_Orth     (Vf,Vf);  
  InnerProduct_Orth     (Wf,Wf);  
  
  Jtransbwd_Orth        (Uf,Uf);
  Jtransbwd_Orth        (Vf,Vf);
  Jtransbwd_Orth        (Wf,Wf);
  
  AddTime(ts[5]);
  
  if(Je-1){
    // store previous flux
    dcopy(U->htot,Uf->base_h,1,Omega->uf[0],1);
    dcopy(V->htot,Vf->base_h,1,Omega->vf[0],1);
    dcopy(W->htot,Wf->base_h,1,Omega->wf[0],1);
    
    // calculate average velocity
    AverageVel (Je-1,Omega->uf,BarrayA[Je-2][0],Uf->base_h,U->htot);
    AverageVel (Je-1,Omega->vf,BarrayA[Je-2][0],Vf->base_h,V->htot);
    AverageVel (Je-1,Omega->wf,BarrayA[Je-2][0],Wf->base_h,W->htot);
    
    EulerStep  (U->htot,u,Uf->base_h,U->base_h,dt);
    EulerStep  (V->htot,v,Vf->base_h,V->base_h,dt);
    EulerStep  (W->htot,w,Wf->base_h,W->base_h,dt);
    AddTime(ts[6]);
    
    // DG convection

    // interpolate 
    ExtField(U->htot,BarrayC[Je-2][0]*dt,Uf->base_h,
	     Omega->u,Omega->dt,min(Omega->step,Je));
    ExtField(U->htot,BarrayC[Je-2][0]*dt,Vf->base_h,
	     Omega->v,Omega->dt,min(Omega->step,Je));
    ExtField(W->htot,BarrayC[Je-2][0]*dt,Wf->base_h,
	     Omega->w,Omega->dt,min(Omega->step,Je));
    Uf->Set_state('p');
    Vf->Set_state('p');
    Wf->Set_state('p');
    
    // fill edges with u
    Fill_sides  (U,U,Omega->Ubc,Omega,time+BarrayC[Je-2][0]*dt);
    Fill_sides  (V,V,Omega->Vbc,Omega,time+BarrayC[Je-2][0]*dt);
    Fill_sides  (W,W,Omega->Wbc,Omega,time+BarrayC[Je-2][0]*dt);
    // fill edges with uf
    Fill_sides   (Uf,Uf,Omega->Ubc,Omega,time+BarrayC[Je-2][0]*dt);
    Fill_sides   (Vf,Vf,Omega->Vbc,Omega,time+BarrayC[Je-2][0]*dt);
    Fill_sides   (Wf,Wf,Omega->Wbc,Omega,time+BarrayC[Je-2][0]*dt);
    AddTime(ts[2]);

#ifdef PARALLEL // Do parallel communication all at once
    DO_PARALLEL  exchange_sides(6,Vcomm);  
#endif  
   
    Upwind_edge_advection (U,V,W,Uf,Vf,Wf);
    Zero_Symm_BCs(Omega->Ubc,U,V,W,Uf,Vf,Wf);
    AddTime(ts[3]);

    // new flux evaluation
    CdgradV (Omega); 
    AddTime(ts[1]);
    
    Add_surface_contrib (U,Uf);
    Add_surface_contrib (V,Vf);
    Add_surface_contrib (W,Wf);
    AddTime(ts[4]);
    
    InnerProduct_Orth  (Uf,Uf);
    InnerProduct_Orth  (Vf,Vf);  
    InnerProduct_Orth  (Wf,Wf);  
    
    Jtransbwd_Orth     (Uf,Uf);
    Jtransbwd_Orth     (Vf,Vf);
    Jtransbwd_Orth     (Wf,Wf);
    
    AddTime(ts[5]);
    
    // store previous flux
    dcopy(U->htot,Uf->base_h,1,Omega->uf[1],1);
    dcopy(V->htot,Vf->base_h,1,Omega->vf[1],1);
    dcopy(W->htot,Wf->base_h,1,Omega->wf[1],1);
    
    AverageVel (Je,Omega->uf,BarrayB[Je-2],Uf->base_h,Uf->htot);
    AverageVel (Je,Omega->vf,BarrayB[Je-2],Vf->base_h,Vf->htot);
    AverageVel (Je,Omega->wf,BarrayB[Je-2],Wf->base_h,Wf->htot);

    // replace previous flux in uf for AB scheme
    dcopy(U->htot,Omega->uf[1],1,uf,1);
    dcopy(V->htot,Omega->vf[1],1,vf,1);
    dcopy(W->htot,Omega->wf[1],1,wf,1);
  }
  
  EulerStep  (U->htot,u,Uf->base_h,u,dt);
  EulerStep  (V->htot,v,Vf->base_h,v,dt);
  EulerStep  (W->htot,w,Wf->base_h,w,dt);
  AddTime(ts[6]);

  // finish off using AB - max of order 2 at present
  if(Iord == 3) fprintf(stderr,"Warning 3rd order is not debugged\n");

  set_order_CNAB(Iord);
  getbeta(Beta);
  for(i = 1; i < N; ++i){
#else
  for(i = 0; i < N; ++i){

    if(i)
      set_order_CNAB(Iord);
    else
      set_order_CNAB(1);
    
    getbeta(Beta);
#endif

    dcopy(U->htot,u,1,U->base_h,1);
    dcopy(V->htot,v,1,V->base_h,1);
    dcopy(W->htot,w,1,W->base_h,1);
    
    // DG convection

    // interpolate 
    ExtField(U->htot,i*dt,Uf->base_h,Omega->u,Omega->dt,min(Omega->step,Je));
    ExtField(U->htot,i*dt,Vf->base_h,Omega->v,Omega->dt,min(Omega->step,Je));
    ExtField(W->htot,i*dt,Wf->base_h,Omega->w,Omega->dt,min(Omega->step,Je));
    Uf->Set_state('p');
    Vf->Set_state('p');
    Wf->Set_state('p');
    
    // fill edges with u
    Fill_sides            (U,U,Omega->Ubc,Omega,time+i*dt);
    Fill_sides            (V,V,Omega->Vbc,Omega,time+i*dt);
    Fill_sides            (W,W,Omega->Wbc,Omega,time+i*dt);
    // fill edges with uf
    Fill_sides            (Uf,Uf,Omega->Ubc,Omega,time+i*dt);
    Fill_sides            (Vf,Vf,Omega->Vbc,Omega,time+i*dt);
    Fill_sides            (Wf,Wf,Omega->Wbc,Omega,time+i*dt);
    AddTime(ts[2]);
    
#ifdef PARALLEL // Do parallel communication all at once
    DO_PARALLEL  exchange_sides(6,Vcomm);  
#endif  

    Upwind_edge_advection (U,V,W,Uf,Vf,Wf);
    Zero_Symm_BCs(Omega->Ubc,U,V,W,Uf,Vf,Wf);
    AddTime(ts[3]);

    // calcualte forcing within element 
    CdgradV (Omega); 
    AddTime(ts[1]);
    
    // put upwinded flux minus local flux (stored in U) into Uf 
    Add_surface_contrib   (U,Uf);
    Add_surface_contrib   (V,Vf);
    Add_surface_contrib   (W,Wf);
    AddTime(ts[4]);
    
    InnerProduct_Orth     (Uf,Uf);
    InnerProduct_Orth     (Vf,Vf);  
    InnerProduct_Orth     (Wf,Wf);  
    
    Jtransbwd_Orth        (Uf,Uf);
    Jtransbwd_Orth        (Vf,Vf);
    Jtransbwd_Orth        (Wf,Wf);

    AddTime(ts[5]);
    
    // integrate using Adams Bashforth
    daxpy(U->htot,dt*Beta[0],Uf->base_h,1,u,1);
    daxpy(V->htot,dt*Beta[0],Vf->base_h,1,v,1);
    daxpy(W->htot,dt*Beta[0],Wf->base_h,1,w,1);
    
    if(Iord-1){ // 2nd order terms
      daxpy(U->htot,dt*Beta[1],uf,1,u,1);
      daxpy(V->htot,dt*Beta[1],vf,1,v,1);
      daxpy(W->htot,dt*Beta[1],wf,1,w,1);
      
      dcopy(U->htot,Uf->base_h,1,uf,1);
      dcopy(V->htot,Vf->base_h,1,vf,1);
      dcopy(W->htot,Wf->base_h,1,wf,1);
    }
    AddTime(ts[6]);
  }

  set_order(Je);
}

// needs to be process into Element_list structure 
/* 
   This function interploltes the value in Bc at the quadrature points
   and the subtracts this value from the field "in" and put the value 
   into field "out"
*/
static void SubtractBC(Element_List *in, Element_List *out, Bndry *Bc,
		       double dt){
  int eid;
  Bndry *B;
  double gamma = getgamma();

  /* set up boundary conditions  */
  for(B = Bc; B; B = B->next)
    if((B->type == 'V')||(B->type == 'W')){
      eid  = B->elmt->id;
      in->flist[eid]->SubtractBC(1/dt,gamma/dt,B->face,out->flist[eid]);
    }
}


/* Sum Uout = sum_N Uf[i]*Barray[i]                                   */
/* Note: Summation is in reverse order due to second call in SubCycle */
static void AverageVel(int N, double  **Uf, double *Barray, 
		       double  *Uout, int ntot){
  int i,k,l;
  
  dsmul(ntot,Barray[N-1],Uf[N-1],1,Uout,1);
    
  for(i = N-2; i >=0; --i)
    daxpy(ntot,Barray[i],Uf[i],1,Uout,1);
}

/* ------------------------------------------------------------------*
 * This function fills the edges array 'edge->h' with the physical
 * values of the state vector. It also fills the edge->link->h array
 * of any boundary value with the appropriate boundary conditions. All
 * values are interpolated to the gauss points of order edge->qedg.
 * This assumes that the state vector Us is in physical space.
 *
 * W - the same as reflexive but with no characteristic b.c. 
 * 
 * S - outflow b.c., is updated every time step 
 * ----------------------------------------------------------------- */

static void update_bndry_h(Bndry *Bc, double t, Domain *Omega);
static void Fill_sides(Element_List *U, Element_List *Uf, 
		       Bndry *Ubc, Domain *Omega, double t){
  Element      *E,*Ef;

  if(U->fhead->state != 'p')
    error_msg(Fill_sides: state vector not in physical space);
  
  /* set interior edge points */
  for(E=U->fhead,Ef = Uf->fhead;E;E=E->next,Ef=Ef->next)
    Ef->fill_edges(E->h_3d[0][0], NULL, NULL);
  
  Fill_BC_sides(Uf,Ubc, Omega,t);
    
#ifdef PARALLELNOT
  DO_PARALLEL  exchange_sides(1,&Uf);  
#endif  
}


static void Fill_BC_sides(Element_List *Uf, Bndry *Ubc, 
			  Domain *Omega, double t){
  Bndry        *B;
  Face         *f;
  int          eid;
  int          tvary = option("tvarying");

  if(Omega->Tfun && Omega->Tfun->type == Womersley)
    SetWomSol(Omega,t,0);

  dparam_set("t",t);// needs to be set for explicit function 

  /* set up boundary conditions  */
  for(B = Ubc; B; B = B->next){
    if(tvary)
      update_bndry_h(B,t,Omega);
    
    eid = B->elmt->id;
    f   = Uf->flist[eid]->face[B->face].link;
    dcopy(f->qface*f->qface,B->phys_val_g,1,f->h,1);
  }
}
// evaluate boundary string and put into gauss space in Bc

static  void update_bndry_h(Bndry *Bc, double time, Domain *Omega){
  int     qt,fq,qa,qb;
  double  *f,**ima,**imb;
  Coord   X;
  Element *E = Bc->elmt;
  
  fq   = E->face[Bc->face].qface;

  qa = E->qa;
  qb = E->qb;
  qt = qa*qb;

  getim(qa,fq,&ima,a2g); 
  if(E->Nfverts(Bc->face) == 3)
    getim(qb,fq,&imb,b2g); 
  else
    getim(qb,fq,&imb,a2g); 
  
  f   = dvector(0,qt-1);

  if(Omega->Tfun && Omega->Tfun->type == Womersley){
    int    cnt,nfv,l,lf,i;
    double *store = dvector(0,QGmax*QGmax-1);
    
    nfv = Bc->elmt->Nfverts(Bc->face);
    lf  = Bc->elmt->face[Bc->face].l;
    if(nfv == 3)
      lf = lf*(lf+1)/2;
    else
      lf = lf*lf;
    
    for(i = 0; i < nfv; ++i)
      store[i] = Bc->bvert[i];
    
    cnt = nfv;
    for(i = 0; i < nfv; ++i){
      l = E->edge[E->ednum(Bc->face,i)].l;
      dcopy(l,Bc->bedge[i],1,store+cnt,1);
      cnt += l;
    }
    
    if(lf)
      dcopy(lf,Bc->bface[0],1,store+cnt,1);

    E->Jbwdfac1(Bc->face,store,f);

    free(store);
  }
  else{
    X.x = dvector(0,3*(qt+qa)-1);
    X.y = X.x + qt+qa;
    X.z = X.y + qt+qa; 
    
    vector_def("x y z", Bc->bstring);
    dparam_set("t",time);
    
    E->GetFaceCoord(Bc->face,&X);
    E->InterpToFace1(Bc->face, X.x, f);
    dcopy(qa*qb,f,1,X.x,1);
    E->InterpToFace1(Bc->face, X.y, f);
    dcopy(qa*qb,f,1,X.y,1);
    E->InterpToFace1(Bc->face, X.z, f);
    dcopy(qa*qb,f,1,X.z,1);
    
    vector_set(qt,X.x,X.y,X.z,f);
    free(X.x);
  }

  /* store gaussian distribution for euler part */
  Interp2d(*ima,*imb, f,qa,qb,Bc->phys_val_g,fq,fq);

  free(f);
    
  return;
}

void Add_surface_contrib(Element_List *Us, Element_List *Uf){
  Element *E, *Ef;
  
  for(E = Us->fhead, Ef = Uf->fhead; E; E = E->next, Ef = Ef->next)
    E->Add_Surface_Contrib(Ef, Ef->h_3d[0][0], 'n');
} 

/* calculate upwinded normal flux normal and put in Uf. The interior
   flux is also put in U */
static void Upwind_edge_advection(Element_List *U,  Element_List *V,
				  Element_List *W,  Element_List *UF, 
				  Element_List *VF, Element_List *WF){
  register     int i,j;
  int          *loc_id,*adj_id,qface,qt;
  Coord        *normal;
  double       sloc,sadj,floc,fadj;
  double       un, un_loc, un_adj;  
  int          *lid, *aid;
  Element      *Us, *Vs, *Ws, *Uf, *Vf, *Wf;
  Face         *LF, *F;
  
  /* only need to calculate the edges which have edge->link defined */
  for(Us = U->fhead,  Vs = V->fhead,  Ws = W->fhead, Uf = UF->fhead, 
      Vf = VF->fhead, Wf = WF->fhead; Us; Us = Us->next, Vs = Vs->next, 
      Ws = Ws->next,  Uf = Uf->next,  Vf = Vf->next, Wf = Wf->next){
    
    for(i = 0; i < Us->Nfaces; ++i){
      F  = Us->face + i;
      LF = Us->face[i].link;
      if(LF)
	if((F->eid > LF->eid)||(F->eid == LF->eid && LF->id >= i)){
	  
	  normal = F->n; 
	  qface  = F->qface;
	  
	  qt = qface*qface;

	  if(Us->Nfverts(i) == 3){
	    lid = Tri_nmap(qface,F->con);
	    aid = Tri_nmap(qface,LF->con);
	  }
	  else{
	    lid = Quad_nmap(qface,F->con);
	    aid = Quad_nmap(qface,LF->con);
	  }

	  /* finally put upwinded flux calculation in here */
	  for(j = 0; j < qt; ++j){

	    // calculate normal velcoity based on avaerage of both sides 
	    un_loc = Uf->face[i].h[lid[j]]*normal->x[lid[j]] + 
	      Vf->face[i].h[lid[j]]*normal->y[lid[j]] +
	      Wf->face[i].h[lid[j]]*normal->z[lid[j]];
	    
	    un_adj = Uf->face[i].link->h[aid[j]]*normal->x[lid[j]] + 
	      Vf->face[i].link->h[aid[j]]*normal->y[lid[j]] + 
	      Wf->face[i].link->h[aid[j]]*normal->z[lid[j]];
	    
	    un = 0.5*(un_loc + un_adj);
	    
	    /* calculate upwinded flux for Uf */
	    if(un < 0){
	      // negate flux to be consistent with RHS term
	      floc =  -un_adj*Us->face[i].link->h[aid[j]]; 
	      fadj =  -floc;
	    }
	    else{      
	      // negate flux to be consistent with RHS term
	      floc = -un_loc*Us->face[i].h[lid[j]]; 
	      fadj = -floc;
	    }
	    
	    Uf->face[i].h[lid[j]]       = floc;
	    Uf->face[i].link->h[aid[j]] = fadj;
	    
	    /* calculate upwinded flux for Vf*/
	    if(un < 0){
	      // negate flux to be consistent with RHS term
	      floc = -un_adj*Vs->face[i].link->h[aid[j]]; 
	      fadj = -floc;
	    }
	    else{      
	      // negate flux to be consistent with RHS term
	      floc = -un_loc*Vs->face[i].h[lid[j]]; 
	      fadj = -floc;
	    }
	    
	    Vf->face[i].h[lid[j]]       = floc;
	    Vf->face[i].link->h[aid[j]] = fadj;
	    
	    /* calculate upwinded flux for Wf*/
	    if(un < 0){
	      // negate flux to be consistent with RHS term
	      floc = -un_adj*Ws->face[i].link->h[aid[j]]; 
	      fadj = -floc;
	    }
	    else{      
	      // negate flux to be consistent with RHS term
	      floc = -un_loc*Ws->face[i].h[lid[j]]; 
	      fadj = -floc;
	    }
	    
	    Wf->face[i].h[lid[j]]       = floc;
	    Wf->face[i].link->h[aid[j]] = fadj;
	    
	    /* store local flux (negative) in Us fields */
	    Us->face[i].h[lid[j]]        *= -un_loc;
	    Us->face[i].link->h[aid[j]]  *=  un_adj;

	    /* store local flux (negative) in Vs fields */
	    Vs->face[i].h[lid[j]]        *= -un_loc;
	    Vs->face[i].link->h[aid[j]]  *=  un_adj;
	    
	    /* store local flux (negative) in Ws fields */
	    Ws->face[i].h[lid[j]]        *= -un_loc;
	    Ws->face[i].link->h[aid[j]]  *=  un_adj;
	  }
	}
    }
  }
}
 
static void Zero_Symm_BCs(Bndry *Ubc, Element_List *U, Element_List *V, 
  Element_List *W, Element_List *Uf, Element_List *Vf, Element_List *Wf){

  int face, eid, qt;
  Bndry *B;

  for(B = Ubc; B; B = B->next){
    if(B->usrtype == 'Z'){
      face = B->face;
      eid = B->elmt->id;
      qt = U->flist[eid]->face[face].qface;
      qt = qt*qt;
      
      dzero(qt,U->flist[eid]->face[face].link->h,1);
      dzero(qt,V->flist[eid]->face[face].link->h,1);
      dzero(qt,W->flist[eid]->face[face].link->h,1);
      dzero(qt,Uf->flist[eid]->face[face].link->h,1);
      dzero(qt,Vf->flist[eid]->face[face].link->h,1);
      dzero(qt,Wf->flist[eid]->face[face].link->h,1);
    }
 }
}
/* Function to project orthogonal bases by setting up matrix and 
   multiplying all element storage. Assumes fixed order and clustering
   of types */
static int SetupProjection(double **P, Element *U);

static void Project(Element_List *U, Element_List *Uf, double *wk){
  static int init = 1, Ptot;
  static double *P; 
  int nel = U->nel;
  
  if(init){
    Ptot = SetupProjection(&P, U->fhead);
    init = 0;
  }
  
  dgemm('n','n',Ptot,nel,Ptot,1.0,P,Ptot,U->base_h,Ptot,0.0,wk, Ptot);
  dcopy(U->htot,wk,1,Uf->base_h,1);
}

static int SetupProjection(double **P, Element *U){
  int i;
  int L = U->dgL;
  int qtot = U->qtot;
  double *p;
  
  P[0] =  dvector(0,qtot*qtot-1); 
  dzero(qtot*qtot,*P,1);

  for(i = 0; i < qtot; ++i){
    p = P[0] + i*qtot;
    p[i] = 1;
    U->Ofwd(p,p,L);
    U->Obwd(p,p,L);
  }

  return qtot;
}

/* extrapolate field using equally time spaced field un,un-1,un-2, (at
   dt intervals) to time n+t at order Ord */
static void ExtField(int htot, double t, double *un, double **up, double dt, 
		   int ord){
  register int i,j;
  double l[4];
		
  // calculate lagrange interpolants
  dfill(4,1,l,1);
  for(i = 0; i <= ord; ++i)
    for(j = 0; j <= ord; ++j)
      if(i != j){
	l[i] *= (j*dt+t);
	l[i] /= (j*dt-i*dt);
      }

  dsmul(htot,l[0],up[0],1,un,1);
  for(i = 1; i <= ord; ++i)
    daxpy(htot,l[i],up[i],1,un,1);
  
  return;
}
