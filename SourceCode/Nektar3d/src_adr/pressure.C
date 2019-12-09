/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $Source: /homedir/cvs/Nektar/Nektar3d/src/pressure.C,v $
 * $Revision: 1.2 $
 * $Date: 2006/05/08 10:01:39 $ 
 * $Author: ssherw $ 
 * $State: Exp $ 
 *---------------------------------------------------------------------------*/
#include "nektar.h"
#include "string.h"

#ifdef PBC_1D
/* define communicators for inlets and outlets */
extern MPI_Comm  Communicator_inlet[26];
extern MPI_Comm  Communicator_outlet[26];
extern MPI_Comm  Communicator_inlet_plus_ROOT[26];
extern MPI_Comm  Communicator_outlet_plus_ROOT[26];
#endif

/* external function , defined in flowrate.C */
double Flux(Domain *omega, char *label);
double Flux(Domain *omega, int *ID_faces_per_outlet, int Nfaces);
double FluxInlet(Domain *omega, char *label);
double PressureMean(Domain *omega,  char *label);
double PressureMeanInlet(Domain *omega,  char *label);

/*
 * Build boundary conditions for the pressure
 */

Bndry *BuildPBCs(Element_List *P, Bndry *temp)
{
  int      nbcs = 0, Je = iparam("INTYPE");
  double   pval;
  Bndry    *Pbc;
  Element  *E;
  register  int i,j;

  /* Count the number of boundaries to create */
  if(!temp) 
    return (Bndry*)0;

  while (temp[nbcs++].next);

  Pbc = (Bndry*)calloc(nbcs, sizeof(Bndry));

  for(i = 0; i < nbcs; ++i) {
    Pbc[i].id          = i;
    Pbc[i].type        = temp[i].type;
    Pbc[i].usrtype     = temp[i].usrtype;
    Pbc[i].elmt        = P->flist[temp[i].elmt->id];
    Pbc[i].face        = temp[i].face;
    Pbc[i].next        = Pbc + i + 1;
    Pbc[i].bstring     = temp[i].bstring;
  }
  Pbc[nbcs-1].next = (Bndry*) NULL;
  
  
  /* Translate the boundary conditions */
  
  for(i = 0; i < nbcs; ++i) {
    switch (Pbc[i].usrtype) {
    case 'O':
      Pbc[i].type = 'o';                       /* fall through */
      E = Pbc[i].elmt;
      E->Surface_geofac(Pbc+i);
      E->MemBndry(Pbc+i,Pbc[i].face,(iparam("EQTYPE") == Rotational)? Je:1);
      Pbc[i].blabel  = temp[i].blabel;      
#ifndef VELINTERP
      vector_def("",Pbc[i].bstring);
      vector_set(1,&pval);
#else
      pval = 0.0;
#endif

      for(j = 0; j < E->Nfverts(Pbc[i].face); ++j) 
	Pbc[i].bvert[j] = pval;
      
      break;     	  
    case 'V': case 'v':          /* allocate additional multi-step storage */
			Pbc[i].usrtype = 'v';
    case 'W': case 'F':
      Pbc[i].type = 'F';
      Pbc[i].elmt->Surface_geofac(Pbc+i);
      Pbc[i].elmt->MemBndry(Pbc+i,Pbc[i].face,Je);

      if (Pbc[i].usrtype == 'v')
        Pbc[i].blabel  = temp[i].blabel;

      break;
    default:
      Pbc[i].elmt->MemBndry(Pbc+i,Pbc[i].face,1);
      break;
    }
  }

  return  bsort(Pbc,nbcs);
  
}

/* this function resets the global mesh and boundary conditions so that
   they are correctly set for the pressure element structure */
void Set_Global_Pressure(Element_List *Mesh, Bndry *Meshbcs){
  register int i,j,id;
  Element *U;
  Bndry   *B;
  int flag = 0;

  /* reset solve mask */
  for(U = Mesh->fhead; U; U = U->next)
    for(i = 0; i < U->Nverts; ++i)
      U->vert[i].solve = 1;

  DO_PARALLEL{ /* set interface solve masks to 2 */
    if((Mesh->fhead->dim() == 3)&&(pllinfo.partition)){
      for(U = Mesh->fhead; U; U = U->next)
	for(i = 0; i < U->Nfaces; ++i)
	  if(U->face[i].link)
	    if(pllinfo.partition[U->id] 
	       != pllinfo.partition[U->face[i].link->eid])
	      for(j = 0; j < U->Nfverts(i); ++j){
		id = U->vnum(i,j);
		U->vert[id].solve *= 2;
		// clamp value for mult calls
		U->vert[id].solve = min(U->vert[id].solve,2); 
	      }
    }
  }

  /* loop through boundary conditions and set outflow boundaries */  
  for(B = Meshbcs; B; B = B->next){
    switch (B->type) {
    case 'O':
      B->type = 'o';                       /* fall through */
      U = B->elmt;
      for(i = 0; i < U->Nfverts(B->face); ++i) 
	U->vert[U->vnum(B->face,i)].solve = 0;
      flag = 1;
      break;     	  
    case 'V': case 'v':
    case 'W':  case 'Z':
      B->type = 'F';
      break;
    }
  }
}

/* replace vertex, edge and face numbering using the global numbering
   scheme specified in GU */

void Replace_Numbering(Element_List *UL, Element_List *GUL){
  Element *U, *GU;
  int i,k;

  for(k = 0; k < UL->nel; ++k){
    U  = UL->flist[k];
    GU = GUL->flist[pllinfo.eloop[k]];
    for(i = 0; i < U->Nverts; ++i){
      U->vert[i].gid   = GU->vert[i].gid;
      U->vert[i].solve = GU->vert[i].solve;
    }

    for(i = 0; i < U->Nedges; ++i)
      U->edge[i].gid   = GU->edge[i].gid;

    for(i = 0; i < U->Nfaces; ++i)
      U->face[i].gid   = GU->face[i].gid;
  }
}


/* replace vertex, edge and face numbering using the numbering
   scheme specified in GU */

void Replace_Local_Numbering(Element_List *UL, Element_List *GUL){
  Element *U, *GU;
  int i,k;

  for(k = 0; k < UL->nel; ++k){
    U  = UL->flist[k];
    GU = GUL->flist[k];
    for(i = 0; i < U->Nverts; ++i){
      U->vert[i].gid   = GU->vert[i].gid;
      U->vert[i].solve = GU->vert[i].solve;
    }

    for(i = 0; i < U->Nedges; ++i)
      U->edge[i].gid   = GU->edge[i].gid;
    
    for(i = 0; i < U->Nfaces; ++i)
      U->face[i].gid   = GU->face[i].gid;
  }
}

/* ------------------------------------------------------------------------- *
 * SetPBCs() -  Set boundary conditions for the pressure                     *
 *                                                                           *
 *                   Je                     1                                *
 *        dP/dn = SUM    beta[q] * [ N(u) - - curl ( curl u ) ] * n          *
 *                   q=0                    R                                *
 *                                                                           *
 * where n is the unit outward normal along the edge, u is the velocity      *
 * field, and N(u) are the non-linear terms in the momentum equation.  This  *
 * routine computes the RHS (already integrated) of the above equation at    *
 * the current time level and saves it in the corresponding Bedge for Press. *
 * ------------------------------------------------------------------------- */

static void CalcPbc (Bndry *B, Element_List *V[3], Element_List *N[3], double nu);
static void IntPbc  (Bndry *B, int Je);
static void outflow (Bndry *B, Domain *omega);
static void AddAcc  (Bndry *B, Domain *Omega);
static void AddDuDt (Bndry *B, Element_List *N[3], double dt);
static void PbcOutflowStable(Bndry *B, Domain *omega);

void SetPBCs(Domain *omega){
  Bndry    *Pbc  =  omega->Pbc;
  Element_List  *V[3], *N[3];
  int       Je   =  iparam("INTYPE");
  double    nu   =  dparam("KINVIS");
  int       tvary = option("tvarying");
  SLVTYPE   slvtype = (SLVTYPE) iparam("SLVTYPE");

  V[0] =  omega->U,  V[1] =  omega->V, V[2] = omega->W;
  N[0] =  omega->Uf, N[1] = omega->Vf, N[2] = omega->Wf;
  
#ifdef PBC_1D
   /* compute and store flow rate though each outlet */

   int i,step;
   int Nout = omega->IMPBC[0].Nout;
   int Ninl = omega->IMPBC[0].Ninl;
   double dtemp;
   double *flowrate = dvector(0,Nout+Ninl-1);
   double *pressure = dvector(0,Nout+Ninl-1);
   double *work = dvector(0,Nout+Ninl-1);

   dtemp = dparam("DELT");
   step = ((int) ((dparam("t")+0.1*dtemp)/dtemp) );

   for (i = 0; i < Nout; ++i)
     flowrate[i] = Flux(omega, omega->IMPBC[0].standard_labels[i]);

   for (i = 0; i < Ninl; ++i)
     flowrate[Nout+i] = FluxInlet(omega, omega->IMPBC[0].standard_labels[i]);

   memset(work,'\0',(Nout+Ninl)*sizeof(double));
   DO_PARALLEL
     gdsum(flowrate,(Nout+Ninl),work);

   if (omega->IMPBC[0].Nimpedance_modes > 0)
      omega->IMPBC[0].ComputePressure();
   else{
			if (iparam("ILPN"))
				omega->IMPBC[0].ComputeCoronaryPressure(flowrate);
			else
	      omega->IMPBC[0].ComputePressureSteady(flowrate);	// Alireza: constant RCR
			//omega->IMPBC[0].ComputePressureRCNonsteady(flowrate);	// Alireza: R(t)C for measured outlet flow waveforms
		}

   if  ( step % iparam("HISSTEP") == 0 ){
     for (i = 0; i < Nout; ++i)
       pressure[i] = PressureMean(omega, omega->IMPBC[0].standard_labels[i]);

     for (i = 0; i < Ninl; ++i)
       pressure[Nout+i] = PressureMeanInlet(omega, omega->IMPBC[0].standard_labels[i]);

     DO_PARALLEL
       gdsum(pressure,(Nout+Ninl),work);

     for (i = 0; i < (Nout+Ninl); ++i)
       pressure[i] /= omega->IMPBC[0].Area[i];

     /* write flowrate and pressure at inlet/outlets */
     ROOTONLY{
       fprintf(omega->flo_file,"%f ",dparam("t"));
       for (i = 0; i < (Nout+Ninl); i++)
         fprintf(omega->flo_file,"%e ",flowrate[i]);
       fprintf(omega->flo_file," \n");
       fflush(omega->flo_file);

       fprintf(omega->pre_file,"%f ",dparam("t"));
       for (i = 0; i < (Nout+Ninl); i++)
         fprintf(omega->pre_file,"%e ",pressure[i]);
       fprintf(omega->pre_file," \n");
       fflush(omega->pre_file);
     }
   }

   free(flowrate); free(pressure); free(work);
#endif

  /* Get the integration coefficients */  
  while (Pbc) {
    switch (Pbc->type) {
    case 'D': case 'N': case 'P':
      break;
      
    case 'o':
#ifdef PBC_1D
     /* use value computed through impedance B.C. */
      double pval;
      pval = omega->IMPBC[0].GetPval(Pbc);
			if (pval == 0.0)
	      for(i = 0; i < Pbc->elmt->Nfverts(Pbc->face); ++i)
  	      Pbc->bvert[i] = pval;
			else
				PbcOutflowStable(Pbc, omega);	// Alireza: stabilization of backflow for the outlet boundary	
#endif
      if (iparam("EQTYPE") == Rotational){
	outflow (Pbc, omega);
	IntPbc  (Pbc, Je);
      }
      break;
      
    case 'F': case 'R': {
      CalcPbc(Pbc,V,N,nu);
      IntPbc (Pbc,Je);
      if(slvtype == SubStep) // this need to be added after integration
	AddDuDt (Pbc, N, omega->dt);
      else
	if(tvary) // add acceleration term 
	  AddAcc (Pbc, omega);
      break;
    }
      
    default:
      error_msg(SetPBCs -- unknown pressure b.c.)
      break;
    }
    
    Pbc = Pbc->next;
  }
  return;
}

static void IntPbc(Bndry *B, int Je){
  register int i,j;
  const    int face = B->face;
  int      ll, nfv;
  double   beta[3];
  Element  *E;
  double   *tmp = dvector(0,QGmax*QGmax-1);

  E = B->elmt;
  getbeta (beta);        
  
  // copy from 3d
  nfv = E->Nfverts(face);

  /* integrate vertices */
  dzero(nfv,tmp,1);
  for(i = 0; i < nfv; ++i)
    for(j = 0; j < Je; ++j)
      tmp[i] += beta[j]*B->bvert[i+j*nfv];
  
  dcopy(nfv * (Je-1), B->bvert, 1, B->bvert+nfv,1);
  dcopy(nfv         ,      tmp, 1, B->bvert,    1);
  
  /* integrate edges */
  for(i = 0; i < nfv; ++i){
    ll = E->edge[E->ednum(face,i)].l;
    dzero(ll,tmp,1);
    for(j = 0; j < Je; ++j)
      daxpy(ll, beta[j], B->bedge[i] + j*ll, 1, tmp, 1);
  
    dcopy(ll * (Je-1), B->bedge[i], 1, B->bedge[i] + ll, 1);
    dcopy(ll         , tmp        , 1, B->bedge[i]     , 1);
  } 
  
  /* integrate face */
  ll = E->face[face].l; 
  ll = (nfv == 3) ? ll*(ll+1)/2 : ll*ll;
  if(ll){
    dzero(ll,tmp,1);
    for(j = 0; j < Je; ++j)
      daxpy(ll, beta[j], *B->bface + j*ll, 1, tmp, 1);
  
    dcopy(ll * (Je-1), *B->bface, 1, *B->bface + ll, 1);
    dcopy(ll         , tmp      , 1, *B->bface     , 1);
  }
  
  free(tmp);
}

/* --------------------------------------------------------------------- *
 * This is a function to caclulate the flux integral of                  *
 *             dP/dn = Int {g, n * [N(u) - 1/Re Curl (Curl U)]}          *
 * Given N and U. The value is then stored in B.                         *
 * --------------------------------------------------------------------- */

static void CalcPbc(Bndry *B, Element_List *V[3], 
		              Element_List *N[3], double nu){


  const    int id  = B->elmt->id;
  Element  *E = B->elmt;
  int       tot = E->qtot;
  double   **Q = dmatrix(0,11,0,tot-1),
            *ux = Q[0], *uy = Q[1], *uz = Q[2],
            *vx = Q[3], *vy = Q[4], *vz = Q[5],
            *wx = Q[6], *wy = Q[7], *wz = Q[8],
            *Qx = Q[9], *Qy = Q[10], *Qz = Q[11];
  static   SLVTYPE   slvtype = (SLVTYPE) iparam("SLVTYPE");

  Element *u  = V[0]->flist[id],  *v = V[1]->flist[id],  *w = V[2]->flist[id];
  Element *Nu = N[0]->flist[id], *Nv = N[1]->flist[id], *Nw = N[2]->flist[id];
  
  
  /* calculate Q = curl(V) */
  u->Grad_d(ux, uy, uz, 'a');
  v->Grad_d(vx, vy, vz, 'a');
  w->Grad_d(wx, wy, wz, 'a');
 
  dvsub  (tot, wy, 1, vz, 1, Qx, 1);
  dvsub  (tot, uz, 1, wx, 1, Qy, 1);
  dvsub  (tot, vx, 1, uy, 1, Qz, 1);

  /* calculate Q = curl(Q) */
  u->Grad_h(Qx, ux, uy, uz, 'a');
  v->Grad_h(Qy, vx, vy, vz, 'a');
  w->Grad_h(Qz, wx, wy, wz, 'a');

  dvsub  (tot, wy, 1, vz, 1, Qx, 1);
  dvsub  (tot, uz, 1, wx, 1, Qy, 1);
  dvsub  (tot, vx, 1, uy, 1, Qz, 1);
  
  if(slvtype == SubStep){ // only do linear term in substep scheme here.
    // just calc -nu*curl(Q)
    dscal (tot, -nu, Qx, 1);
    dscal (tot, -nu, Qy, 1);
    dscal (tot, -nu, Qz, 1);    
  }
  else{ 
    /* calc  Non-linear terms - nu*curl(Q)  */
    dsvtvp (tot, -nu, Qx, 1, **Nu->h_3d, 1, Qx, 1);
    dsvtvp (tot, -nu, Qy, 1, **Nv->h_3d, 1, Qy, 1);
    dsvtvp (tot, -nu, Qz, 1, **Nw->h_3d, 1, Qz, 1);
  }
  /* take dot product with normal */
  u->GetFace(Qx, B->face, ux);
  u->GetFace(Qy, B->face, uy);
  u->GetFace(Qz, B->face, uz);
  
  u->InterpToFace1(B->face,ux,vx);
  u->InterpToFace1(B->face,uy,vy);
  u->InterpToFace1(B->face,uz,vz);
  
  if(u->Nfverts(B->face)==3){
    switch(u->identify()){
    case Nek_Prism: case Nek_Pyr:
      tot = u->qa*u->qc;
      break;
    case Nek_Tet:
      tot = u->qa*u->qb;
      break;
    }
  }
  else
    tot = u->qa*u->qb;

  if(u->curvX){
    dvmul  (tot, B->nx.p, 1, vx, 1, Qx, 1);
    dvvtvp (tot, B->ny.p, 1, vy, 1, Qx, 1, Qx, 1);
    dvvtvp (tot, B->nz.p, 1, vz, 1, Qx, 1, Qx, 1);
  }
  else{
    dsmul  (tot, B->nx.d, vx, 1, Qx, 1);
    dsvtvp (tot, B->ny.d, vy, 1, Qx, 1, Qx, 1);
    dsvtvp (tot, B->nz.d, vz, 1, Qx, 1, Qx, 1);
  }
  
  u->MakeFlux(B,0,Qx);

  free_dmatrix(Q,0,0); 
}

/*
 * PI = p + 1/2 U.U   (outflow boundary conditions)
 */

static void outflow (Bndry *B, Domain *omega){
  const    int id = B->elmt->id, face = B->face;
  Element  *U,*V,*W;
  double   *u = dvector(0,QGmax*QGmax-1),   *v = dvector(0,QGmax*QGmax-1);
  double   *w = dvector(0,QGmax*QGmax-1), *tmp = dvector(0,QGmax*QGmax-1),pval;
  int       q, nfv,vn,i;

  // figure from 3d

  U = omega->U->flist[id];
  V = omega->V->flist[id];
  W = omega->W->flist[id];

  nfv = U->Nfverts(face);

#ifdef PBC_1D
  /* use value computed through impedance B.C. */
  pval = omega->IMPBC[0].GetPval(B);
#else
  vector_def("",B->bstring);
  vector_set(1,&pval);
#endif

  for(i=0;i<nfv;++i){
    vn = U->vnum(face, i);
    B->bvert[i] = pval + 0.5*(U->vert[vn].hj[0]*U->vert[vn].hj[0]+
			      V->vert[vn].hj[0]*V->vert[vn].hj[0]+
			      W->vert[vn].hj[0]*W->vert[vn].hj[0]);
  }
  
  U->GetFace(**U->h_3d, face, u);
  V->GetFace(**V->h_3d, face, v);
  W->GetFace(**W->h_3d, face, w);

  if(nfv==3){
    switch(U->identify()){
    case Nek_Prism: case Nek_Pyr:
      q = U->qa*U->qc;
      break;
    case Nek_Tet:
      q = U->qa*U->qb;
      break;
    }
  }
  else
    q = U->qa*U->qb;

  U->InterpToFace1(face, u, tmp);
  dvmul(q, tmp, 1, tmp, 1, u, 1);
  
  V->InterpToFace1(face, v, tmp);
  dvvtvp(q, tmp, 1, tmp, 1, u, 1, u, 1);
  
  W->InterpToFace1(face, w, tmp);
  dvvtvp(q, tmp, 1, tmp, 1, u, 1, u, 1);
  
  dscal  (q, 0.5,  u, 1);
  dsadd  (q, pval, u, 1, u, 1);

  B->face = 0; 
  U->JtransFace(B,u);
  B->face = face;

  free(u); free(v); free(w); free(tmp);
  return;
}

//  dp     -du . n                           
//  -- +=   --     += - sum_{i=0}^{Je} Delta^i u^i  
//  dn      dt                               
//
// -(u^{n+1} - u^n)/dt   only first order accurate

static double Delta_Int[4] = { 1.0, -1.0, 0.0, 0.0};

static double  Delta_SS[][4] = {
  { 1.0, -1.0, 0.0, 0.0},
  { 1.5, -2.0, 0.5, 0.0},
  { 11/6., -3.0, 3/2., -1/3.}};

void set_delta(int Je){ 
  dcopy(4, Delta_SS[Je-1], 1, Delta_Int, 1);
}

static void AddAcc(Bndry *B, Domain *Omega){
  char   type;
  int    trip = 0;
  static int startup = 1, *PtoVbc;;
  Bndry  *Bc, *Vbc;
  double invdt = 1/Omega->dt;
  int    Je = iparam("INTYPE");

  if(startup){
    int count;
    /* search through pressure BC and find matching Ubc */
    /* count PBC's */

    count = 0;
    for(Bc = Omega->Pbc; Bc; Bc = Bc->next)
      ++count;

    PtoVbc = ivector(0,count-1);

    for(Bc = Omega->Pbc; Bc; Bc = Bc->next)
      for(Vbc = Omega->Ubc; Vbc; Vbc = Vbc->next){
	if((Vbc->elmt->id == Bc->elmt->id)&&(Vbc->face == Bc->face)){
	  PtoVbc[Bc->id] = Vbc->id;
	  break;
	}
      }
    startup = 0;
  }

  type = Omega->Ubc[PtoVbc[B->id]].type;

  if((type  == 'V')||(type == 'M')){ // might need other types here
    register int i,j,n;
    int     eid = B->elmt->id,vbid, nfv;
    int     l,lf,dim,tot,cnt;
    double  **tmp  = dmatrix(0,3,0,QGmax*QGmax-1);
    double  *store = dvector(0,LGmax*LGmax-1);
    Bndry   *VelBc[MAXDIM];
    Element *E;

    E  = B->elmt;

    nfv  = B->elmt->Nfverts(B->face);

    if(nfv==3){
      switch(E->identify()){
      case Nek_Prism: case Nek_Pyr:
	tot = E->qa*E->qc;
	break;
      case Nek_Tet:
	tot = E->qa*E->qb;
	break;
    }
    }
    else
      tot = E->qa*E->qb;
    
    VelBc[0] = Omega->Ubc;
    VelBc[1] = Omega->Vbc;
    VelBc[2] = Omega->Wbc;
    dim      = 3;
    
    dzero(3*QGmax*QGmax,tmp[0],1);
    
    /* get velocity boundary condition matching this pressure bc */
    vbid = PtoVbc[B->id];

    lf = B->elmt->face[B->face].l;
    if(nfv == 3)
      lf = lf*(lf+1)/2;
    else
      lf = lf*lf;

    for(n = 0; n < dim; ++n){
      Vbc = VelBc[n];
      
      /* fill difference of boundary conditions into store */
      for(i = 0; i < nfv; ++i){
	store [i] = 0;
	for(j =0; j < Je+1; ++j)
	  store[i] -= invdt*Delta_Int[j]*Vbc[vbid].bvert[j*nfv+i];
      }
      
      cnt = nfv;
      for(i = 0; i < nfv; ++i){
	l = E->edge[E->ednum(B->face,i)].l;
	dsmul(l, -invdt*Delta_Int[0], Vbc[vbid].bedge[i],1, store+cnt,1);
	for(j =1; j < Je+1; ++j){
	  dsvtvp(l,-invdt*Delta_Int[j], Vbc[vbid].bedge[i]+j*l,1,
		 store+cnt,1,store+cnt,1);
	}
	cnt += l;
      }

      if(lf){
	dsmul(lf, -invdt*Delta_Int[0], Vbc[vbid].bface[0],1, store+cnt,1);
	for(j =1; j < Je+1; ++j)
	  dvsub(lf,Vbc[vbid].bface[0],1,Vbc[vbid].bface[0]+j*lf,1,store+cnt,1);
      }

      /* backward transform */
      E->Jbwdfac1(B->face,store,tmp[0]);
    
      switch(n){
      case 0:
	if(!E->curvX)
	  daxpy (tot,B->nx.d,tmp[0],1,tmp[2],1);
	else
	  dvvtvp(tot,B->nx.p,1,tmp[0],1,tmp[2],1,tmp[2],1);
	break;
      case 1:
	if(!E->curvX)
	  daxpy (tot,B->ny.d,tmp[0],1,tmp[2],1);
	else
	  dvvtvp(tot,B->ny.p,1,tmp[0],1,tmp[2],1,tmp[2],1);
	break;
      case 2:
	if(!E->curvX)
	  daxpy (tot,B->nz.d,tmp[0],1,tmp[2],1);
	else
	  dvvtvp(tot,B->nz.p,1,tmp[0],1,tmp[2],1,tmp[2],1);
	break;
      }
    } 

    /* store existing flux */
    cnt = 0;
    for(i = 0; i < nfv; ++i)
      store[cnt++] = B->bvert[i];
    
    for(i = 0; i < nfv; ++i){
      dcopy((l=E->edge[E->ednum(B->face,i)].l),B->bedge[i],1,
	    store+cnt,1);
      cnt +=l;
    }
	
    if(lf) dcopy(lf,B->bface[0],1,store + cnt,1);
    
    E->MakeFlux(B,0,tmp[2]);
    
    /* add back original flux */
    cnt = 0;
    for(i = 0; i < nfv; ++i)
      B->bvert[i] += store[cnt++];
    
    for(i = 0; i < nfv; ++i){
      dvadd((l=E->edge[E->ednum(B->face,i)].l),B->bedge[i],1,
	    store+cnt,1,B->bedge[i],1);
      cnt += l;
    }
    
    if(lf) dvadd(lf,B->bface[0],1,store + cnt,1,B->bface[0],1);
    free(store); free_dmatrix(tmp,0,0);
  }  
}


static void AddDuDt(Bndry *B, Element_List *N[3], double dt){
  const    int id  = B->elmt->id;
  int      tot,nfv,cnt,lf,l,i,fq; 
  double  *store = dvector(0,LGmax*LGmax-1);
  double   **Q = dmatrix(0,4,0,QGmax*QGmax-1);
  Element  *Nu = N[0]->flist[id], *Nv = N[1]->flist[id], *Nw = N[2]->flist[id];
  int      qa = Nu->qa, qb = Nu->qb, qc = Nu->qc;
  double   **ima,**imb,*f,*fi;
  double  gamma = getgamma();

  fi = Q[3];

  Nu->GetFace(**Nu->h_3d, B->face, Q[1]);
  Nv->GetFace(**Nv->h_3d, B->face, Q[2]);
  Nw->GetFace(**Nw->h_3d, B->face, Q[3]);

  Nu->InterpToFace1(B->face,Q[1],Q[0]);
  Nv->InterpToFace1(B->face,Q[2],Q[1]);
  Nw->InterpToFace1(B->face,Q[3],Q[2]);

  nfv  = B->elmt->Nfverts(B->face);
  fq   = Nu->face[B->face].link->qface;

  if(nfv==3){
    switch(Nu->identify()){
    case Nek_Prism: case Nek_Pyr:
      tot = qa*qc;
      break;
    case Nek_Tet:
      tot = qa*qb;
      break;
    }
    getim   (fq,qa,&ima,g2a);
    getim   (fq,qb,&imb,g2b);
  }
  else{
    tot = Nu->qa*Nu->qb;
    getim(fq,qa,&ima,g2a);
    getim(fq,qb,&imb,g2a);
  }

  // Get u^n+1 from edge structure and interp to face 1
  f    = Nu->face[B->face].link->h;
  Interp2d(*ima,*imb, f,fq,fq,fi,qa,qb);
  dscal  (qa*qb,gamma/dt,fi,1);
  dsvtvm (qa*qb,1/dt,Q[0],1,fi,1,Q[0],1);

  f    = Nv->face[B->face].link->h;
  Interp2d(*ima,*imb, f,fq,fq,fi,qa,qb);
  dscal  (qa*qb,gamma/dt,fi,1);
  dsvtvm (qa*qb,1/dt,Q[1],1,fi,1,Q[1],1);

  f    = Nw->face[B->face].link->h;
  Interp2d(*ima,*imb, f,fq,fq,fi,qa,qb);
  dscal  (qa*qb,gamma/dt,fi,1);
  dsvtvm (qa*qb,1/dt,Q[2],1,fi,1,Q[2],1);

  if(Nu->curvX){
    dvmul  (tot, B->nx.p, 1, Q[0], 1, Q[3], 1);
    dvvtvp (tot, B->ny.p, 1, Q[1], 1, Q[3], 1, Q[3], 1);
    dvvtvp (tot, B->nz.p, 1, Q[2], 1, Q[3], 1, Q[3], 1);
  }
  else{
    dsmul  (tot, B->nx.d, Q[0], 1, Q[3], 1);
    dsvtvp (tot, B->ny.d, Q[1], 1, Q[3], 1, Q[3], 1);
    dsvtvp (tot, B->nz.d, Q[2], 1, Q[3], 1, Q[3], 1);
  }

  lf = B->elmt->face[B->face].l;
  if(nfv == 3)
    lf = lf*(lf+1)/2;
  else
    lf = lf*lf;

  /* store existing flux */
  cnt = 0;
  for(i = 0; i < nfv; ++i)
    store[cnt++] = B->bvert[i];
  
  for(i = 0; i < nfv; ++i){
    dcopy((l=Nu->edge[Nu->ednum(B->face,i)].l),B->bedge[i],1,
	  store+cnt,1);
    cnt +=l;
  }
	
  if(lf) dcopy(lf,B->bface[0],1,store + cnt,1);
    
  B->elmt->MakeFlux(B,0,Q[3]);
    
  /* add back original flux */
  cnt = 0;
  for(i = 0; i < nfv; ++i)
    B->bvert[i] += store[cnt++];
  
  for(i = 0; i < nfv; ++i){
    dvadd((l=Nu->edge[Nu->ednum(B->face,i)].l),B->bedge[i],1,
	  store+cnt,1,B->bedge[i],1);
    cnt += l;
  }
  
  if(lf) dvadd(lf,B->bface[0],1,store + cnt,1,B->bface[0],1);
  
  free(store);  free_dmatrix(Q,0,0); 
}

/*
 * PI = p - beta/2 U.n(U.n)_   (outflow boundary conditions)
 */
static void PbcOutflowStable(Bndry *B, Domain *omega){
  const    int id = B->elmt->id, face = B->face;
  Element  *U, *V, *W;
  double   *u = dvector(0,QGmax*QGmax-1), *v = dvector(0,QGmax*QGmax-1);
  double   *w = dvector(0,QGmax*QGmax-1), *f = dvector(0,QGmax*QGmax-1);
  double   *tmp = dvector(0,QGmax*QGmax-1);
	double		pval, fun;
  int       i, j, k, q, nfv, vn;
  static double beta = dparam("BACKFLOW_STABLE"); // stabilization coefficient
	if (beta == 0.0) beta = 0.2;

  U = omega->U->flist[id];
  V = omega->V->flist[id];
  W = omega->W->flist[id];

  nfv = U->Nfverts(face);

#ifdef PBC_1D
  /* use value computed through impedance B.C. */
  pval = omega->IMPBC[0].GetPval(B);
#endif

	// not implemented for curved elements yet
  for (i = 0; i < nfv; ++i){
    vn = U->vnum(face, i);
		fun = (U->vert[vn].hj[0]*B->nx.d) + 
					(V->vert[vn].hj[0]*B->ny.d) + (W->vert[vn].hj[0]*B->nz.d);
    B->bvert[i] = pval - 0.5 * beta * fun * (fun-fabs(fun));
  }

  if (nfv == 3){
    switch (U->identify()){
    case Nek_Prism: case Nek_Pyr:
      q = U->qa*U->qc;
      break;
    case Nek_Tet:
      q = U->qa*U->qb;
      break;
    }
  }
  else
    q = U->qa*U->qb;

  U->GetFace(**U->h_3d, face, tmp);
  U->InterpToFace1(face, tmp, u);

  V->GetFace(**V->h_3d, face, tmp);
  V->InterpToFace1(face, tmp, v);

  W->GetFace(**W->h_3d, face, tmp);
  W->InterpToFace1(face, tmp, w);

  //  tmp = u*nx+v*ny+w*nz
  if (U->curvX){
    dvmul  (q, B->nx.p, 1, u, 1, tmp, 1);
    dvvtvp (q, B->ny.p, 1, v, 1, tmp, 1, tmp, 1);
    dvvtvp (q, B->nz.p, 1, w, 1, tmp, 1, tmp, 1);
  }
  else{
    dsmul  (q, B->nx.d, u, 1, tmp, 1);
    dsvtvp (q, B->ny.d, v, 1, tmp, 1, tmp, 1);
    dsvtvp (q, B->nz.d, w, 1, tmp, 1, tmp, 1);
  }

  dvabs(q, tmp, 1, f, 1);
  dvneg(q, f, 1, f, 1);
  dvadd(q, tmp, 1, f, 1, f, 1);
  dscal(q, 0.5*beta, f, 1);
  dvmul(q, tmp, 1, f, 1, f, 1);
  dvneg(q, f, 1, f, 1);
  dsadd(q, pval, f, 1, f, 1);

	U->InterpToFace1(face, f, tmp);
  B->face = 0;
  U->JtransFace(B, tmp);
  B->face = face;

  free(u); free(v); free(w); free(f); free(tmp);

	return;
}
