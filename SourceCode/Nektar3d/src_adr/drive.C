/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $Source: 
 * $Revision:
 * $Date: 
 * $Author:
 * $State: 
 *---------------------------------------------------------------------------*/
// tcew done 
#include <stdio.h>
#include <time.h>
#include "nektar.h"
#include <cstdarg>

static int    Je;            /* Externals */
static double dt, Re;
#ifdef ADR
static int nspecs;
#endif

void MakeF   (Domain *omega, ACTION act);
static void StartUp(Domain *);
static void AddVarSource(Domain *);
static void PrintValues(Domain *);

void Bsystem_mem_free(Bsystem *Ubsys, Element_List *U);
void P_solve(Element_List *U, Element_List *Uf,Bndry *Ubc,Bsystem *Ubsys,
	     SolveType Stype,int step); 
void V_solve(Element_List *U, Element_List *Uf,Bndry *Ubc,Bsystem *Ubsys,
	     SolveType Stype, Domain *Omega, int step, ...);

#if DEBUG
#define Timing(s) \
 ROOTONLY  { fprintf(stdout,"%s Took %g seconds\n",s,(clock()-st)/cps); \
							       st = clock();}
#else
#define Timing(s) \
 /* Nothing */
#endif

int main(int argc, char *argv[]){
  Domain      *Omega;             /* Solution Domain            */
  int          step, nsteps;      /* Discrete time (step)       */
  double       time = 0.0;        /* Continuous time            */
  ACTION       WaveProp;          /* Convective form            */
  SLVTYPE      SolveType;         /* Solver type                */
  double       begin_clock,t,cfl,full_cfl,cfl_max,cfl_min;
  double       st, cps = (double)CLOCKS_PER_SEC;
	
  int          Torder, eid_max;
  
  init_comm(&argc, &argv);
  
#ifdef DEBUG
/*  mallopt(M_DEBUG,1); */
  init_debug();
#endif
  
  Omega     = PreProcess(argc,argv);
 
  time      = dparam("STARTIME");
  nsteps    = iparam("NSTEPS");
  WaveProp  = (ACTION) iparam("EQTYPE");
  SolveType = (SLVTYPE) iparam("SLVTYPE");
  Je        = iparam("INTYPE");
  dt        = dparam("DELT");
  step      = 0;
  Re        = 1./dparam("KINVIS");  
#ifdef ADR
  nspecs    = iparam("NSPEC");
#endif
  
  dparam_set("t",time);
 
	int	*wk = ivector(0,pllinfo.nprocs-1);
	int	Ndof[0];
	Ndof[0] = Omega->U->hjtot;
	gisum(Ndof,1,wk);
  ROOTONLY printf("rank No. %d: DOF is %d DOFTOT is %d\n", pllinfo.procid, Omega->U->hjtot, Ndof[0]);
  free(wk);
 
#ifndef WANNIER
  Omega->U->Terror(Omega->soln[0]); 
  Omega->V->Terror(Omega->soln[1]);
  Omega->W->Terror(Omega->soln[2]);
  Omega->P->Terror(Omega->soln[3]);
#endif

  Omega->U->Trans(Omega->U,J_to_Q);
  Omega->V->Trans(Omega->V,J_to_Q);
  Omega->W->Trans(Omega->W,J_to_Q);
  Omega->U->Set_state('p');
  Omega->V->Set_state('p');
  Omega->W->Set_state('p');
#ifdef ADR
	for (int i = 0; i < nspecs; i++){
  	Omega->T[i]->Trans(Omega->T[i],J_to_Q);	  
  	Omega->T[i]->Set_state('p');
	}
#endif

  StartUp (Omega);
  Omega->dt = dt;

  if (option ("SurForce") == 1)
   printf("Yes, I am dumping sectional forces\n");

  if(!step){ // impose boundary conditions 
    SetBCs(Omega->U,Omega->Ubc,Omega->Usys);
    SetBCs(Omega->V,Omega->Vbc,Omega->Vsys);
    SetBCs(Omega->W,Omega->Wbc,Omega->Wsys);
#ifdef ADR
		for (int i = 0; i < nspecs; i++)
			SetBCs(Omega->T[i],Omega->Tbc[i],Omega->Tsys[i]);
#endif  
  }

  forces(Omega,step,time);

  /*------------------------------------------------------------------*
   *                    Main Integration Loop                         *
   *------------------------------------------------------------------*/
  begin_clock = (double) clock();
  
  switch(SolveType){
  case Splitting:
    while (step < nsteps) {
      st = clock();
      MakeF (Omega, Prep);
      cfl = cfl_checker(Omega,dt);
      ROOTONLY fprintf(stdout,"CFL: %lf\n",cfl); 
      full_cfl = full_cfl_checker(Omega,dt,&eid_max);
      ROOTONLY fprintf(stdout,"Full CFL: %lf in eid %d\n",full_cfl,eid_max+1); 

      set_order(((step+1 < Je)? step+1: Je));

      Timing("Prep.......");
      
      MakeF (Omega, WaveProp);
      Timing("WavePropo..");
      
      // set v^{n+1} before pressure BC's determined 
      if(Omega->Tfun && Omega->Tfun->type == Womersley)
	SetWomSol(Omega,time+dt,1);
      
 if(option("tvarying") == 1){
	Bndry *Bc;
	t = dparam("t");
	dparam_set("t",t+dt);
	for(Bc=Omega->Ubc;Bc;Bc=Bc->next)
	  Bc->elmt->update_bndry(Bc,1);
	for(Bc=Omega->Vbc;Bc;Bc=Bc->next)
	  Bc->elmt->update_bndry(Bc,1);
	for(Bc=Omega->Wbc;Bc;Bc=Bc->next)
	  Bc->elmt->update_bndry(Bc,1);
 }
#ifdef ADR
 if(option("tvarying") != 0){
	Bndry *Bc;
	t = dparam("t");
	dparam_set("t",t+dt);
	for (int i = 0; i < nspecs; i++){
	  for(Bc=Omega->Tbc[i];Bc;Bc=Bc->next)
	    Bc->elmt->update_bndry(Bc,1);}
	}
#endif
      
      SetPBCs(Omega);
#ifdef ADR
      SetTBCs(Omega);
#endif

      Integrate_SS(Je, dt, Omega->U, Omega->Uf, Omega->u, Omega->uf);
      Integrate_SS(Je, dt, Omega->V, Omega->Vf, Omega->v, Omega->vf);
      Integrate_SS(Je, dt, Omega->W, Omega->Wf, Omega->w, Omega->wf);
#ifdef ADR
			for (int i = 0; i < nspecs; i++)
      	Integrate_SS(Je, dt, Omega->T[i], Omega->Tf[i], Omega->t[i], Omega->tf[i]);
#endif

      Timing("Integrate..");
      
      MakeF (Omega, Pressure);
      Timing("Pressure...");

      P_solve (Omega->P,Omega->Pf,Omega->Pbc,Omega->Pressure_sys,
	       Helm,step);
      Timing("P_solve....");
      
      MakeF (Omega, Viscous);
      Timing("Viscous....");
      
      V_solve (Omega->U,Omega->Uf,Omega->Ubc,Omega->Usys,Helm,Omega,step);
      V_solve (Omega->V,Omega->Vf,Omega->Vbc,Omega->Vsys,Helm,Omega,step);
      V_solve (Omega->W,Omega->Wf,Omega->Wbc,Omega->Wsys,Helm,Omega,step);
#ifdef ADR
			for (int i = 0; i < nspecs; i++)
    	  V_solve (Omega->T[i],Omega->Tf[i],Omega->Tbc[i],Omega->Tsys[i],Helm,Omega,step,i);
#endif
      
      Timing("V_solve....");
      
      MakeF (Omega, Post);
      Timing("Post.......");
      
      Analyser(Omega, ++step, (time += dt));
      Timing("Analyser...");
    }
    break;
  case SubStep:{	// Alireza: Not set up for ADR yet
    while (step < nsteps) {
      Omega->step = step;

      st = clock();

      MakeF (Omega, Prep);
      Timing("Prep.......");

      Torder = (step+1 < Je)? step+1: Je;
      set_order(Torder);
      

      SubCycle(Omega,Torder);
      Timing("SubCycle...");

      MakeF (Omega, Pressure);
      Timing("Pressure...");
      
      P_solve (Omega->P,Omega->Pf,Omega->Pbc,Omega->Pressure_sys,
	       Helm,step);
      Timing("P_solve....");
      
      MakeF (Omega, Viscous);
      Timing("Viscous....");
      
      V_solve (Omega->U,Omega->Uf,Omega->Ubc,Omega->Usys,Helm,Omega,step);
      V_solve (Omega->V,Omega->Vf,Omega->Vbc,Omega->Vsys,Helm,Omega,step);
      V_solve (Omega->W,Omega->Wf,Omega->Wbc,Omega->Wsys,Helm,Omega,step);
      
      Timing("V_solve....");
      
      MakeF (Omega, Post);
      Timing("Post.......");
      
      Analyser(Omega, ++step, (time += dt));
      Timing("Analyser...");
    }
  }
  break;
  }
  printf("User time of solver (seconds):  %lf \n",
	 (clock()-begin_clock)/(double)CLOCKS_PER_SEC);
  
#ifdef WOMERR
  if(Omega->Tfun && Omega->Tfun->type == Womersley)
    WomError(Omega,time);
#else
  Omega->U->Terror(Omega->soln[0]); 
  Omega->V->Terror(Omega->soln[1]);
  Omega->W->Terror(Omega->soln[2]);
  Omega->P->Terror(Omega->soln[3]);
#endif


  PostProcess(Omega, step, time);

  exit_comm();

  return 0;
}
  

//********************************************************
//                EXTRA ROUTINES                        //
//********************************************************

void P_solve(Element_List *U, Element_List *Uf, 
	   Bndry *Ubc, Bsystem *Ubsys, SolveType Stype, int ){

  SetBCs (U,Ubc,Ubsys);
  Solve  (U,Uf,Ubc,Ubsys,Stype);

}

void V_solve(Element_List *U, Element_List *Uf, 
     Bndry *Ubc, Bsystem *Ubsys, SolveType Stype, Domain *Omega, int step, ...){

	va_list(args);
	va_start(args, step);
	int spec = va_arg(args, int);

  if(step && step < Je){
    int k;
#ifdef ADR
    if(U->fhead->type == 't')
	    for(k=0;k<U->nel;++k)
  	    Ubsys->lambda[k].d = Re*Omega->Pr[spec]*getgamma()/dt;
		else
#endif
    	for(k=0;k<U->nel;++k)
      	Ubsys->lambda[k].d = Re*getgamma()/dt;

    switch(U->fhead->type){
		case 'u': case 't':
      Bsystem_mem_free(Ubsys, U);
      goto ReCalcMat;
    case 'v':
      if(option("REFLECT1")||option("REFLECT0")){
	if(Ubsys->smeth == direct) 
	  Bsystem_mem_free(Ubsys, U);
	else
	  Ubsys->Gmat = Omega->Usys->Gmat;
	goto ReCalcMat;
      }
      else 
	break;
    case 'w':
      if(option("REFLECT2")||
	 (option("REFLECT0")&&option("REFLECT1")&&(!option("REFLECT2")))){
	if(Ubsys->smeth == direct) 
	  Bsystem_mem_free(Ubsys, U);
	else
	  Ubsys->Gmat = Omega->Usys->Gmat;
	goto ReCalcMat;
      }
      else 
	break;
			
    ReCalcMat:
      if(Ubsys->Precon == Pre_LEnergy)
			option_set("ReCalcPrecon",0); // don't recalc precon at present
      double *save_hj = dvector(0, U->hjtot-1);
      ROOTONLY fprintf(stdout,"Regenerating %c-Matrix\n",U->fhead->type);
      dcopy(U->hjtot, U->base_hj, 1, save_hj, 1);
      GenMat(U,Ubc,Ubsys,Ubsys->lambda,Helm);
      dcopy(U->hjtot, save_hj, 1, U->base_hj, 1);
      free(save_hj);
      break;
    }

    if(Ubc&&Ubc->DirRHS){
      free(Ubc->DirRHS);
      Ubc->DirRHS = (double*) 0;
    }

#ifndef PARALLEL
    if(!option("tvarying"))
      DirBCs(U,Ubc,Ubsys,Helm);
#endif
  }

  SetBCs (U,Ubc,Ubsys);
  Solve  (U,Uf,Ubc,Ubsys,Stype);

	va_end(args);
}

void MakeF(Domain *omega, ACTION act){

  Element_List  *U    =  omega->U,  *V    =  omega->V,   *W   = omega->W,
                *Uf   =  omega->Uf, *Vf   =  omega->Vf,  *Wf  = omega->Wf,
                *P    =  omega->P,  *Pf   =  omega->Pf;
#ifdef ADR
	Element_List	**T    =  omega->T,  **Tf   =  omega->Tf;
#endif

  int Nmodes = U->hjtot, Nquad = U->htot;

  switch (act) {
  case Prep: /* put fields in physical space for waveprop */
    U->Trans(U,J_to_Q);
    V->Trans(V,J_to_Q);
    W->Trans(W,J_to_Q);

    U->Set_state('p');
    V->Set_state('p');
    W->Set_state('p');
    
    dcopy(U->htot,U->base_h,1,omega->u[0],1);
    dcopy(V->htot,V->base_h,1,omega->v[0],1);
    dcopy(W->htot,W->base_h,1,omega->w[0],1);
#ifdef ADR
		for (int i = 0; i < nspecs; i++){
	    T[i]->Trans(T[i],J_to_Q);	  
  	  T[i]->Set_state('p');
  	  dcopy(T[i]->htot,T[i]->base_h,1,omega->t[i][0],1);
		}
#endif

    break;
    
  case Rotational:
    VxOmega (omega);
    goto AddForcing; 

  case Convective:
    VdgradV (omega); 
    goto AddForcing;

  case Stokes:
    //    StokesBC (omega);
    Uf->zerofield();
    Vf->zerofield();
    Wf->zerofield();
    Uf->Set_state('p');
    Vf->Set_state('p');
    Wf->Set_state('p');
    goto AddForcing; 

  AddForcing: {
    double fx = dparam("FFX");
    double fy = dparam("FFY");
    double fz = dparam("FFZ");

    if(fx)
      dsadd(Nquad, fx, Uf->base_h, 1, Uf->base_h, 1);
    if(fy)
      dsadd(Nquad, fy, Vf->base_h, 1, Vf->base_h, 1);
    if(fz)
      dsadd(Nquad, fz, Wf->base_h, 1, Wf->base_h, 1);

    if(omega->ForceFuncs){
      dvadd(Nquad, omega->ForceFuncs[0], 1, Uf->base_h, 1, Uf->base_h, 1);
      dvadd(Nquad, omega->ForceFuncs[1], 1, Vf->base_h, 1, Vf->base_h, 1);
      dvadd(Nquad, omega->ForceFuncs[2], 1, Wf->base_h, 1, Wf->base_h, 1);
    }

		// Alireza: bodyforces + variable heat/mass source terms
    AddVarSource( omega );

#ifdef ADR
    double gx = dparam("GX");		// Alireza: buoyancy-driven flow and natural convection
    double gy = dparam("GY");
    double gz = dparam("GZ");
		double Ra = dparam("Ra");
		for (int i = 0; i < nspecs; i++){
			double RaPr = Ra * omega->Pr[i];
	    if(gx)
		  	dsvtvp (Nquad, RaPr, T[i]->base_h, 1, Uf->base_h, 1, Uf->base_h, 1);
			if(gy)
	   		dsvtvp (Nquad, RaPr, T[i]->base_h, 1, Vf->base_h, 1, Vf->base_h, 1);
			if(gz)
	   		dsvtvp (Nquad, RaPr, T[i]->base_h, 1, Wf->base_h, 1, Wf->base_h, 1);
		}

#if defined (CCHF) || (CCLF)
		// Alireza: coagulation cascade
		AddCCSource( omega );
#ifdef CCLF
		// Alireza: Brinkman friction term for a porous region
		AddBrinkman( omega );
#endif
#endif
#endif

    break;
  }

  case Pressure: { 
    double    dtinv = 1./dt;
    double *nul= (double*)0;
    
      
    U->Grad_d(Pf->base_h, nul, nul, 'x');

    V->Grad_d(nul, P->base_h, nul, 'y');
    dvadd(Nquad, P->base_h, 1, Pf->base_h, 1, Pf->base_h, 1);

    W->Grad_d(nul,  nul,  P->base_h, 'z');
    dvadd(Nquad, P->base_h, 1, Pf->base_h, 1, Pf->base_h, 1);

    Pf->Set_state('p');

    Pf->Iprod(Pf);

    dscal(Nmodes, dtinv, Pf->base_hj, 1);
    Pf->Set_state('t');

    break;
  }
    
  case Viscous: {
    double dtinv = 1/dt;

    P->Trans(P, J_to_Q); 
    P->Set_state('p');
    P->Grad_d (Uf->base_h,Vf->base_h,Wf->base_h,'a');
    P->Set_state('t');

    daxpy(Nquad, -dtinv, U->base_h, 1, Uf->base_h, 1);
    Uf->Set_state('p');
    Uf->Iprod(Uf);    
    dscal(Nmodes,    Re, Uf->base_hj, 1);
    Uf->Set_state('t');

    daxpy(Nquad, -dtinv, V->base_h, 1, Vf->base_h, 1);
    Vf->Set_state('p'); 
    Vf->Iprod(Vf);
    dscal(Nmodes,    Re, Vf->base_hj, 1);
    Vf->Set_state('t');

    daxpy(Nquad, -dtinv, W->base_h, 1, Wf->base_h, 1);
    Wf->Set_state('p');
    Wf->Iprod(Wf);
    dscal(Nmodes,    Re, Wf->base_hj, 1);
    Wf->Set_state('t');

    if (option("SVVE") && iparam("NCutM")){ //Alireza: explicit SVV filtering
        SVVExplicit(U, Uf);
        SVVExplicit(V, Vf);
        SVVExplicit(W, Wf);
		}

#ifdef ADR
		for (int i = 0; i < nspecs; i++){
			dcopy(Nquad, T[i]->base_h, 1, Tf[i]->base_h, 1);
			dscal(Nquad, -dtinv, Tf[i]->base_h, 1);  
			Tf[i]->Set_state('p');
			Tf[i]->Iprod(Tf[i]);
			double RePr = Re * omega->Pr[i];
			dscal(Nmodes,    RePr, Tf[i]->base_hj, 1);
			Tf[i]->Set_state('t');
			
			if (option("SVVE") && iparam("NCutM"))
				SVVExplicit(T[i], Tf[i]);	//Alireza: explicit SVV filtering
		}
#endif
		
    break;
  }
  case Post: {

    break;
    
  }
  default:
    error_msg(MakeF--unknown type of action);
    break;
  }

  return;
}

/* Do initial time step assuming for case where startup field is in physical *
 * space or copy V field to Vs if it is a restart                           */

static void StartUp(Domain *omega){

  if(omega->U->fhead->state == 'p'){
    ROOTONLY fprintf(stdout,"Locally transforming initial conditions\n");
    omega->U->Trans(omega->U,Q_to_J);
    omega->V->Trans(omega->V,Q_to_J);
    omega->W->Trans(omega->W,Q_to_J);  
    omega->U->Set_state('t');
    omega->V->Set_state('t');
    omega->W->Set_state('t');
    omega->P->Set_state('t');
#ifdef ADR
		for (int i = 0; i < nspecs; i++){
			omega->T[i]->Trans(omega->T[i],Q_to_J);
			omega->T[i]->Set_state('t');
		}
#endif
  }
}

static void AddVarSource(Domain *omega){
	int nfields;
  Coord X;
  X.x = dvector(0, QGmax*QGmax*QGmax - 1);
  X.y = dvector(0, QGmax*QGmax*QGmax - 1);
  X.z = dvector(0, QGmax*QGmax*QGmax - 1);
  double *func = dvector(0, QGmax*QGmax*QGmax - 1);
#ifdef ADR
	nfields = 3 + nspecs;
#else
	nfields = 3;
#endif
	Element **E = (Element **) malloc(nfields*sizeof(Element *));

	for (int i = 0; i < omega->U->nel; ++i){
  	E[0] = omega->Uf->flist[i]; E[1] = omega->Vf->flist[i]; E[2] = omega->Wf->flist[i];
#ifdef ADR
		for (int k = 0; k < nspecs; k++)
			E[3 + k] = omega->Tf[k]->flist[i];
#endif
		E[0]->coord(&X);
		for (int j = 0; j < nfields; ++j){
			if (strcmp(omega->sources[j],"null") == 0) continue;
			vector_def("x y z", omega->sources[j]); // evaluate each source
		  vector_set(E[j]->qtot, X.x, X.y, X.z, func); // source values in func vector
			dvadd(E[j]->qtot, func, 1, E[j]->h_3d[0][0], 1, E[j]->h_3d[0][0], 1);
		}
	}

  free(X.x); free(X.y); free(X.z);
  free(func);
}

double setQfactor(int ind, int limit, int low, int high){ //SVV
 double num   = ind - (limit - high),
        denom = ind - low + 1 ,
        Q;

 if(ind>=limit-high)
   num = 0.0;

 Q = num*num/(denom*denom);
 return Q;
}

void SVVExplicit(Element_List *U, Element_List *Uf){
	int   Nm= U->fhead->Nmodes,qt= U->fhead->qtot;
  Basis *b,*db;

  static double *store,*u1,*u2,*u3;
  static double *str;
  static double *ftmp1;
  static double *ftmp2;

  int i,j,k,cnt;
  double Qa=0,Qb=0,Qc=0;
  static double *filter;
  static double  Epsilon;
  static int NCutMT, NCutM, init=0;

  if(init==0){
  	store=dvector(0,4*QGmax*QGmax*QGmax-1);
    u1    = store+QGmax*QGmax*QGmax;
    u2    = store+2*QGmax*QGmax*QGmax;
    u3    = store+3*QGmax*QGmax*QGmax;
    str = dvector(0,3*qt-1);
    ftmp1 = dvector(0,QGmax*QGmax*QGmax-1);
    ftmp2 = dvector(0,QGmax*QGmax*QGmax-1);
    filter = dvector(0,Nm-1);

    init_ortho_basis();
    Epsilon = dparam("EPSILON");
    double SVVEfactor = dparam("SVVEfactor");
    if (SVVEfactor > 1.0) SVVEfactor = 1.0;

    NCutMT = iparam("MNT");
    NCutM = (iparam("NCutM")?iparam("NCutM"):((int) ((U->fhead->lmax-2)*(1.0-SVVEfactor))));
    ROOTONLY
    	fprintf(stdout,"SVVE:: in viscous step: NCutM = %d, NCutMT = %d, EPSILON = %f SVVEfactor = %f\n",
                     NCutM, NCutMT, Epsilon, SVVEfactor);

    init = 1;
    // set up filter vector to scale weak Laplacian
    dzero(Nm,filter,1);
    cnt = 0;

    if(U->fhead->identify() == Nek_Tet)  {
      for (i = 0; i < U->fhead->lmax; i++){
        for (j = 0; j < U->fhead->lmax-i; j++){
          for (k = 0; k<U->fhead->lmax-i-j;k++,++cnt){
            if(i+j+k >= NCutM)  {
              Qb = setQfactor (i+j+k, U->fhead->lmax, NCutM, NCutMT);
              filter[cnt] = exp(-Qb);
            }
          }
        }
      }
    }
    else if(U->fhead->identify() == Nek_Hex)  {
      for (i = 0; i < U->fhead->lmax; i++){
        Qa=0;
        if( i >= NCutM)
          Qa = setQfactor ( i, U->fhead->lmax, NCutM, NCutMT);
          for (j = 0; j < U->fhead->lmax; j++){
         		Qb=0;
           	if( j>= NCutM)
           		Qb = setQfactor ( j, U->fhead->lmax, NCutM, NCutMT);
             	for (k=0; k< U->fhead->lmax;k++,++cnt){
               	Qc=0;
               	if( k>= NCutM)
               		Qc = setQfactor ( k, U->fhead->lmax, NCutM, NCutMT);
                 	if( i >= NCutM || j>= NCutM || k>= NCutM)  {
                  	Qb= (i>j?Qa:Qb);
                   	Qb= (k>j?Qc:Qb);
                   	filter[cnt] = exp(-(Qb));
                 	}
            	}
          }
      }
    }
  }//if(init)

  //back to u^n
  U->Trans(U,J_to_Q);
  for(i=0;i<U->nel;i++)  {
   	dzero(4*QGmax*QGmax*QGmax, store, 1);
    b  = U->flist[i]->getbasis();
    db = U->flist[i]->derbasis();

    U->flist[i]->Grad_d(u1, u2, u3, 'a');

    dcopy(qt,u1,1,ftmp1,1);
    U->flist[i]->Ofwd(ftmp1,ftmp2,U->fhead->lmax); //physical->modal
    dvmul(Nm,filter,1,ftmp2,1,ftmp2,1);                  //apply filter
    U->flist[i]->Obwd(ftmp2,ftmp1,U->fhead->lmax); //modal->physical
    dcopy(qt,ftmp1,1,str,1);

    dcopy(qt,u2,1,ftmp1,1);
    U->flist[i]->Ofwd(ftmp1,ftmp2,U->fhead->lmax);
    dvmul(Nm,filter,1,ftmp2,1,ftmp2,1);
    U->flist[i]->Obwd(ftmp2,ftmp1,U->fhead->lmax);
    dcopy(qt,ftmp1,1,str+qt,1);

    dcopy(qt,u3,1,ftmp1,1);
    U->flist[i]->Ofwd(ftmp1,ftmp2,U->fhead->lmax);
    dvmul(Nm,filter,1,ftmp2,1,ftmp2,1);
    U->flist[i]->Obwd(ftmp2,ftmp1,U->fhead->lmax);
    dcopy(qt,ftmp1,1,str+2*qt,1);

    U->flist[i]->form_diprod(str,str+qt,str+2*qt,b->vert+1);

    dsmul(qt,Epsilon,str,1,U->flist[i]->h_3d[0][0],1);
    U->flist[i]->Iprod_3d(U->flist[i],db,b,b);
    dvadd(Nm,U->flist[i]->vert->hj,1,store,1,store,1);

    dsmul(qt,Epsilon,str+qt,1,U->flist[i]->h_3d[0][0],1);
    U->flist[i]->Iprod_3d(U->flist[i],b,db,b);
    dvadd(Nm,U->flist[i]->vert->hj,1,store,1,store,1);

    dsmul(qt,Epsilon,str+2*qt,1,U->flist[i]->h_3d[0][0],1);
    U->flist[i]->Iprod_3d(U->flist[i],b,b,db);
    dvadd(Nm,U->flist[i]->vert->hj,1,store,1,store,1);

    dcopy(Nm,store,1,U->flist[i]->vert->hj,1);
    U->flist[i]->state = 't';
	}//end of for(i=0;i<U->nel;i++)

	dvsub(U->hjtot,Uf->base_hj,1,U->base_hj,1,Uf->base_hj,1);// for Je>1, needs to sum_{q=0}^{q<Je} *^{n-1}
}

static void PrintValues(Domain *omega){
  FILE *pFile;
  char fname[FILENAME_MAX];
  int i,j;
  int face;
  int      ll, nfv;
  Element  *E;
  Bndry    *B;

  sprintf(fname, "bvert_%d.dat", pllinfo.procid);
  pFile = fopen(fname, "w");
#ifdef ADR
  for (i = 0; i < nspecs; i++){
    B = omega->Tbc[i];
    while (B){
      face = B->face;
      E = B->elmt;
      nfv = E->Nfverts(face);
      for(j = 0; j < nfv; ++j)
        fprintf(pFile, "%d %d %f \n", E->id, face, B->bvert[j]);
      B = B->next;
    }
  }
#endif
  fclose(pFile);
}
