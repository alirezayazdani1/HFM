/*
 * Womersley Flow 
 *
 * This file sets up the solution for Womersley type flow, an exact
 * solution to the Navier-Stokes equations for a osscilatory flow
 * 
 * ------------------------------------------------------------------------- */
#include <stdio.h>
#include <ctype.h>
#include <stdarg.h>
#if defined (_CRAY) // this definition must be before nektar.h
#define zbesj_ ZBESJ
#endif

#if defined (_AIX) // this definition must be before nektar.h
#define zbesj_ zbesj
#endif
#include "nektar.h"

#define SQRT2 1.41421356237309504880

#if DIM == 3
extern   int vnum[][3], ednum[][3];
#else
extern   int vnum[][2];
#endif

#define TANTOL 1e-10
/* new atan2 function to stop Nan on atan(0,0)*/
static double atan2_proof (double x, double y)
{
  if (fabs(x) + fabs(y) > TANTOL) return (atan2(x,y));
  else return (0.);
}
#define atan2 atan2_proof

typedef struct {
  double re;
  double im;
} Cmplx;

#ifdef WOMERSLEY
static void velocity(int k, int qt, double ccos, double csin, double *r, 
         double R, double wnum, double *re, double *im);


static char     dir;
#if DIM == 2
static void womersley(int n, double *x, double *y, double *u, char dir);
#else
static void womersley(int n, double *x, double *y, double *z,
		      double *u, char dir);
#endif
/* 
 * Redefine vec_init() and set_vector() 
 * ---------------------------------------------------------------------- */
 
void vector_def(char *s, char *f)
{
  while (isspace(*f++));
  dir = *--f;
  return;
}

void vector_set(int n, ...)
{
#if DIM == 2
  double  *x, *y, *f, t;
#else
  double  *x, *y, *z, *f, t;
#endif
  va_list ap;

  va_start(ap, n);
  x = va_arg(ap, double *);
  y = va_arg(ap, double *);
#if DIM == 3
  z = va_arg(ap, double *);
#endif
  f = va_arg(ap, double *);
  va_end  (ap);

#if DIM == 2
  womersley(n, x, y, f, dir);
#else
  womersley(n, x, y, z, f, dir);
#endif
  return;
}

/* ---------------------------------------------------------------------- */

#if DIM == 2
static void womersley(int n, double *x, double *y, double *u, char dir)
#else
static void womersley(int n, double *x, double *y, double *z, double *u, 
		      char dir)
#endif
{
  register int i;

  switch (dir) {
  case 'u':
    {
#if DIM == 2
      double a,w,h,kinvis,t;
      double delta,k,fac;
      double Rev,Imv;

      a      = dparam("WOMA");
      w      = dparam("WOMW");
      kinvis = dparam("KINVIS");
      h      = dparam("WOMHEIGHT");
      t      = dparam("TIME");

      delta  = sqrt(2*kinvis/fabs(w));
      k      = 1/delta;

      fac = 1.0/(cos(0.5*h*k)*cosh(0.5*h*k)*cos(0.5*h*k)*cosh(0.5*h*k)
		 + sin(0.5*h*k)*sinh(0.5*h*k)*sin(0.5*h*k)*sinh(0.5*h*k));
      for (i = 0; i < n; i++)
	u[i]  =  a/w*cos(w*t)*(1-(cos(0.5*h*k)*cos(y[i]*k)*cosh(0.5*h*k)*
			       cosh(y[i]*k) + sin(0.5*h*k)*sin(y[i]*k)*
			       sinh(0.5*h*k)*sinh(y[i]*k))*fac) - 
	         a/w*sin(w*t)*(cos(y[i]*k)*cosh(y[i]*k)*sin(0.5*h*k)*
			       sinh(0.5*h*k) - cos(0.5*h*k)*cosh(0.5*h*k)*
			       sin(y[i]*k)*sinh(y[i]*k))*fac;
#else
    for (i = 0; i < n; i++) 
      u[i]  = 0.0;
#endif
    }
    break;
  case 'v':
    for (i = 0; i < n; i++) 
      u[i]  = 0.0;
    break;
#if DIM == 3
  case 'w':
    {
      double a,w,h,kinvis,t;
      double rad,re,im;
      Cmplx pulse[2];

      a      = dparam("WOMA");
      w      = dparam("WOMW");
      kinvis = dparam("KINVIS");
      h      = dparam("WOMHEIGHT");
      t      = dparam("TIME");

      pulse[0].re = 0.;
      pulse[0].im = 0.;
      pulse[1].re = a;
      pulse[1].im = 0.;

      for (i=0;i<n;++i) {
          rad = sqrt(y[i]*y[i] + x[i]*x[i]);
          velocity(0, 1,pulse[0].re,pulse[0].im ,&rad, 0.5,w,&re,&im);
          u[i] = re;
          velocity(1, 1,pulse[1].re,pulse[1].im ,&rad, 0.5,w,&re,&im);
          u[i] += re*cos(w*t) -  im*sin(w*t);
      }

      rad=0.;
      velocity(0, 1,pulse[0].re,pulse[0].im ,&rad, 0.5,w,&re,&im);
      h = re;u[i] = re;
      velocity(1, 1,pulse[1].re,pulse[1].im ,&rad, 0.5,w,&re,&im);
      h += re*cos(w*t) -  im*sin(w*t);
      printf("%lf \n",h);

    }
    break;
#endif
  case 'p':
    for (i = 0; i < n; i++) 
      u[i] = 0.0;
    break;
  case '0':
    for (i = 0; i < n; i++) 
      u[i] = 0.0;
    break;
  default:
    printf("unknown direction in womersley -- %c\n", dir);
    break;
  }
  
  return;
}

static void velocity(int k, int qt, double ccos, double csin, double *r, 
         double R, double wnum, double *re, double *im){
   int i;
   double alpha,mu=dparam("KINVIS"),rho=1.0;
   Cmplx jc1,jc2,a,b,c;
   int  KODE=1,NZ, IERR,N=1;
   double FNU=0.;


   if (k) {
      alpha = R*sqrt(wnum/mu);
      jc2.re = -alpha/SQRT2;
      jc2.im = alpha/SQRT2;
      for (i=0;i<qt;++i) {

         jc1.re = jc2.re*r[i];
         jc1.im = jc2.im*r[i];

         zbesj_(&jc1.re, &jc1.im, &FNU, &KODE, &N, &c.re, &c.im, &NZ, &IERR);
         if (IERR) {printf("A %d error\n",IERR);exit(-1);}
	 
	 


         zbesj_(&jc2.re, &jc2.im, &FNU, &KODE, &N, &b.re, &b.im, &NZ, &IERR);
         if (IERR) {printf("B %d error\n",IERR);exit(-1);}
         a.re = b.re - c.re;
         a.im = b.im - c.im;
         c.re = csin / (rho * wnum);
         c.im = ccos / (rho * wnum);

         re[i] = (a.re*b.re*c.re + a.im*b.im*c.re - a.im*b.re*c.im + a.re*b.im*c.im) /
             (b.re*b.re + b.im*b.im);
         im[i] = (a.re*b.re*c.im + a.im*b.im*c.im + a.im*b.re*c.re - a.re*b.im*c.re) /
             (b.re*b.re + b.im*b.im);
      }
   } else {
      for (i=0;i<qt;++i) {
         re[i] = ccos*R*R*(r[i]*r[i]-1.)/(4.*mu);
         im[i] = 0.;
      }
   }
}


void WomInit(Domain *Omega){
}

void SetWomSol(Domain *Omega, double time, int save){
}

#else


static void   velocity(int k, int qt,  double ccos, double csin, double *r, 
		       double R, double wnum, double *cw, double *sw);
static void   massflow(int k, double *ccos, double *csin, double R,
			double wnum);
static void   encode(int n, double *f, int m, double *ccos, double *csin, double T);
static double decode(int n,  double *ccos, double *csin, double t, double T);

typedef struct wominfo {
  int nbcs;
  int *bcid;
  double ***U_re;
  double ***U_im;
  double ***V_re;
  double ***V_im;
  double ***W_re;
  double ***W_im;
} Wominfo;
Wominfo Winfo;


static void WomMem(Domain *omega){
  Bndry *Ubc;
  int nbcs,nmodes,Lsize,l;

  Lsize = 0;
  for(Ubc = omega->Ubc,nbcs=0;Ubc;Ubc= Ubc->next){
    switch(Ubc->type){
    case 'V': case 'v': case 'I': case 'i':
      nbcs++;
      l = Ubc->elmt->lmax;
      if(Ubc->elmt->Nfverts(Ubc->face) == 3)
	Lsize += l*(l+1)/2;
      else
	Lsize += l*l;
      break;
    }
  }
   
  if(!nbcs) return;

  nmodes = omega->Tfun->info.wexpand.nmodes;
  
  Winfo.nbcs   = nbcs;
  Winfo.bcid   = ivector(0,nbcs-1);
  Winfo.U_re   = dtarray(0,nbcs-1,0,nmodes-1,0,Lsize-1);
  Winfo.U_im   = dtarray(0,nbcs-1,0,nmodes-1,0,Lsize-1);
  Winfo.V_re   = dtarray(0,nbcs-1,0,nmodes-1,0,Lsize-1);
  Winfo.V_im   = dtarray(0,nbcs-1,0,nmodes-1,0,Lsize-1);
  Winfo.W_re   = dtarray(0,nbcs-1,0,nmodes-1,0,Lsize-1);
  Winfo.W_im   = dtarray(0,nbcs-1,0,nmodes-1,0,Lsize-1);

  for(Ubc = omega->Ubc,nbcs=0;Ubc;Ubc= Ubc->next){
    switch(Ubc->type){
        case 'V': case 'v': case 'I': case 'i':
	  Winfo.bcid[nbcs++] = Ubc->id;
	  break;
	}
  }  
}
   
void WomInit(Domain *omega){
  double R, *axs, nrm[3], *r;
  double phi, cp, sp, theta, ct, st, *x, *y, *z;
  double wnum, *ccos, *csin, *f, ccos0;
  int i, j,l,qt,nbc,eid,face,cnt,nfverts;
  Coord   X;
  Element *E;
  Bndry  *Ubc,*Vbc;
  
  option_set("tvarying",2);

  R    = omega->Tfun->info.wexpand.radius;
  axs  = omega->Tfun->info.wexpand.axispt;
  wnum = omega->Tfun->info.wexpand.wnum;
  ccos = omega->Tfun->info.wexpand.ccos;
  csin = omega->Tfun->info.wexpand.csin;
  f    = omega->Tfun->info.wexpand.raw;
  
  // leave nrm untouched for recalling function
  dcopy(3,omega->Tfun->info.wexpand.axisnm,1,nrm,1);
  
  /* deal with encoding and mass flow conversions */
  if(omega->Tfun->info.wexpand.type == RAW) {
    omega->Tfun->info.wexpand.nmodes = iparam("NFOU");
    ccos = dvector(0,omega->Tfun->info.wexpand.nmodes-1);
    csin = dvector(0,omega->Tfun->info.wexpand.nmodes-1);
    encode(omega->Tfun->info.wexpand.nraw,f,
	   omega->Tfun->info.wexpand.nmodes,ccos,csin,2*M_PI/wnum);
  }
  
  if (omega->Tfun->info.wexpand.form == MASSFLOW)
    for (i=0;i<omega->Tfun->info.wexpand.nmodes;++i)
      massflow( i, &ccos[i], &csin[i], R, i*wnum);
  
  WomMem(omega);
  
  X.x = dvector(0,max(QGmax*(QGmax+1),LGmax*LGmax-1));
  X.y = dvector(0,max(QGmax*(QGmax+1),LGmax*LGmax-1));
  X.z = dvector(0,max(QGmax*(QGmax+1),LGmax*LGmax-1));
  r   = dvector(0,max(QGmax*(QGmax+1),LGmax*LGmax-1));
  
  phi  = atan2(nrm[1], nrm[0]);
  cp  = cos(phi); sp = sin(phi);
  drot(1,&nrm[0],1,&nrm[1],1,cp,sp);
  theta = atan2(nrm[0],nrm[2]);
  ct    = cos(theta); st = sin(theta);
  
  /* find coordinates of all quadrature points on region */
  for(Ubc = omega->Ubc,nbc=0;Ubc;Ubc= Ubc->next){
    switch(Ubc->type){
    case 'V': case 'v': case 'I': case 'i':
      eid  = Ubc->elmt->id;
      face = Ubc->face;

      Vbc     = omega->Vbc+Ubc->id;
      E       = omega->U->flist[eid];
      nfverts = E->Nfverts(face);
      
      if(E->identify() == Nek_Tet){
	if(face == 0)
	  qt = E->qa*E->qb;
	else 
	  if(face == 1)
	    qt = E->qa*E->qc;
	  else 
	    qt = (E->qb+1)*E->qc;
      }
      else{
	if(face == 0)
	  qt = E->qa*E->qb;
	else if ((face == 1)||(face == 3))
	  qt = E->qa*E->qc;
	else
	  qt = E->qb*E->qc;
      }
      
      E->GetFaceCoord(face,&X);

      X.x[qt] = E->vert[E->vnum(face,2)].x;
      X.y[qt] = E->vert[E->vnum(face,2)].y;
      X.z[qt] = E->vert[E->vnum(face,2)].z;
      
      if(nfverts == 4){
	X.x[qt+1] = E->vert[E->vnum(face,3)].x;
	X.y[qt+1] = E->vert[E->vnum(face,3)].y;
	X.z[qt+1] = E->vert[E->vnum(face,3)].z;
      }
      else{
	X.x[qt+1] = 0.0;
	X.y[qt+1] = 0.0;
	X.z[qt+1] = 0.0;
      }
      
      /* subtract offset */
      dsadd(qt+2,-axs[0],X.x,1,X.x,1);
      dsadd(qt+2,-axs[1],X.y,1,X.y,1);
      dsadd(qt+2,-axs[2],X.z,1,X.z,1);

      /* rotate to align normal with z axis */
      drot(qt+2,X.x,1,X.y,1,cp,sp);
      drot(qt+2,X.z,1,X.x,1,ct,st);
      
      /* calculate scaled radius */
      dvmul (qt+2,X.x,1,X.x,1,r,1);
      dvvtvp(qt+2,X.y,1,X.y,1,r,1,r,1);
      dvsqrt(qt+2,r,1,r,1);
      dscal (qt+2,1./R,r,1);
      
      for (i = 0; i < omega->Tfun->info.wexpand.nmodes; ++i){
	/* find velocity in terms of real and imaginary parts 
	   for a single fourier mode for every quadrature point */
	velocity(i, qt+2, ccos[i], csin[i], r, R, i*wnum, X.x, X.y);
	  

	for (j=0; j<qt+2; ++j) {
	  X.x[j] *= omega->Tfun->info.wexpand.scal0;
	  X.y[j] *= omega->Tfun->info.wexpand.scal0;
	}
	
	
	if(nfverts == 3){
	  Ubc->bvert[2] = X.x[qt];
	  Vbc->bvert[2] = X.y[qt];

	  E->JtransFace(Ubc,X.x);	
	  E->JtransFace(Vbc,X.y);

	}
	else{
	  // set vertices 
	  Ubc->bvert[0] = X.x[0];
	  Ubc->bvert[1] = X.x[E->qa-1];
	  Ubc->bvert[2] = X.x[qt];
	  Ubc->bvert[3] = X.x[qt+1];

	  E->InterpToFace1(face,X.x,X.z);
	  E->JtransFace(Ubc,X.z);	

	  Vbc->bvert[0] = X.y[0];
	  Vbc->bvert[1] = X.y[E->qa-1];
	  Vbc->bvert[2] = X.y[qt];
	  Vbc->bvert[3] = X.y[qt+1];

	  E->InterpToFace1(face,X.y,X.z);
	  E->JtransFace(Vbc,X.z);
	}
	
	/* copy real transformed values into X.z such that values
	   are stored contiguously in X.z rather than uncontiguously
	   in Ubc */
	dcopy(nfverts,Ubc->bvert,1,X.z,1);
	cnt = nfverts;
	for(j = 0; j < nfverts; ++j){
	  l = E->edge[E->ednum(face,j)].l;
	  dcopy(l,Ubc->bedge[j],1,X.z+cnt,1);
	  cnt += l;
	}
	
	l = E->face[face].l;
	if(nfverts == 3)
	  l = l*(l+1)/2;
	else
	  l = l*l;

	if (l)
	  dcopy(l,Ubc->bface[0],1,X.z+cnt,1);
	cnt += l;
	
	dzero(cnt,X.x,1);
	dzero(cnt,X.y,1);
	drot(cnt,X.z,1,X.x,1,ct,-st);
	drot(cnt,X.x,1,X.y,1,cp,-sp);
	
	dcopy(cnt,X.x,1,Winfo.U_re[nbc][i],1);
	dcopy(cnt,X.y,1,Winfo.V_re[nbc][i],1);
	dcopy(cnt,X.z,1,Winfo.W_re[nbc][i],1);
	
	/* copy imaginary transformed values into X.z */
	dcopy(nfverts,Vbc->bvert,1,X.z,1);
	cnt = nfverts;
	for(j = 0; j < nfverts; ++j){
	  l = E->edge[E->ednum(face,j)].l;
	  dcopy(l,Vbc->bedge[j],1,X.z+cnt,1);
	  cnt +=l;
	}

	l = E->face[face].l;
	if(nfverts == 3)
	  l = l*(l+1)/2;
	else
	  l = l*l;

	if (l)
	  dcopy(l,Vbc->bface[0],1,X.z+cnt,1);
	cnt += l;
	
	dzero(cnt,X.x,1);
	dzero(cnt,X.y,1);
	drot(cnt,X.z,1,X.x,1,ct,-st);
	drot(cnt,X.x,1,X.y,1,cp,-sp);
	
	dcopy(cnt,X.x,1,Winfo.U_im[nbc][i],1);
	dcopy(cnt,X.y,1,Winfo.V_im[nbc][i],1);
	dcopy(cnt,X.z,1,Winfo.W_im[nbc][i],1);
      }
      
      nbc++;
      break;
    default:
      break;
    }
  }
  
  free(X.x); free(X.y); free(X.z); free(r); 
}

static double **wom_r, **wom_i;

static void InitWomField(Domain *omega){

  if(!wom_r){
    double R, *axs, nrm[3], *r;
    double phi, cp, sp, theta, ct, st, *x, *y, *z,wn;
    double wnum, *ccos, *csin, *f, ccos0, *ReSol, *ImSol;
    int i, j,l,qt,nbc,eid,face,cnt,nmodes;
    Coord   X;
    Element *E;
    Bndry  *Ubc,*Vbc;
    
    R      = omega->Tfun->info.wexpand.radius;
    axs    = omega->Tfun->info.wexpand.axispt;
    wnum   = omega->Tfun->info.wexpand.wnum;
    ccos   = omega->Tfun->info.wexpand.ccos;
    csin   = omega->Tfun->info.wexpand.csin;
    f      = omega->Tfun->info.wexpand.raw;
    nmodes = omega->Tfun->info.wexpand.nmodes;

    // need to leave original untouched 
    dcopy(3,omega->Tfun->info.wexpand.axisnm,1,nrm,1);
    
    //Note init should be called before womfield 
    
    X.x   = dvector(0,QGmax*QGmax*QGmax-1);
    X.y   = dvector(0,QGmax*QGmax*QGmax-1);
    X.z   = dvector(0,QGmax*QGmax*QGmax-1);
    ReSol = dvector(0,QGmax*QGmax*QGmax-1);
    ImSol = dvector(0,QGmax*QGmax*QGmax-1);
    r     = dvector(0,QGmax*QGmax*QGmax-1);
  
    phi  = atan2(nrm[1], nrm[0]);
    cp   = cos(phi); sp = sin(phi);
    drot(1,&nrm[0],1,&nrm[1],1,cp,sp);
    theta = atan2(nrm[0],nrm[2]);
    ct    = cos(theta); st = sin(theta);
    
    wom_r = dmatrix(0,nmodes-1,0,omega->U->htot-1);
    wom_i = dmatrix(0,nmodes-1,0,omega->U->htot-1);

    dzero(nmodes*omega->U->htot,wom_r[0],1);
    dzero(nmodes*omega->U->htot,wom_i[0],1);

    /* find coordinates of all quadrature points on region */
    cnt = 0;
    for(E = omega->U->fhead; E; E = E->next,cnt+=qt){
      qt = E->qa*E->qb*E->qc;
      
      E->coord(&X);
      
      /* subtract offset */
      dsadd(qt,-axs[0],X.x,1,X.x,1);
      dsadd(qt,-axs[1],X.y,1,X.y,1);
      dsadd(qt,-axs[2],X.z,1,X.z,1);
      
      /* rotate to align normal with z axis */
      drot(qt,X.x,1,X.y,1,cp,sp);
      drot(qt,X.z,1,X.x,1,ct,st);
      
      /* calculate scaled radius */
      dvmul (qt,X.x,1,X.x,1,r,1);
      dvvtvp(qt,X.y,1,X.y,1,r,1,r,1);
      dvsqrt(qt,r,1,r,1);
      dscal (qt,1./R,r,1);
      
      for (i = 0; i < omega->Tfun->info.wexpand.nmodes; ++i){
	/* find velocity in terms of real and imaginary parts 
	   for a single fourier mode for every quadrature point */
	velocity(i, qt, ccos[i], csin[i], r, R, i*wnum, 
		 wom_r[i]+cnt, wom_i[i]+cnt);
	
	dscal(qt,omega->Tfun->info.wexpand.scal0,wom_r[i]+cnt,1);
	dscal(qt,omega->Tfun->info.wexpand.scal0,wom_i[i]+cnt,1);
      } 
    }
    
    free(ReSol); free(ImSol); free(X.x); free(X.y); free(X.z); free(r); 
    free(ccos); free(csin);
  }
}

void SetWomField(Domain *omega, double *u, double *v, double *w, double time){
  register int i,qt;
  int nmodes,cnt;
  double nrm[3], wnum, phi, cp, sp, theta, ct, st, *ccos, *csin;
  double *umod,*vmod,*wmod;
  Element *E;

  if(!wom_r)
    InitWomField(omega);

  umod = dvector(0,omega->U->htot-1);
  vmod = dvector(0,omega->U->htot-1);
  wmod = dvector(0,omega->U->htot-1);

  wnum   = omega->Tfun->info.wexpand.wnum;
  ccos   = omega->Tfun->info.wexpand.ccos;
  csin   = omega->Tfun->info.wexpand.csin;
  nmodes = omega->Tfun->info.wexpand.nmodes;
  
  // need to leave original untouched 
  dcopy(3,omega->Tfun->info.wexpand.axisnm,1,nrm,1);
  
  phi   = atan2(nrm[1], nrm[0]);
  cp    = cos(phi);   sp = sin(phi);
  drot(1,&nrm[0],1,&nrm[1],1,cp,sp);

  theta = atan2(nrm[0],nrm[2]);
  ct    = cos(theta); st = sin(theta);

  dzero(omega->U->htot,u,1);
  dzero(omega->U->htot,v,1);
  dzero(omega->U->htot,w,1);

  cnt = 0;
  for(E = omega->U->fhead; E; E = E->next,cnt+=qt){
    qt      = E->qa*E->qb*E->qc;    
    for (i = 0; i < nmodes; ++i){
	
      // Real component
      dzero(qt,umod,1);
      dzero(qt,vmod,1);
      dcopy(qt,wom_r[i]+cnt,1,wmod,1);
      drot(qt,wmod,1,umod,1,ct,-st);
      drot(qt,umod,1,vmod,1,cp,-sp);	
	
      daxpy(qt,cos(i*wnum*time),umod,1,u+cnt,1);
      daxpy(qt,cos(i*wnum*time),vmod,1,v+cnt,1);
      daxpy(qt,cos(i*wnum*time),wmod,1,w+cnt,1);
	
      // Imaginary component
      dzero(qt,umod,1);
      dzero(qt,vmod,1);
      dcopy(qt,wom_i[i]+cnt,1,wmod,1);
      drot(qt,wmod,1,umod,1,ct,-st);
      drot(qt,umod,1,vmod,1,cp,-sp);
	
      daxpy(qt,-sin(i*wnum*time),umod,1,u+cnt,1);
      daxpy(qt,-sin(i*wnum*time),vmod,1,v+cnt,1);
      daxpy(qt,-sin(i*wnum*time),wmod,1,w+cnt,1);
    } 
  }
  
  free(umod); 
  free(vmod); 
  free(wmod);
}

void WomError(Domain *omega, double time){
  int htot = omega->U->htot;

  double *u,*v,*w;
  
  u = dvector(0,htot-1);
  v = dvector(0,htot-1);
  w = dvector(0,htot-1);
  
  if(!wom_r)
    InitWomField(omega);
  
  SetWomField(omega,u,v,w,time);
  
  if(omega->U->fhead->state == 't'){
    omega->U->Trans(omega->U,J_to_Q);
    omega->V->Trans(omega->V,J_to_Q);
    omega->W->Trans(omega->W,J_to_Q);
  }

#if 0

  dvsub(htot,omega->U->base_h,1,u,1,omega->U->base_h,1);
  dvsub(htot,omega->V->base_h,1,v,1,omega->V->base_h,1);
  dvsub(htot,omega->W->base_h,1,w,1,omega->W->base_h,1);

  fprintf(stdout,"Linf Error (u,v,w): %lg %lg %lg \n",
	  fabs(omega->U->base_h[idamax(htot,omega->U->base_h,1)]),
	  fabs(omega->V->base_h[idamax(htot,omega->V->base_h,1)]),
	  fabs(omega->W->base_h[idamax(htot,omega->W->base_h,1)]));

  //dcopy(htot,u,1,omega->U->base_h,1);
  //dcopy(htot,v,1,omega->V->base_h,1);
  //dcopy(htot,w,1,omega->W->base_h,1);

  omega->U->Trans(omega->U,Q_to_J);
  omega->V->Trans(omega->V,Q_to_J);
  omega->W->Trans(omega->W,Q_to_J);

#else
  dvsub(htot,omega->U->base_h,1,u,1,u,1);
  dvsub(htot,omega->V->base_h,1,v,1,v,1);
  dvsub(htot,omega->W->base_h,1,w,1,w,1);

  fprintf(stdout,"Linf Error (u,v,w): %lg %lg %lg \n",
	  fabs(u[idamax(htot,u,1)]),
	  fabs(v[idamax(htot,v,1)]),fabs(w[idamax(htot,w,1)]));

  omega->U->Set_state('t');
  omega->V->Set_state('t');
  omega->W->Set_state('t');

#endif

  free(u); free(v); free(w);
}

void SetWomSol(Domain *omega, double time, int save){
   int i,j,k,m,n,nmod,l,eid,face,cnt,dim,nfv;
   double wn,w;
   Bndry *Ubc[MAXDIM];
   Element **U;
   static double Store_Time = -99999.99;
   static int Je = iparam("INTYPE");

   if(Store_Time == time) return;  // fail safe against multiple calls
   Store_Time = time;

   nmod = omega->Tfun->info.wexpand.nmodes;
   wn   = omega->Tfun->info.wexpand.wnum;
   U    = omega->U->flist;
   dim  = U[0]->dim();
   
   /* for each boundary face */
   for(i = 0; i < Winfo.nbcs; ++i){
     Ubc[0] = omega->Ubc + Winfo.bcid[i];
     Ubc[1] = omega->Vbc + Winfo.bcid[i];
     Ubc[2] = omega->Wbc + Winfo.bcid[i];
     eid    = Ubc[0]->elmt->id;
     face   = Ubc[0]->face;
     nfv    = U[eid]->Nfverts(face);

     if(save){
       nfv = U[eid]->Nfverts(face);
       for(j = Je; j > 0; --j){
	 for(n = 0; n < nfv; ++n)
	   for(k = 0; k < dim; ++k)
	     Ubc[k]->bvert[j*nfv+n] = Ubc[k]->bvert[(j-1)*nfv+n];
	 
	 for(n = 0; n < nfv; ++n){
	   l = Ubc[0]->elmt->edge[U[eid]->ednum(face,n)].l;
	   for(k = 0; k < dim; ++k)
	     dcopy(l,Ubc[k]->bedge[n]+(j-1)*l,1,Ubc[k]->bedge[n]+j*l,1);
	 }
	 
	 l  =  U[eid]->face[face].l;
	 if(nfv == 3)
	   l = l*(l+1)/2;
	 else
	   l = l*l;
	 for(k = 0; k < dim; ++k)
	   dcopy(l,Ubc[k]->bface[0]+(j-1)*l,1,Ubc[k]->bface[0]+j*l,1);
       }
     }
     
     /* copy velocity soln into the next time step level 
	the first time round, garbage is copied */
     for (m = 0; m < dim; ++m)
       dcopy(nfv,Ubc[m]->bvert,1,Ubc[m]->bvert+nfv,1);

     /* sum wom to get the velocity on VERTICES for time t */
     for (j = 0; j < nfv; ++j) {
      /* fill in the first time step level starting with zeroth fourier mode */
       Ubc[0]->bvert[j] = Winfo.U_re[i][0][j];
       Ubc[1]->bvert[j] = Winfo.V_re[i][0][j];
       Ubc[2]->bvert[j] = Winfo.W_re[i][0][j];
       
       /* sum high fourier modes */
       for (k = 1; k < nmod; ++k){
	 Ubc[0]->bvert[j] += Winfo.U_re[i][k][j]*cos(k*wn*time) - 
	   Winfo.U_im[i][k][j]*sin(k*wn*time);
	 Ubc[1]->bvert[j] += Winfo.V_re[i][k][j]*cos(k*wn*time) - 
	   Winfo.V_im[i][k][j]*sin(k*wn*time);
	 Ubc[2]->bvert[j] += Winfo.W_re[i][k][j]*cos(k*wn*time) - 
	   Winfo.W_im[i][k][j]*sin(k*wn*time);
       }
     }
     
     cnt = nfv;
     /* sum wom to get the velocity on dim number of EDGEs for time t */
     for (j = 0; j < nfv; ++j) {
       l = U[eid]->edge[U[eid]->ednum(face,j)].l;
       /* copy velocity soln into the next time step level */
       for (m = 0; m < dim; ++m)
	 dcopy(l,Ubc[m]->bedge[j],1,Ubc[m]->bedge[j]+l,1);

       for (k = 0; k < l; ++k) {
	 /* fill in the first time step level with zeroth fourier mode */
	 Ubc[0]->bedge[j][k] = Winfo.U_re[i][0][cnt+k];
	 Ubc[1]->bedge[j][k] = Winfo.V_re[i][0][cnt+k];
	 Ubc[2]->bedge[j][k] = Winfo.W_re[i][0][cnt+k];
	 for (m = 1; m<nmod; ++m) {
	   Ubc[0]->bedge[j][k] += Winfo.U_re[i][m][cnt+k]*cos(m*wn*time) - 
	     Winfo.U_im[i][m][cnt+k]*sin(m*wn*time);
	   Ubc[1]->bedge[j][k] += Winfo.V_re[i][m][cnt+k]*cos(m*wn*time) - 
	     Winfo.V_im[i][m][cnt+k]*sin(m*wn*time);
	   Ubc[2]->bedge[j][k] += Winfo.W_re[i][m][cnt+k]*cos(m*wn*time) - 
	     Winfo.W_im[i][m][cnt+k]*sin(m*wn*time);
	 }
       }
       cnt += l;
     }  

     /* sum wom to get the velocity on FACES for time t */
     l = U[eid]->face[face].l;
     if(nfv == 3)
       l = l*(l+1)/2;
     else
       l = l*l;
     
     if(l){
       /* copy velocity soln into the next time step level */
       for(m = 0; m < dim; ++m)
	 dcopy(l,Ubc[m]->bface[0],1,&Ubc[m]->bface[0][l],1);
       for (k = 0; k < l; ++k) {
	 Ubc[0]->bface[0][k] = Winfo.U_re[i][0][cnt+k];
	 Ubc[1]->bface[0][k] = Winfo.V_re[i][0][cnt+k];
	 Ubc[2]->bface[0][k] = Winfo.W_re[i][0][cnt+k];
	 for (j = 1; j < nmod; ++j){
	   Ubc[0]->bface[0][k] += Winfo.U_re[i][j][cnt+k]*cos(j*wn*time) - 
	     Winfo.U_im[i][j][cnt+k]*sin(j*wn*time);
	   Ubc[1]->bface[0][k] += Winfo.V_re[i][j][cnt+k]*cos(j*wn*time) - 
	     Winfo.V_im[i][j][cnt+k]*sin(j*wn*time);
	   Ubc[2]->bface[0][k] += Winfo.W_re[i][j][cnt+k]*cos(j*wn*time) - 
	     Winfo.W_im[i][j][cnt+k]*sin(j*wn*time);
	 }
       }
     }
   }
}

void SaveWomSol(Domain *omega){
   int i,j,k,m,n,l,eid,face,cnt,dim,nfv;
   Bndry *Ubc[MAXDIM];
   Element **U;
   static int Je = iparam("INTYPE");

   U   = omega->U->flist;
   dim = U[0]->dim();
   
   /* for each boundary face */
   for(n = 0; n < Winfo.nbcs; ++n){
     Ubc[0] = omega->Ubc + Winfo.bcid[n];
     Ubc[1] = omega->Vbc + Winfo.bcid[n];
     Ubc[2] = omega->Wbc + Winfo.bcid[n];
     eid    = Ubc[0]->elmt->id;
     face   = Ubc[0]->face;
     nfv    = U[eid]->Nfverts(face);

     for(j = Je; j > 0; --j){
       for(i = 0; i < nfv; ++i)
	 for(k = 0; k < dim; ++k)
	   Ubc[k]->bvert[j*nfv+i] = Ubc[k]->bvert[(j-1)*nfv+i];
       
       for(i = 0; i < nfv; ++i){
	 l = U[eid]->edge[U[eid]->ednum(face,i)].l;
	 for(k = 0; k < dim; ++k)
	   dcopy(l,Ubc[k]->bedge[i]+(j-1)*l,1,Ubc[k]->bedge[i]+j*l,1);
       }
       
       l  =  U[eid]->face[face].l;
       if(nfv == 3)
	 l = l*(l+1)/2;
       else
	 l = l*l;
       for(k = 0; k < dim; ++k)
	 dcopy(l,Ubc[k]->bface[0]+(j-1)*l,1,Ubc[k]->bface[0]+j*l,1);
     }
   }
}

static void velocity(int k, int qt, double ccos, double csin, double *r,
		     double R, double wnum, double *re, double *im){
  int i;
  double alpha,mu=dparam("KINVIS"),rho=1.0;
  Cmplx jc1,jc2,a,b,c;
  int  KODE=2,NZ, IERR,N=1;
  double FNU=0.;
  
  if (k) {
    alpha = R*sqrt(wnum/mu);
    jc2.re = -alpha/SQRT2;
    jc2.im = alpha/SQRT2;
    for (i=0;i<qt;++i) {
      
      jc1.re = jc2.re*r[i];
      jc1.im = jc2.im*r[i];
      
      zbesj_(&jc1.re, &jc1.im, &FNU, &KODE, &N, &c.re, &c.im, &NZ, &IERR);
//      if (IERR) {printf("A %d error\n",IERR);exit(-1);}
      zbesj_(&jc2.re, &jc2.im, &FNU, &KODE, &N, &b.re, &b.im, &NZ, &IERR);
//      if (IERR) {printf("B %d error\n",IERR);exit(-1);}
      a.re = b.re - c.re;
      a.im = b.im - c.im;
      c.re = csin/(rho * wnum);
      c.im = ccos/(rho * wnum);
      
      re[i] = (a.re*b.re*c.re + a.im*b.im*c.re -
	       a.im*b.re*c.im + a.re*b.im*c.im)/(b.re*b.re + b.im*b.im);
      im[i] = (a.re*b.re*c.im + a.im*b.im*c.im + 
	       a.im*b.re*c.re - a.re*b.im*c.re)/(b.re*b.re + b.im*b.im);
    }
  } 
  else {
    for (i=0;i<qt;++i) {
      re[i] = ccos*R*R*(r[i]*r[i]-1.)/(4.*mu);
      im[i] = 0.;
    }
  }
}

/* calculates the value of the real and imaginary parts of the
   mass flow equation */
static void massflow(int k, double *ccos, double *csin, double R, double wnum){
  double alpha,mu=dparam("KINVIS"),rho=1.0;
  Cmplx jc,a,b,c;
  int  KODE=2,NZ, IERR,N=2;
  double FNU=0., fac, R2 = R*R, cyr[2], cyi[2];
  
  if (k) {
    alpha = R*sqrt(wnum/mu);
    jc.re = -alpha/SQRT2;
    jc.im = alpha/SQRT2;
    zbesj_(&jc.re, &jc.im, &FNU, &KODE, &N, cyr, cyi, &NZ, &IERR);
//    if (IERR) {printf("massflow: %d error\n",IERR);exit(-1);}
    b.re = cyr[0]; b.im = cyi[0];
    a.re = cyr[1]; a.im = cyi[1];
    fac  = SQRT2 / alpha;
    c.re = fac;
    c.im = fac;
    
    jc.re = (a.re*b.re*c.re + a.im*b.im*c.re - 
	     a.im*b.re*c.im + a.re*b.im*c.im) / (b.re*b.re + b.im*b.im);
    jc.im = (a.re*b.re*c.im + a.im*b.im*c.im + 
	     a.im*b.re*c.re - a.re*b.im*c.re) / (b.re*b.re + b.im*b.im);
    
    jc.re = jc.re + 1.;
    fac = (rho*wnum)/((jc.re*jc.re+jc.im*jc.im)*M_PI*R2);
      
    a.re = -fac * (jc.re * (*csin) + jc.im * (*ccos) );
    a.im =  fac * (jc.re * (*ccos) - jc.im * (*csin) );
    *ccos = a.re;
    *csin = a.im;
    
  } else {
    *ccos = - 8. * mu * (*ccos) / (M_PI * R2 * R2);
    *csin = 0.;
  }
}

static double integ(int m, double a, double b, double *f);

/* calculates the Fourier coefficients for the pressure gradient signal */
static void encode(int m, double *f, int n, double *ccos, 
		   double *csin, double T){
  int j, k;
  double *tmp;
  double dt = T/(m-1);
  
  tmp = (double *)calloc(m,sizeof(double));
  ccos[0] = integ(m, 0., T, f)/T;
  csin[0] = 0.;
  
  for (j=1;j<n;++j) {
    for (k=0;k<m;++k)
      tmp[k] = f[k]*cos(2.*M_PI*j*k*dt/T);
    ccos[j] = integ(m, 0., T, tmp)*2./T;
  }
  
  for (j=1;j<n;++j) {
    for (k=0;k<m;++k)
      tmp[k] = f[k]*sin(2.*M_PI*j*k*dt/T);
    csin[j] = integ(m, 0., T, tmp)*2./T;
  }
  
  free(tmp);
}

static double decode(int n, double *ccos, double *csin, double t, double T){
  int k;
  double w;
  
  w = ccos[0];
  for(k=1;k<n;++k)
    w += ccos[k]*cos(2*M_PI*k*t/T) + csin[k]*sin(2*M_PI*k*t/T);
  
  return w;
}

static double simpson(int m, double a, double b, double *f);
static double trapezium(int m, double a, double b, double *f);

static double integ(int m, double a, double b, double *f){
  double ans;
  
  ans = trapezium(m, a, b, f);
  
  return ans;
}

static double simpson(int m, double a, double b, double *f){
  int i;
  double u=0.;
  
  u = f[0] + f[m-1];
  for (i=1;i<m-1;i+=2)
    u += 4.*f[i];
  for (i=2;i<m-2;i+=2)
    u += 2.*f[i];
  
  return (b-a)*u/((m-1)*3.);
}

static double trapezium(int m, double a, double b, double *f){
  int i;
  double u=0.;
  
  for (i=0;i<m-1;++i)
    u += f[i]+f[i+1];
  
  return 0.5*(b-a)*u/(m-1);
}
#endif
