/*
 * Velocity interpolation routines for inflow
 * 
 * ------------------------------------------------------------------------- */
#include <stdio.h>
#include <ctype.h>
#include <stdarg.h>
#include <math.h>
#include "nektar.h"

/* 
 * Redefine vec_init() and set_vector() 
 * ---------------------------------------------------------------------- */
 
#ifdef VELINTERP
static char     dir;

typedef  struct {
  int    nr;
  int    nt;
  double offset;
  double r;
  double centre[3];
  double normal[3];
  double **array;
} VInt;
VInt Vdata;

static void VInterpInit(char *file);


#define TANTOL 1e-10
/* new atan2 function to stop Nan on atan(0,0)*/
static double atan2_proof (double x, double y)
{
  if (fabs(x) + fabs(y) > TANTOL) return (atan2(x,y));
  else return (0.);
}
#define atan2 atan2_proof

void vector_def(char *s, char *f)
{
  static int init=1;

  if(init){
    VInterpInit("inflow");
    init = 0;
  }
  while (isspace(*f++));
  dir = *--f;
  return;
}

static void VInterp(int n, double *x, double *y, double *z, 
		    double *f, char dir);
void vector_set(int n, ...)
{
  double  *x, *y, *z, *f, t;

  va_list ap;

  va_start(ap, n);
  x = va_arg(ap, double *);
  y = va_arg(ap, double *);
  z = va_arg(ap, double *);
  f = va_arg(ap, double *);
  va_end  (ap);

  VInterp(n, x, y, z, f, dir);
  return;
}

/* setup Vdata from data in file inflow */
static void VInterpInit(char *file){
  register int i,j;
  FILE     *fp;
  double   fac;
  char     buf[BUFSIZ];
  
  fp   = fopen(file,"r");
  fgets  (buf,BUFSIZ,fp);
  sscanf (buf,"%lf%lf%lf",Vdata.centre,Vdata.centre+1,Vdata.centre+2);

  fgets  (buf,BUFSIZ,fp);
  sscanf (buf,"%lf%lf%lf",Vdata.normal,Vdata.normal+1,Vdata.normal+2);
  fac = sqrt(Vdata.normal[0]*Vdata.normal[0] + Vdata.normal[1]*Vdata.normal[1]
	     + Vdata.normal[2]*Vdata.normal[2]);

  Vdata.normal[0] /= fac;
  Vdata.normal[1] /= fac;
  Vdata.normal[2] /= fac;
  
  fgets  (buf,BUFSIZ,fp);
  sscanf (buf,"%lf%lf\n",&Vdata.r,&Vdata.offset);

  fgets  (buf,BUFSIZ,fp);
  sscanf (buf,"%d%d",&Vdata.nr,&Vdata.nt);
  
  Vdata.array = dmatrix(0,Vdata.nr-1,0,Vdata.nt-1);
  
  for(j = 0; j < Vdata.nt; ++j)
    for(i = 0; i < Vdata.nr; ++i)
      fscanf(fp,"%lf",Vdata.array[i]+j);
  /*Vdata.array[i][j] = i/(double)(Vdata.nr-1)*(1.0-i/(double)(Vdata.nr-1))*
	(sin(2*M_PI*j/(double)(Vdata.nt-1))); */
  
}

/* Note this routine overwrites x,y and z */
static void VInterp(int n, double *X, double *Y, double *Z, 
		    double *f, char dir){
  register int i,j;
  int    nt,nr;
  double theta, r, xn, yn,dr,dt;
  double phi,cp,sp,ct,st;
  double *x,*y,*z;

  /* copy x,y,z so they are not overwritten */
  x = dvector(0,n-1);
  y = dvector(0,n-1);
  z = dvector(0,n-1);
  dcopy(n,X,1,x,1);
  dcopy(n,Y,1,y,1);
  dcopy(n,Z,1,z,1);
  
  /* translate co-ordinates so that given point is origin */
  dsadd(n, -Vdata.centre[0], x, 1, x, 1);
  dsadd(n, -Vdata.centre[1], y, 1, y, 1);
  dsadd(n, -Vdata.centre[2], z, 1, z, 1);
    
  /* rotate so that normal is  aligned with z axis */
  phi = atan2(Vdata.normal[1],Vdata.normal[0]);
  cp  = cos(phi); sp = sin(phi);
  xn = Vdata.normal[0]; yn = Vdata.normal[1];
  drot(1,&xn,1,&yn,1,cp,sp);
  theta = atan2(xn,Vdata.normal[2]);
  ct = cos(theta); st = sin(theta);

  drot(n,x,1,y,1,cp,sp);
  drot(n,z,1,x,1,ct,st);

  /* determine axial velocity in local coordinate */
  nr = 0; nt = 0;
  for(i = 0; i < n; ++i){
    r     = sqrt (x[i]*x[i] + y[i]*y[i]);
    theta = atan2(y[i],x[i]) + Vdata.offset;
    if(theta < 0.0)    theta += 2*M_PI;
    if(theta > 2*M_PI) theta -= 2*M_PI;

    /* assume theta is equisapced between 0 and PI and r is 
       equispaced between 0 and Vdata.r; */
    dr     = Vdata.r/(double)(Vdata.nr-1);
    dt     = 2*M_PI /(double)(Vdata.nt-1);
    nr     = min((int) (r/dr),Vdata.nr-2);
    nt     = min((int) (theta/dt),Vdata.nt-2);
    
    dr     = (r - nr*dr)/dr;
    dt     = (theta - nt*dt)/dt;

    /* put local axial component in the z direction */
    z[i] = Vdata.array[nr][nt]*(1.0-dr)*(1.0-dt) + 
      Vdata.array[nr][nt+1]*(1.0-dr)*dt + 
	Vdata.array[nr+1][nt]*dr*(1.0-dt) + 
	  Vdata.array[nr+1][nt+1]*dr*dt;
    /* zero out x,  y */
    x[i] = y[i] = 0.0;
  }
  /* rotate solution back solution */
  drot(n,z,1,x,1,ct,-st);
  drot(n,x,1,y,1,cp,-sp);

  switch(dir){
  case 'u':
    dcopy(n,x,1,f,1);
    break;
  case 'v':
    dcopy(n,y,1,f,1);
    break;
  case 'w':
    dcopy(n,z,1,f,1);
    break;
  }
  
  free(x); free(y); free(z);
}
#endif
