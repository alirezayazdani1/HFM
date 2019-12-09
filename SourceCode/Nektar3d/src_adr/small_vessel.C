/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $Source: 
 * $Revision:
 * $Date: 
 * $Author:
 * $State: 
 *---------------------------------------------------------------------------*/

#include "nektar.h"
#include "pbc_1d.h"
#include <math.h>
#include <string.h>
#include "veclib.h"

#ifdef PBC_1D

void ComputeFj(double Ws_n, double *Re_part, double *Im_part);
my_dcmplx ComputeFj(double Ws_n);
my_dcmplx ComputeWaveSpeed(double Ws_n, double radius, double Compliance);
void M2Ptransform(int Nmodes, double om, my_dcmplx *MODES, int N, double Time_period, double *f);
void M2Ptransform(int Nmodes, double om, my_dcmplx* MODES , int N, double Time_period, int index_start, int index_end, double* f);
void Convolve(int N, double* f1, double* f2, double* ans);
void Convolve(int N, double* f1, double* f2, double* ans, int time_index);
void ParallelConvolve(int N, double* f1, double* f2, double* ans, int it);
void ParallelConvolve(int N, double* f1, double* f2, double* ans, int index_start, int index_end, int it);
void parallelConvolveTest(int N, double* f1, double* f2, double* ans, int index_start, int index_end, int it);
void FilterConv(int N, int Nmodes, double  *f);
void FilterDFT(int N, int Nmodes, double  *f);
void FilterDFT(int N, int Nmodes, double  *f, double *Acos_sin, int FLAG_averag);
void FilterSmooth(int N, double *f);
void Linspace(int N, double *x);
void GetBesselj0Local(double r, double *real_part, double *imag_part);
void GetBesselj1Local(double r, double *real_part, double *imag_part);


SMALL_VESSEL::SMALL_VESSEL(double ro)
  : radius(ro), Length(Lrr_small_atree*ro),
    Zo(0), Zl(0),
    type('\0'),
    Compliance(1.5*(M_PI*radius*radius)/(k1_small_atree*exp(k2_small_atree*radius)+k3_small_atree)),
    Ws(sqrt(radius*radius*omega_small_atree/ni_small_atree)),
    child_1(0), child_2(0)
{
  if(ro < Rmin_small_atree)
    type = 'T';
  else
  {
    type = 'C';
    child_1 = new SMALL_VESSEL(radius*0.9);
    child_2 = new SMALL_VESSEL(radius*0.6);
  }
}


void my_dcmplx::set_val(double re, double im){
    real = re;
    imag = im;
}

void my_dcmplx::my_dcmplx_axpy(double a, my_dcmplx x, my_dcmplx y){
    real = a*x.real+y.real;
    imag = a*x.imag+y.imag;
}

void my_dcmplx::my_dcmplx_axpy(double a, my_dcmplx x, double b, double y){
    real = a*x.real+b*y;
    imag = a*x.imag;
}
void my_dcmplx::my_dcmplx_dot(my_dcmplx x, my_dcmplx y){
    
    real = x.real*y.real-x.imag*y.imag;
    imag = x.real*y.imag+x.imag*y.real;
}

void my_dcmplx::my_dcmplx_inv_dot(my_dcmplx x, my_dcmplx y){
    
  /* ans = x/y */

    double den;
    den = y.real*y.real+y.imag*y.imag;

    real = ( x.real*y.real+x.imag*y.imag)/den;
    imag = (-x.real*y.imag+x.imag*y.real)/den;
}

void my_dcmplx::my_dcmplx_sin(my_dcmplx x){
/* sin(a+ib) = sin(a)*cosh(b)+i*cos(a)*sinh(b) */
    real = sin(x.real)*cosh(x.imag);
    imag = cos(x.real)*sinh(x.imag);

}

void my_dcmplx::my_dcmplx_cos(my_dcmplx x){
/* cos(a+ib) = cos(a)*cosh(b)-i*sin(a)*sinh(b) */
    real =  cos(x.real)*cosh(x.imag);
    imag = -sin(x.real)*sinh(x.imag);
}

void my_dcmplx::my_dcmplx_sqrt(my_dcmplx x){
/*    sqrt(x) = sqrt(r)*exp(i*teta/2)
      r = x.real*x.real+x.imag*x.imag
      teta = atan(x.imag/x.real)
*/
    double r,teta;
    r = sqrt(x.real*x.real+x.imag*x.imag);
    if (r < (1.0e-15)) {
	real = 0.0;
	imag = 0.0;
	return;
    }

    if (fabs(x.real) < (1.0e-15))
	teta = M_PI/2.0*x.imag/fabs(x.imag);
    else
	teta = atan(x.imag/x.real);

    r=sqrt(r);

    real = r*cos(teta*0.5);
    imag = r*sin(teta*0.5);
}

void my_dcmplx::my_dcmplx_conj(my_dcmplx x){
    real = x.real;
    imag = -x.imag;
}

my_dcmplx SMALL_VESSEL::GetZL(int mode){

/* compute impedance at the exit of the vessel */

  my_dcmplx ans,num,den;
  ans.set_val(0.0,0.0);//terminal resistance

  if (type == 'T')
    return ans;
  else{
    my_dcmplx Z0_child1 = child_1->GetZ0(mode);
    my_dcmplx Z0_child2 = child_2->GetZ0(mode);

    num.my_dcmplx_dot(Z0_child1,Z0_child2);
    den.my_dcmplx_axpy(1.0, Z0_child1, Z0_child2);
    ans.my_dcmplx_inv_dot(num, den);

    return ans;
  }
}


my_dcmplx SMALL_VESSEL::GetZ0(int mode){


/* compute impedance at the entrance to the vessel */
  my_dcmplx ZL,ans;

  if (mode == 0){
      ZL = GetZL(mode);

    ans.real=8.0*(ni_small_atree*density_small_atree)*
           Lrr_small_atree/M_PI/(radius*radius*radius)+ZL.real;
    ans.imag = ZL.imag;

    return ans;
  }

  my_dcmplx a,b,II,arg; 
  
  my_dcmplx wave_speed, g;

  my_dcmplx temp1,temp2;

  II.set_val(0.0,1.0);
  wave_speed = ComputeWaveSpeed(Ws*sqrt(1.0*mode), radius, Compliance);
  g.set_val(wave_speed.real*Compliance,wave_speed.imag*Compliance);

  ans.set_val(omega_small_atree*mode*Length,0.0);
  arg.my_dcmplx_inv_dot(ans,wave_speed);
  ZL = GetZL(mode);

  a.my_dcmplx_dot(II,g);
  b.my_dcmplx_dot(a,ZL);  //b = II*g*ZL;

  temp1.my_dcmplx_sin(arg);
  temp2.my_dcmplx_cos(arg);

  ans.my_dcmplx_dot(b,temp1);        //ans = II*g*ZL*sin(arg);
  a.my_dcmplx_axpy(1.0,temp2,ans);  //a = cos(arg)+II*g*ZL*sin(arg);
    
  ans.my_dcmplx_dot(b,temp2);        //ans = II*g*ZL*cos(arg);
  b.my_dcmplx_axpy(-1.0,ans,temp1); //b = sin(arg)-II*g*ZL*cos(arg);

  temp1.my_dcmplx_inv_dot(II,g);
  temp2.my_dcmplx_inv_dot(b,a);

  ans.my_dcmplx_dot(temp1,temp2);


  return ans; 
}

my_dcmplx ComputeFj(double Ws_n){

  double Re1,Im1,Re2,Im2;
  
  my_dcmplx II,num,denum;
  my_dcmplx temp1,temp2;

  II.set_val(0.0,1.0);
 
  GetBesselj1Local(Ws_n, &Re2, &Im2);
  GetBesselj0Local(Ws_n, &Re1, &Im1);

  num.set_val(Re2,Im2);
  denum.set_val(Re1,Im1);

  temp1.my_dcmplx_inv_dot(num,denum);
  temp2.my_dcmplx_axpy(-1.0,II,1.0,1.0);

  temp2.real = temp2.real/sqrt(2.0)*Ws_n;
  temp2.imag = temp2.imag/sqrt(2.0)*Ws_n;

  num.my_dcmplx_inv_dot(temp1,temp2);
  num.real *= 2.0;
  num.imag *= 2.0;


  return num;
}

my_dcmplx ComputeWaveSpeed(double Ws_n, double radius, double Compliance){

  my_dcmplx ans,temp;
  ans.set_val(0.0,0.0);  

  temp = ComputeFj(Ws_n);
  ans.my_dcmplx_axpy(-1.0,temp,1.0,1.0);
  ans.real =  ans.real*M_PI*radius*radius/(density_small_atree*Compliance);
  ans.imag =  ans.imag*M_PI*radius*radius/(density_small_atree*Compliance);

  temp.my_dcmplx_sqrt(ans);

  return temp;
}


void M2Ptransform(int Nmodes, double om, my_dcmplx* MODES , int N, double Time_period, double* f){
/*transform from fourier space into physical  */
  
  register int it, n,k;
  double t, Dt;
  my_dcmplx ans,temp1,temp2;
  

  Dt = Time_period/(-1.0+N);

  for (it = 0; it < N; it++){

    ans.set_val(0.0,0.0);
    k = 0;
    t = Dt*it;
    for (n = -Nmodes; n <= Nmodes; n++){
        temp1.real = cos(om*((double) n)*t);
        temp1.imag = sin(om*((double) n)*t);

	temp2.my_dcmplx_dot(MODES[k],temp1);
        ans.real += temp2.real;
        ans.imag += temp2.imag;
	k++;
    }

    f[it] = ans.real;
  }
}


void M2Ptransform(int Nmodes, double om, my_dcmplx* MODES , int N, double Time_period, int index_start, int index_end, double* f){
/*transform from fourier space into physical  */

  register int it, n,k;
  double t, Dt;
  my_dcmplx ans,temp1,temp2;


  Dt = Time_period/(-1.0+N);

  for (it = index_start; it < index_end; it++){

    ans.set_val(0.0,0.0);
    k = 0;
    t = Dt*it;
    for (n = -Nmodes; n <= Nmodes; n++){
        temp1.real = cos(om*((double) n)*t);
        temp1.imag = sin(om*((double) n)*t);

        temp2.my_dcmplx_dot(MODES[k],temp1);
        ans.real += temp2.real;
        ans.imag += temp2.imag;
        k++;
    }

    f[it-index_start] = ans.real;
  }
}

void ParallelConvolve(int N, double* f1, double* f2, double* ans, int index_start, int index_end, int it ){

  /*  perform piecewise convolution 
     array f1 is replicated at each processore but array f2 is saved piecewise 
  */


  /* assume f[0] = f[N]  */
  int itau,index;
  register double sum;
  int it2;
  sum = 0.0;

  if (index_end <= (N-it))
    it2 = index_end;
  else
    it2 = N-it;

  for (itau = index_start; itau < it2; itau++)
    sum += f1[it+itau]*f2[itau-index_start];


  if (it2 < index_start)
   it2 = index_start; 

  for (itau = it2; itau < index_end; itau++)
    sum += f1[it-N+itau]*f2[itau-index_start];

  ans[0] = sum/N;

}


void ParallelConvolve(int N, double* f1, double* f2, double* ans, int it ){


  /* assume f[0] = f[N]  */
  int itau,index;
  register double sum;

  sum = 0.0;

  for (itau = 0; itau <(N-it); itau++)
    sum += f1[it+itau]*f2[itau];

  for (itau = (N-it); itau < N; itau++)
    sum += f1[it-N+itau]*f2[itau];

  ans[0] = sum/N;
}


void Convolve(int N, double* f1, double* f2, double* ans, int time_index ){

  /* assume f[0] = f[N]  */
  int it,itau,index;
  register double sum;

  it = time_index;

  sum = 0.0;

  for (itau = 0; itau < N; itau++){

    index = it+itau;
    if (index  >= N) index = index-N;

    sum += f1[index]*f2[itau];
  }
  ans[0] = sum/N;
}

void ParallelConvolveTest(int N, double* f1, double* f2, double* ans, int index_start, int index_end, int it ){

  /* assume f[0] = f[N]  */
  int itau,index;
  register double sum;

  sum = 0.0;

  for (itau = index_start; itau < index_end; itau++){

    index = it+itau;
    if (index  >= N) index = index-N;

    sum += f1[index]*f2[itau-index_start];
  }
  ans[0] = sum/N;

  sum = 0.0;
  gdsum(ans,1,&sum);

}


void Convolve(int N, double* f1, double* f2, double* ans){

  int it,itau,index;
  double sum;

  for (it = 0; it < (N-1); it++){
    sum = 0.0;

    for (itau = 0; itau < N; itau++){

      index = it-itau;
      if (index < 0) index += N;

      sum += f1[index]*f2[itau];
    }
    ans[it] = sum/(-1.0+N);
  }
}

void FilterConv(int N, int Nmodes, double  *f){

  /* assume f is periodic function f[0] = f[N] */
  /*  */

  register int i,j,k;
  double *temp;
  double filter,dt,arg;

  dt = 1.0/N;
  temp = new double[N];
  memset(temp,'\0',N*sizeof(double));

  for (j = 0; j < N; j++){

    filter = 1.0;
    for (k = 1; k <= Nmodes; k++){
      arg = omega_small_atree*k*dt*j;
      filter += cos(arg)+sin(arg);
    }
    for (i = 0; i < N; i++){

      k = i-j;
      if (k < 0) k += N;

      temp[i] += f[k]*filter*dt;
    }
  }

  memcpy(f,temp,N*sizeof(double));

  delete[] temp;

}

void FilterDFT(int N, int Nmodes, double  *f){

  register int i,j,k;
  double *Acos,*Asin;
  double filter,dt,arg;
  Acos = new double[Nmodes+1];
  Asin = new double[Nmodes+1];  
  memset(Asin,'\0',(Nmodes+1)*sizeof(double));
  memset(Acos,'\0',(Nmodes+1)*sizeof(double));

  dt = 1.0/N;

  for (j = 0; j < N; j++){
    arg = omega_small_atree*j*dt;

    Acos[0] += f[j];
    for (k=1; k <= Nmodes; k++){
      Acos[k] += f[j]*cos(arg*k);
      Asin[k] += f[j]*sin(arg*k);
    }
  }

  Acos[0] = Acos[0]/N;
  for (k=1; k <= Nmodes; k++){
    Acos[k] = Acos[k]*2.0/N;
    Asin[k] = Asin[k]*2.0/N;
  }

  for (j = 0; j < N; j++){
    arg = omega_small_atree*j*dt;
    f[j] = Acos[0];
    for (k=1; k <= Nmodes; k++)
      f[j] += Acos[k]*cos(arg*k) + Asin[k]*sin(arg*k);
  }

  delete[] Asin;
  delete[] Acos;

}

void FilterDFT(int N, int Nmodes, double  *f, double *Acos_sin, int FLAG_average){

  register int i,j,k;
  double filter,dt,arg;
  double *Acos_sin_temp;
  double omega = 2.0*M_PI;

  if (FLAG_average == 1){
    Acos_sin_temp = new double[Nmodes*2+1];
    memcpy(Acos_sin_temp,Acos_sin,(Nmodes*2+1)*sizeof(double));
  }

  memset(Acos_sin,'\0',(Nmodes*2+1)*sizeof(double));

  dt = 1.0/N;

  for (j = 0; j < N; j++){
    arg = omega*j*dt;

    Acos_sin[0] += f[j];
    for (k=1; k <= Nmodes; k++){
      Acos_sin[k*2-1] += f[j]*cos(arg*k);
      Acos_sin[k*2] += f[j]*sin(arg*k);
    }
  }

  Acos_sin[0] = Acos_sin[0]*dt;
  for (k=1; k <= Nmodes; k++){
    Acos_sin[k*2-1] = Acos_sin[k*2-1]*2.0*dt;
    Acos_sin[k*2]   = Acos_sin[k*2]*2.0*dt;
  }

#if 1 
  /* dump modes with short wave-length */
  int n_start = Nmodes/2;
  for (k=n_start; k <= Nmodes; k++){
    arg = ((double) (k-n_start))/( (double) (Nmodes+1-n_start));
    arg = cos(0.5*M_PI*arg); 
    Acos_sin[k*2-1] *= arg;
    Acos_sin[k*2]   *= arg;
  }
#endif

  if (FLAG_average == 0){
    for (j = 0; j < N; j++){
      arg = omega*j*dt;
      f[j] = Acos_sin[0];
      for (k=1; k <= Nmodes; k++)
        f[j] += Acos_sin[k*2-1]*cos(arg*k) + Acos_sin[k*2]*sin(arg*k);
    }
  }
  else{
    double theta_N, theta_Nm1; 
    theta_N = 0.5;
    theta_Nm1 = 1.0-theta_N;
    for (j = 0; j < N; j++){
      arg = omega*j*dt;
      f[j] = (theta_N*Acos_sin[0]+theta_Nm1*Acos_sin_temp[0]);
      for (k=1; k <= Nmodes; k++){
        f[j] += (theta_N*Acos_sin[k*2-1]+theta_Nm1*Acos_sin_temp[k*2-1])*cos(arg*k) +
                (theta_N*Acos_sin[k*2]  +theta_Nm1*Acos_sin_temp[k*2])*sin(arg*k);
      }
    }
  }

  if (FLAG_average == 1)
    delete[] Acos_sin_temp;
}


void FilterSmooth(int N, double *f){

  register int i;
  double *temp;
  temp = new double[N]; 
  double a,b; 
  
  a = 1.0/6.0;
  b = 4.0/6.0;

  i = 0;
  temp[i] = a*(f[N-1]+f[i+1])+b*f[i];

  for (i = 1; i < (N-1); i++)
    temp[i] = a*(f[i-1]+f[i+1])+b*f[i];
  
  memcpy(f,temp,N*sizeof(double));

  delete[] temp;
}


void Linspace(int N, double*x){

  register int i;
  double dx,x0;
  dx = (x[N-1]-x[0])/(-1.0+N);
  x0 = x[0];
  for (i = 1; i < (N-1); i++)
    x[i] = x0+dx*i;

}


void GetBesselj0Local(double r, double *Ar, double *Br){

// computes BesselJ(0,(1-i)/sqrt(2)*r)
// Ar = real      part
// Br = imaginary part

    int k, kc, ks;
    double arg,den;

    double Rp = 0.0;
    double Ip = 0.0;

    double fk=1.0;         //fk      = k!
    double f2pk=1.0;       //f2pk    = (2+k)!
    double four_pk=1.0;    //four_pk = 4^k
    double rp2 = r*r;      //rp2     = r^2
    double rp2pk=1.0;      //rp2pk   = (r^2)^k

    for(k = 0; k <= 14; k=k+2){
        kc  = k;
        ks  = k+1;
        arg = M_PI*kc*0.5;
        den = four_pk*fk*f2pk;
        Rp += rp2pk*cos(arg)/den;

        four_pk *= 4.0;
        rp2pk   *= rp2;
        fk      *= (1.0+k);
        f2pk    *= (1.0+k);

        den = four_pk*fk*f2pk;
        arg = M_PI*ks*0.5;
        Ip  = Ip-rp2pk*sin(arg)/den;

        four_pk *= 4.0;
        rp2pk   *= rp2;
        fk      *= (1.0+ks);
        f2pk    *= (1.0+ks);
    }

    Ar[0] = Rp;
    Br[0] = -Ip;
}


void GetBesselj1Local(double r, double *Ar, double *Br){

// computes BesselJ(1,(1-i)/sqrt(2)*r)
// Ar = real      part
// Br = imaginary part

    int k, kc, ks;
    double arg,den;

    double Rp = 0.0;
    double Ip = 0.0;
    double fk=1.0;         //fk      = k!
    double f2pk=2.0;       //f2pk    = (1+k+1)! , k=0
    double four_pk=1.0;    //four_pk = 4^k
    double rp2 = r*r;      //rp2     = r^2
    double rp2pk=1.0;      //rp2pk   = (r^2)^k

    for(k = 0; k <= 14; k=k+2){
        kc  = k;
        ks  = k+1;
        arg = M_PI*kc*0.5;
        den = four_pk*fk*f2pk;
        Rp += rp2pk*cos(arg)/den;

        four_pk *= 4.0;
        rp2pk   *= rp2;
        fk      *= (1.0+k);
        f2pk    *= (2.0+k);

        den = -four_pk*fk*f2pk;
        arg = M_PI*ks*0.5;
        Ip  = Ip-rp2pk*sin(arg)/den;
        

        four_pk *= 4.0;
        rp2pk   *= rp2;
        fk      *= (1.0+ks);
        f2pk    *= (2.0+ks);
    }

    fk=Rp;  //temporary usage of variable "fk"
    Rp=Ip;
    Ip=fk;

    Ar[0] =  (Rp+Ip)*r/sqrt(2.0);
    Br[0] =  (Rp-Ip)*r/sqrt(2.0);
}

#endif
