/* ------------------------------------------------------------------------- *
 * VdgradV() - Calculate convective form of the nonlinear terms              *
 *                                                                           *
 * The following function calculates the nonlinear portion of the Navier-    *
 * Stokes equation in convective form in order to minimize aliasing errors.  *
 * In this form, the non-linear part looks like:                             *
 *                                                                           *
 *                                      dV_i                                 *
 *                         N_i = -( V_j ---- )                               *
 *                                      dx_j                                 *
 *                                                                           *
 * RCS Information                                                           *
 * ---------------                                                           
 * $Author:
 * $Date: 
 * $Source: 
 * $Revision:
 * ------------------------------------------------------------------------- */

#include "nektar.h"

static int cnt=0;

#ifdef ADR
	static int nspecs;
#endif

void VdgradV(Domain *omega)
{
  FILE   *fp;
  Element *E;

  Element_List *U    =  omega->U,  *V    =  omega->V,   *W  = omega->W,
               *Ux   =  omega->Uf, *Uy   =  omega->Vf,  *Uz = omega->Wf,
               *Vx   =  omega->Vf, *Vy   =  omega->Wf,  *Vz = omega->P,
		           *Wx   =  omega->Wf, *Wy   =  omega->P,   *Wz = omega->Pf;
#ifdef ADR
	nspecs = iparam("NSPEC");
	for (int i = 0; i < nspecs; i++){
		Element_List *T = omega->T[i], 
							   *Tx   =  omega->Tf[i], *Ty   =  omega->P,   *Tz = omega->Pf;
	
	  T->Grad(Tx,Ty,Tz,'a');
  	dvmul  (T->htot, U->base_h,1,Tx->base_h,1,Tx->base_h,1);
  	dvvtvp (T->htot, V->base_h,1,Ty->base_h,1,Tx->base_h,1,Tx->base_h,1);
  	dvvtvp (T->htot, W->base_h,1,Tz->base_h,1,Tx->base_h,1,Tx->base_h,1);
  	dvneg  (T->htot,Tx->base_h,1,Tx->base_h,1);

		set_state(Tx->fhead,'p');
	}
#endif	

  U->Grad(Ux,Uy,Uz,'a');
  dvmul  (U->htot, U->base_h,1,Ux->base_h,1,Ux->base_h,1);
  dvvtvp (U->htot, V->base_h,1,Uy->base_h,1,Ux->base_h,1,Ux->base_h,1);
  dvvtvp (U->htot, W->base_h,1,Uz->base_h,1,Ux->base_h,1,Ux->base_h,1);
  dvneg  (U->htot,Ux->base_h,1,Ux->base_h,1);
  
  V->Grad(Vx,Vy,Vz,'a');
  dvmul  (V->htot, U->base_h,1,Vx->base_h,1,Vx->base_h,1);
  dvvtvp (V->htot, V->base_h,1,Vy->base_h,1,Vx->base_h,1,Vx->base_h,1);
  dvvtvp (V->htot, W->base_h,1,Vz->base_h,1,Vx->base_h,1,Vx->base_h,1);
  dvneg  (V->htot,Vx->base_h,1,Vx->base_h,1);

  W->Grad(Wx,Wy,Wz,'a');
  dvmul  (W->htot, U->base_h,1,Wx->base_h,1,Wx->base_h,1);
  dvvtvp (W->htot, V->base_h,1,Wy->base_h,1,Wx->base_h,1,Wx->base_h,1);
  dvvtvp (W->htot, W->base_h,1,Wz->base_h,1,Wx->base_h,1,Wx->base_h,1);
  dvneg  (W->htot,Wx->base_h,1,Wx->base_h,1);

  set_state(Ux->fhead,'p'); set_state(Vx->fhead,'p'); set_state(Wx->fhead,'p');

  return;
}

void CdgradV(Domain *omega)
{
  FILE   *fp;
  Element *E;

  Element_List *U    =  omega->U,  *V    =  omega->V,   *W  = omega->W,
               *Ux   =  omega->Uf, *Uy   =  omega->Vf,  *Uz = omega->Wf,
               *Vx   =  omega->Vf, *Vy   =  omega->Wf,  *Vz = omega->P,
	       			 *Wx   =  omega->Wf, *Wy   =  omega->P,   *Wz = omega->Pf;  
  double       *cx,*cy,*cz;
  cx = dvector(0,3*U->htot-1);
  cy = cx + U->htot;
  cz = cy + U->htot;

  // store convective velocity from Uf,Vf, Wf
  dcopy  (U->htot,Ux->base_h,1,cx,1);
  dcopy  (U->htot,Uy->base_h,1,cy,1);
  dcopy  (U->htot,Uz->base_h,1,cz,1);

  U->Grad(Ux,Uy,Uz,'a');
  dvmul  (U->htot, cx,1,Ux->base_h,1,Ux->base_h,1);
  dvvtvp (U->htot, cy,1,Uy->base_h,1,Ux->base_h,1,Ux->base_h,1);
  dvvtvp (U->htot, cz,1,Uz->base_h,1,Ux->base_h,1,Ux->base_h,1);
  dvneg  (U->htot,Ux->base_h,1,Ux->base_h,1);
  
  V->Grad(Vx,Vy,Vz,'a');
  dvmul  (V->htot, cx,1,Vx->base_h,1,Vx->base_h,1);
  dvvtvp (V->htot, cy,1,Vy->base_h,1,Vx->base_h,1,Vx->base_h,1);
  dvvtvp (V->htot, cz,1,Vz->base_h,1,Vx->base_h,1,Vx->base_h,1);
  dvneg  (V->htot,Vx->base_h,1,Vx->base_h,1);

  W->Grad(Wx,Wy,Wz,'a');
  dvmul  (W->htot, cx,1,Wx->base_h,1,Wx->base_h,1);
  dvvtvp (W->htot, cy,1,Wy->base_h,1,Wx->base_h,1,Wx->base_h,1);
  dvvtvp (W->htot, cz,1,Wz->base_h,1,Wx->base_h,1,Wx->base_h,1);
  dvneg  (W->htot,Wx->base_h,1,Wx->base_h,1);

  set_state(Ux->fhead,'p'); 
  set_state(Vx->fhead,'p'); 
  set_state(Wx->fhead,'p');

  free(cx);
  return;
}
