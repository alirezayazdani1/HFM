/* ------------------------------------------------------------------------- *
 * StokesBC() - Calculate the high-order boundary conditions for Stokes flow *
 *                                                                           *
 * This routine simply sets the non-linear terms to zero.                    *
 *                                                                           *
 * RCS Information                                  
 * ---------------
 * $Author: ssherw $
 * $Date: 2006/05/08 10:01:39 $
 * $Source: /homedir/cvs/Nektar/Nektar3d/src/stokes.C,v $
 * $Revision: 1.2 $
 * ------------------------------------------------------------------------- */

#include "nektar.h"

void StokesBC(Domain *omega){
  int eDIM = omega->Uf->fhead->dim();
  /* Set the non-linear terms to zero */

  if(eDIM == 2){
    dzero(omega->Uf->htot, omega->Uf->base_h, 1);
    dzero(omega->Vf->htot, omega->Vf->base_h, 1);
    
    dzero(omega->Uf->hjtot, omega->Uf->base_hj, 1);
    dzero(omega->Vf->hjtot, omega->Vf->base_hj, 1);

    omega->Uf->Set_state('p');
    omega->Vf->Set_state('p');
  }
  else{
    dzero(omega->Uf->htot, omega->Uf->base_h, 1);
    dzero(omega->Vf->htot, omega->Vf->base_h, 1);
    dzero(omega->Wf->htot, omega->Vf->base_h, 1);
    
    dzero(omega->Uf->hjtot, omega->Uf->base_hj, 1);
    dzero(omega->Vf->hjtot, omega->Vf->base_hj, 1);
    dzero(omega->Wf->hjtot, omega->Vf->base_hj, 1);
    
    set_state(omega->Wf->fhead,'p');
  }
  return;
}

