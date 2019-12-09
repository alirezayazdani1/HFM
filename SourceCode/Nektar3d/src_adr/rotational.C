/* ------------------------------------------------------------------------- *
 * VxOmega() - Calculate the nonlinear terms in rotational form              *
 *                                                                           *
 * The following function calculates the nonlinear portion of the Navier-    *
 * Stokes equation in rotational form to minimize work.  It may be expressed *
 * as:                                                                       *
 *                                                                           *
 *     N(V) = [ v.wz - w.wy ] i + [ w.wx - u.wz ] j + [ u.wy - v.wx ] k      *
 *                                                                           *
 * where V = (u,v,w) are the three velocity components and (wx,wy,wz) are    *
 * the corresponding components of the vorticty vector.                      *
 *                                                                           *
 * RCS Information
 * ---------------
 * $Author: ssherw $
 * $Date: 2006/05/08 10:01:39 $
 * $Source: /homedir/cvs/Nektar/Nektar3d/src/rotational.C,v $
 * $Revision: 1.2 $
 * ------------------------------------------------------------------------- */

#include "nektar.h"

/* -----------------  2D Nonlinear Terms  ------------------- */

void VxOmega(Domain *omega)
{
  fprintf(stderr, "VxOmege not installed yet\n");
  exit(-1);

  return;
}



