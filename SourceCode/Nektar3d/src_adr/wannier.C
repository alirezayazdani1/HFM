/*
 * Wannier Flow 
 *
 * This is benchmark problem for Prism.  This file sets up the solution
 * for Wannier flow, an exact solution to the Stokes equations for a cylinder
 * spinning above a wall.  The boundary conditions are complicated enough
 * to be hard-coded in this file.
 * 
 * ------------------------------------------------------------------------- */
#ifdef WANNIER
#include <stdio.h>
#include <ctype.h>
#include <stdarg.h>
#include "nektar.h"

static char     dir;

#if DIM == 2
static void wannier(int n, double *x, double *y, double *u, char dir);
#else
static void wannier(int n, double *x, double *y, double *z, double *u,
		    char dir);
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
  double  *x, *y, *f;
#else
  double  *x, *y, *z, *f;
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
  wannier(n, x, y, f, dir);
#else
  wannier(n, x, y, z, f, dir);
#endif
  return;
}

/* ---------------------------------------------------------------------- */

#if DIM == 2
static void wannier(int n, double *x, double *y, double *u, char dir)
#else
static void wannier(int n, double *x, double *y, double *z, double *u,
		    char dir)
#endif
{
  double K1, K2, Y1, Y2;

  static double U     = 1.0,
                V     = 0.5, 
                R     = 0.25,
                d     = 0.5,
                s     = 0.0;
  static double A, B, C, D, F, G;
  
  register int i;

  if (s == 0.0) {          /* initialization */
    s = sqrt(d*d - R*R);
    G = (d+s) / (d-s);
    A = -(U*d / log(G) + 0.5 * d*R*V/s);
    B = 2.0*(d+s)*U/log(G) + (d+s)*R*V/s;
    C = 2.0*(d-s)*U/log(G) + (d-s)*R*V/s;
    D = -U;
    F = U / log(G);
  }

  switch (dir) {
  case 'u':
    for (i = 0; i < n; i++) {
      Y1    = y[i] + d;
      Y2    = 2.0 * Y1;
      K1    = x[i]*x[i] + (s+Y1)*(s+Y1);
      K2    = x[i]*x[i] + (s-Y1)*(s-Y1);
      
      u[i]  = -2.0/K1*(A + F*Y1) * (s+Y1 + K1/K2*(s-Y1)) - F*log(K1/K2)
	      - B/K1*( s + Y2 - Y2 * (s+Y1)*(s+Y1)/K1)
	      - C/K2*( s - Y2 + Y2 * (s-Y1)*(s-Y1)/K2) - D;
    }
    break;

  case 'v':
    for (i = 0; i < n; i++) {
      Y1    = y[i] + d;
      Y2    = 2.0 * Y1;
      K1    = x[i]*x[i] + (s+Y1)*(s+Y1);
      K2    = x[i]*x[i] + (s-Y1)*(s-Y1);

      u[i]  = 2.0*x[i]/(K1*K2) * (A + F*Y1) * (K2 - K1) -
              x[i]*B*Y2*(s+Y1)/(K1*K1) - C*x[i]*Y2*(s-Y1)/(K2*K2);
    }
    break;
#if DIM == 3
  case 'w':
    for (i = 0; i < n; i++) 
      u[i] = 0.0;
    break;
#endif
  case 'p':
    break;
  default:
    printf("unknown direction in wannier -- %c\n", dir);
    break;
  }
  
  return;
}
#endif






