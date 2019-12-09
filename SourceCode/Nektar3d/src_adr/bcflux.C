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

static void CalcTbc(Bndry *B, Element_List *V);

static double Re;
static int Je, nspecs;

void SetTBCs(Domain *omega){
  Bndry    **Tbc  =  omega->Tbc;
	Bndry			*bndry;
  Element_List **T = omega->T;
	Re 		 = 1.0/dparam("KINVIS");
  Je   	 = iparam("INTYPE");
	nspecs = iparam("NSPEC");

	for (int i = 0; i < nspecs; i++){
		bndry = Tbc[i];
  	while (bndry) {
  	  switch (bndry->type) {
			case 'O': case 'V': case 'W': case 'T': case 'R':
				break;
  	  case 'F': {
  	    CalcTbc(bndry, T[i]);
#ifdef CCHF
				if (bndry->usrtype == 'f')
					CalcTbc(bndry, T, Re*omega->Pr[i], i);
#endif
				break;
			}
	    default:
  	    error_msg(SetTBCs -- unknown transport b.c.)
  	    break;
  	  }
  	  bndry = bndry->next;
  	}
	}

  return;
}

static void CalcTbc(Bndry *B, Element_List *V){
  const    int id  = B->elmt->id;
  const 	 int face = B->face;
  char 				 *func_string = B->bstring; // function string
	int						qa, qb;
  Element  *E = B->elmt;
  Element  *t = V->flist[id];

  // get order of face quad 
  if (E->identify() == Nek_Tet || E->identify() == Nek_Hex){
	  qa   = E->qa;
  	qb   = E->qb;
 	}
  else if (E->identify() == Nek_Prism){
  	if(face == 1 || face == 3){
  		qa = E->qa;
  		qb = E->qc;
  	}
  	else{
  		qa   = E->qa;
  		qb   = E->qb;
  	}
  }
  else{
  	fprintf(stderr,"CalcTbc needs setting up for pyramids...\n");
  	exit(1);
  }

	Coord X;
  X.x = dvector(0, qa*qb - 1);
  X.y = dvector(0, qa*qb - 1);
  X.z = dvector(0, qa*qb - 1);
  double *func = dvector(0, QGmax*QGmax - 1);
  double *tmp = dvector(0, QGmax*QGmax - 1);

  E->GetFaceCoord(face, &X);

  vector_def("x y z", func_string); // evaluate flux
  vector_set(qa*qb, X.x, X.y, X.z, func); // flux values in func vector

  E->InterpToFace1(face, func, tmp);

  t->MakeFlux(B, face, tmp);

	free(X.x); free(X.y); free(X.z);
	free(func); free(tmp);
}
