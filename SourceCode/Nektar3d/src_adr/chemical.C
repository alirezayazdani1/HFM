/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $Source: 
 * $Revision:
 * $Date: 
 * $Author:
 * $State: 
 *---------------------------------------------------------------------------*/
#include <stdio.h>
#include "nektar.h"

#ifdef ADR
#if defined (CCHF) || (CCLF)

static int 		nspecs, Nquad;
static double TFVIIaW, Linj;
static void 	ChkFld(Domain *omega);

#ifdef CCLF
	static double indicator(double d);
	static double bind_affinity(double eta);
	static void 	GetParticleNumber(Domain *omega, Element *E, int elem);
	static double	GetAgonistSource(double **timeAct, int num, int elem);
	static double PLT0, visc, alphamax, dt;
#endif

void AddCCSource(Domain *omega){
	nspecs = iparam("NSPEC");

  Element **E = (Element **) malloc(nspecs*sizeof(Element *));
  Element **Ef = (Element **) malloc(nspecs*sizeof(Element *));

  //Alireza: high-fidelity coagulation cascade 
#ifdef CCHF
	ChkFld( omega );

  double *den = dvector(0, QGmax*QGmax*QGmax - 1);
  double *num = dvector(0, QGmax*QGmax*QGmax - 1);
  double *first = dvector(0, QGmax*QGmax*QGmax - 1);
  double *second = dvector(0, QGmax*QGmax*QGmax - 1);
  double *third = dvector(0, QGmax*QGmax*QGmax - 1);
  double *tenase = dvector(0, QGmax*QGmax*QGmax - 1);
  double *proth = dvector(0, QGmax*QGmax*QGmax - 1);
  double **sources = dmatrix(0, nspecs-1, 0, QGmax*QGmax*QGmax-1);
  for (int i = 0; i < omega->U->nel; ++i){
		dzero(QGmax*QGmax*QGmax, den, 1);
		dzero(QGmax*QGmax*QGmax, num, 1);
		dzero(QGmax*QGmax*QGmax, first, 1);
		dzero(QGmax*QGmax*QGmax, second, 1);
		dzero(QGmax*QGmax*QGmax, third, 1);
		dzero(QGmax*QGmax*QGmax, tenase, 1);
		dzero(QGmax*QGmax*QGmax, proth, 1);
    for (int k = 0; k < nspecs; k++){
      E[k] = omega->T[k]->flist[i];
      Ef[k] = omega->Tf[k]->flist[i];
			dzero(QGmax*QGmax*QGmax, &sources[k][0], 1);
		}
    int qtot = E[0]->qtot;

    dsadd(qtot, K9M, E[1]->h_3d[0][0], 1, den, 1);
    dvmul(qtot, E[12]->h_3d[0][0], 1, E[1]->h_3d[0][0], 1, num, 1);
    dsmul(qtot, k9, num, 1, num, 1);
    dvdiv(qtot, num, 1, den, 1, first, 1);
    dvmul(qtot, E[14]->h_3d[0][0], 1, E[0]->h_3d[0][0], 1, second, 1);
    dsmul(qtot, -h9, second, 1, second, 1);
    dvadd(qtot, first, 1, second, 1, &sources[0][0], 1); // IXa

    dsmul(qtot, -1.0, first, 1, &sources[1][0], 1); // IX

    dsadd(qtot, K8M, E[3]->h_3d[0][0], 1, den, 1);
    dvmul(qtot, E[8]->h_3d[0][0], 1, E[3]->h_3d[0][0], 1, num, 1);
    dsmul(qtot, k8, num, 1, num, 1);
    dvdiv(qtot, num, 1, den, 1, first, 1);
    dsmul(qtot, -h8, E[2]->h_3d[0][0], 1, second, 1);
    dsadd(qtot, HC8M, E[2]->h_3d[0][0], 1, den, 1);
    dvmul(qtot, E[16]->h_3d[0][0], 1, E[2]->h_3d[0][0], 1, num, 1);
    dsmul(qtot, -hC8, num, 1, num, 1);
    dvdiv(qtot, num, 1, den, 1, third, 1);
    dvadd(qtot, first, 1, second, 1, first, 1);
    dvadd(qtot, first, 1, third, 1, &sources[2][0], 1); // VIIIa

    dsmul(qtot, -1.0, first, 1, &sources[3][0], 1); // VIII

    dsadd(qtot, K5M, E[5]->h_3d[0][0], 1, den, 1);
    dvmul(qtot, E[8]->h_3d[0][0], 1, E[5]->h_3d[0][0], 1, num, 1);
    dsmul(qtot, k5, num, 1, num, 1);
    dvdiv(qtot, num, 1, den, 1, first, 1);
    dsmul(qtot, -h5, E[4]->h_3d[0][0], 1, second, 1);
    dsadd(qtot, HC5M, E[4]->h_3d[0][0], 1, den, 1);
    dvmul(qtot, E[16]->h_3d[0][0], 1, E[4]->h_3d[0][0], 1, num, 1);
    dsmul(qtot, -hC5, num, 1, num, 1);
    dvdiv(qtot, num, 1, den, 1, third, 1);
    dvadd(qtot, first, 1, second, 1, first, 1);
    dvadd(qtot, first, 1, third, 1, &sources[4][0], 1); // Va

    dsmul(qtot, -1.0, first, 1, &sources[5][0], 1); // V

    dvmul(qtot, E[0]->h_3d[0][0], 1, E[2]->h_3d[0][0], 1, tenase, 1);
    dsmul(qtot, 1.0/KdZ, tenase, 1, tenase, 1); // Z

    dsadd(qtot, K10M, E[7]->h_3d[0][0], 1, den, 1);
    dvmul(qtot, tenase, 1, E[7]->h_3d[0][0], 1, num, 1);
    dsmul(qtot, k10, num, 1, num, 1);
    dvdiv(qtot, num, 1, den, 1, first, 1);
    dvmul(qtot, E[6]->h_3d[0][0], 1, E[14]->h_3d[0][0], 1, second, 1);
    dsmul(qtot, -h10, second, 1, second, 1);
    dvmul(qtot, E[15]->h_3d[0][0], 1, E[6]->h_3d[0][0], 1, third, 1);
    dsmul(qtot, -hTFPI, third, 1, third, 1);
    dvadd(qtot, first, 1, second, 1, first, 1);
    dvadd(qtot, first, 1, third, 1, &sources[6][0], 1); // Xa

    dsmul(qtot, -1.0, first, 1, &sources[7][0], 1); // X

    dvmul(qtot, E[4]->h_3d[0][0], 1, E[6]->h_3d[0][0], 1, proth, 1);
    dsmul(qtot, 1.0/KdW, proth, 1, proth, 1); // W

    dsadd(qtot, K2M, E[9]->h_3d[0][0], 1, den, 1);
    dvmul(qtot, proth, 1, E[9]->h_3d[0][0], 1, num, 1);
    dsmul(qtot, k2, num, 1, num, 1);
    dvdiv(qtot, num, 1, den, 1, first, 1);
    dvmul(qtot, E[8]->h_3d[0][0], 1, E[14]->h_3d[0][0], 1, second, 1);
    dsmul(qtot, -h2, second, 1, second, 1);
    dvadd(qtot, first, 1, second, 1, &sources[8][0], 1); // IIa

    dsmul(qtot, -1.0, first, 1, &sources[9][0], 1); // II

    dsadd(qtot, K1M, E[11]->h_3d[0][0], 1, den, 1);
    dvmul(qtot, E[8]->h_3d[0][0], 1, E[11]->h_3d[0][0], 1, num, 1);
    dsmul(qtot, k1, num, 1, num, 1);
    dvdiv(qtot, num, 1, den, 1, &sources[10][0], 1); // Ia

    dsmul(qtot, -1.0, first, 1, &sources[11][0], 1); // I

    dsadd(qtot, K11M, E[13]->h_3d[0][0], 1, den, 1);
    dvmul(qtot, E[8]->h_3d[0][0], 1, E[13]->h_3d[0][0], 1, num, 1);
    dsmul(qtot, k11, num, 1, num, 1);
    dvdiv(qtot, num, 1, den, 1, first, 1);
    dvmul(qtot, E[12]->h_3d[0][0], 1, E[14]->h_3d[0][0], 1, second, 1);
    dsmul(qtot, -h11A3, second, 1, second, 1);
    dvmul(qtot, E[12]->h_3d[0][0], 1, E[18]->h_3d[0][0], 1, third, 1);
    dsmul(qtot, -h11L1, third, 1, third, 1);
    dvadd(qtot, first, 1, second, 1, first, 1);
    dvadd(qtot, first, 1, third, 1, &sources[12][0], 1); // XIa

    dsmul(qtot, -1.0, first, 1, &sources[13][0], 1); // XI

    dsmul(qtot, h9, E[0]->h_3d[0][0], 1, first, 1);
    dsmul(qtot, h10, E[6]->h_3d[0][0], 1, second, 1);
    dvadd(qtot, first, 1, second, 1, first, 1);
    dsmul(qtot, h2, E[8]->h_3d[0][0], 1, second, 1);
    dvadd(qtot, first, 1, second, 1, first, 1);
    dsmul(qtot, h11A3, E[12]->h_3d[0][0], 1, second, 1);
    dvvpvt(qtot, first, 1, second, 1, E[14]->h_3d[0][0], 1, third, 1);
    dsmul(qtot, -1.0, third, 1, &sources[14][0], 1); // ATIII

    dvmul(qtot, E[15]->h_3d[0][0], 1, E[6]->h_3d[0][0], 1, first, 1);
    dsmul(qtot, -hTFPI, first, 1, &sources[15][0], 1); // TFPI

    dsadd(qtot, KPCM, E[17]->h_3d[0][0], 1, den, 1);
    dvmul(qtot, E[8]->h_3d[0][0], 1, E[17]->h_3d[0][0], 1, num, 1);
    dsmul(qtot, kPC, num, 1, num, 1);
    dvdiv(qtot, num, 1, den, 1, first, 1);
    dvmul(qtot, E[16]->h_3d[0][0], 1, E[18]->h_3d[0][0], 1, second, 1);
    dsmul(qtot, -hPC, second, 1, second, 1);
    dvadd(qtot, first, 1, second, 1, &sources[16][0], 1); // APC

    dsmul(qtot, -1.0, first, 1, &sources[17][0], 1); // PC

    dvmul(qtot, E[16]->h_3d[0][0], 1, E[18]->h_3d[0][0], 1, first, 1);
    dsmul(qtot, -hPC, first, 1, first, 1);
    dvmul(qtot, E[12]->h_3d[0][0], 1, E[18]->h_3d[0][0], 1, second, 1);
    dsmul(qtot, -h11L1, second, 1, second, 1);
    dvadd(qtot, first, 1, second, 1, &sources[18][0], 1); // L1AT

    for (int k = 0; k < nspecs; k++)
      dvadd(qtot, Ef[k]->h_3d[0][0], 1, &sources[k][0], 1, Ef[k]->h_3d[0][0], 1);
  }
  free(den); free(num); free(first); free(second); free(third); free(tenase); free(proth); free_dmatrix(sources,0,0);
#endif

//Alireza: low-fidelity coagulation cascade 
#ifdef CCLF
	ChkFld( omega );

  double *first = dvector(0, QGmax*QGmax*QGmax - 1);
  double *second = dvector(0, QGmax*QGmax*QGmax - 1);
  double *third = dvector(0, QGmax*QGmax*QGmax - 1);
  double **sources = dmatrix(0, nspecs-1, 0, QGmax*QGmax*QGmax-1);
	int offset = 0, Nagon = agonists[0];
	double phi, convol;
	SPECIES spec;
	PLT0 = dparam("PLT0");
	dt	 = dparam("DELT");
  for (int i = 0; i < omega->U->nel; ++i){
    dzero(QGmax*QGmax*QGmax, first, 1);
    dzero(QGmax*QGmax*QGmax, second, 1);
    dzero(QGmax*QGmax*QGmax, third, 1);
    for (int k = 0; k < nspecs; k++){
      E[k] = omega->T[k]->flist[i];
      Ef[k] = omega->Tf[k]->flist[i];
      dzero(QGmax*QGmax*QGmax, &sources[k][0], 1);
    }
    int qtot = E[0]->qtot;

  	for (int j = 0; j < Nagon; ++j){
    	spec = static_cast<SPECIES>(agonists[j+1]);
    	daxpy(qtot, weight[j]/thresh[j], E[spec]->h_3d[0][0], 1, first, 1);
  	}
		for (int j = 0; j < qtot; ++j)
			if (first[j] < 1.0) first[j] = 0.0;
			else	first[j] /= -tact;
		dvmul(qtot, first, 1, E[0]->h_3d[0][0], 1, &sources[0][0], 1);	// PLTmu

		dsmul(qtot, -1.0, &sources[0][0], 1, &sources[1][0], 1);
		dsadd(qtot, -PLT0, E[7]->h_3d[0][0], 1, first, 1);
		for (int j = 0; j < qtot; ++j){
			omega->ind_func[j+offset] = phi = indicator( first[j] );
			second[j] = bind_affinity( phi );
		}
		offset += qtot;
		dsmul(qtot, -kcohxPLTmax, second, 1, second, 1);
		dvmul(qtot, second, 1, E[1]->h_3d[0][0], 1, second, 1);
		dvadd(qtot, second, 1, &sources[1][0], 1, &sources[1][0], 1);
		
		dsmul(qtot, -1.0, second, 1, second, 1);
		dvadd(qtot, second, 1, &sources[7][0], 1, &sources[7][0], 1);  // PLTba

		dsvvpt(qtot, kiiAP, E[1]->h_3d[0][0], 1, E[7]->h_3d[0][0], 1, first, 1);
		dsadd(qtot, ksurf, first, 1, first, 1);
		dvmul(qtot, first, 1, E[3]->h_3d[0][0], 1, first, 1);
		dvadd(qtot, first, 1, &sources[2][0], 1, &sources[2][0], 1);
		dsmul(qtot, -kin, E[2]->h_3d[0][0], 1, second, 1);
		dvadd(qtot, second, 1, &sources[2][0], 1, &sources[2][0], 1);	// IIa

		dsmul(qtot, -1.0, first, 1, first, 1);
		dvadd(qtot, first, 1, &sources[3][0], 1, &sources[3][0], 1);	// II

		GetParticleNumber( omega, E[7], i );

		if (omega->number[i] != 0){
			convol = GetAgonistSource( omega->timeAct, omega->number[i], i );
			convol *= PLTmax * omega->volume[i] / (double) Ndis;

			dsadd(qtot, convol*relCon[1]/(double) qtot, &sources[4][0], 1, &sources[4][0], 1);	// ADP
			
			dsadd(qtot, convol*relCon[2]/(double) qtot, &sources[5][0], 1, &sources[5][0], 1);	// TxA2
		}

		daxpy(qtot, k1, E[2]->h_3d[0][0], 1, &sources[6][0], 1);	// Ia

    for (int k = 0; k < nspecs; k++)
      dvadd(qtot, Ef[k]->h_3d[0][0], 1, &sources[k][0], 1, Ef[k]->h_3d[0][0], 1);
	}
  free(first); free(second); free(third); free_dmatrix(sources,0,0);
#endif

	return;
}

static void ChkFld(Domain *omega){
	Nquad = omega->T[0]->htot;

	for (int i = 0; i < nspecs; ++i)
		for (int j = 0; j < Nquad; ++j){
			if (omega->T[i]->base_h[j] < 0.0) omega->T[i]->base_h[j] = 0.0;
#ifdef CCLF
			if ( (i == 0 || i == 1 || i == 7) && omega->T[i]->base_h[j] > PLTmax) omega->T[i]->base_h[j] = PLTmax;
#endif
		}
}

#ifdef CCLF
void AddBrinkman(Domain *omega){

  Element_List	*Uf   =  omega->Uf, *Vf   =  omega->Vf,  *Wf  = omega->Wf,
  							*U    =  omega->U,  *V		=  omega->V,   *W		= omega->W,
  							*T    =  omega->T[7];	// PLTba
  int Nquad = U->htot;
	double *num = dvector(0, Nquad - 1);
	double *den = dvector(0, Nquad - 1);
	double *coeff = dvector(0, Nquad - 1);
	visc 		 = dparam("KINVIS");
	alphamax = dparam("AMAX");

	dsmul(Nquad, 1.0/PLTmax, T->base_h, 1, den, 1);
	dvmul(Nquad, den, 1, den, 1, num, 1);
	dsadd(Nquad, phi0B*phi0B, num, 1, den, 1);
	dsmul(Nquad, -alphamax*visc, num, 1, num, 1);
	dvdiv(Nquad, num, 1, den, 1, coeff, 1);
	
	dvmul(Nquad, coeff, 1, U->base_h, 1, Uf->base_h, 1);
	dvmul(Nquad, coeff, 1, V->base_h, 1, Vf->base_h, 1);
	dvmul(Nquad, coeff, 1, W->base_h, 1, Wf->base_h, 1);

	free(num); free(den); free(coeff);
	return;
}

static double indicator(double d){
  return 0.5 * ( tanh(-d/ksi) + 1.0 );
}

static double bind_affinity(double eta){
	double dum = (eta-etat);
	double exp1 = dum*dum*dum;
	return g0 * std::max( 0.0, exp1/(etastar+exp1) );
}

static void GetParticleNumber(Domain *omega, Element *E, int elem){
  int i, j, qa = E->qa, qb = E->qb, qc = E->qc, qab = qa*qb, qbc = qb*qc, qtot = E->qtot;
	double sum, time = dparam("t");
  double *wabc, *z, *wa, *wb, *wc;
  wabc = dvector(0, QGmax*QGmax*QGmax - 1);

  E->GetZW(&z, &wa, &z, &wb, &z, &wc);
  for(i = 0; i < qbc; ++i)
  	dcopy(qa,wa,1,wabc+i*qa,1);
  for(i = 0; i < qc; ++i)
  	for(j = 0; j < qa; ++j)
    	dvmul(qb,wb,1,wabc+i*qab+j,qa,wabc+i*qab+j,qa);
  for(i = 0; i < qab; ++i)
    dvmul(qc,wc,1,wabc+i,qab,wabc+i,qab);
  if(E->curvX)
  	dvmul(qtot, E->geom->jac.p, 1, wabc, 1, wabc, 1);
  else
    dsmul(qtot, E->geom->jac.d, wabc, 1, wabc, 1);

	if (omega->volume[elem] == 0.0) omega->volume[elem] = dsum(qtot, wabc, 1);
 	dvmul(qtot, wabc, 1, E->h_3d[0][0], 1, wabc, 1);
	sum = dsum (qtot, wabc, 1);
	int num = (int) (sum * (double) Ndis / PLTmax );

	if ( (num  > omega->number[elem]) && (num <= Ndis) ){
		for (i = omega->number[elem]+1; i <= num; ++i) omega->timeAct[elem][i] = time;
		omega->number[elem] = num;
	}

	free(wabc);
	return;
}

static double GetAgonistSource(double **timeAct, int num, int elem){
	double time = dparam("t");
	double released = 0.0;

	for (int i = 1; i <= num; ++i)
    released += exp( -0.5 * (time-timeAct[elem][i]-mu)*(time-timeAct[elem][i]-mu) / (stddev*stddev) );

	return released;
}

#endif

#ifdef CCHF
void CalcTbc(Bndry *B, Element_List **V, double Pec, int index){
  const    int id  = B->elmt->id;
  const    int face = B->face;
  int     qa, qb;
  double *wk = dvector(0, QGmax*QGmax - 1);
  double *tmp = dvector(0, QGmax*QGmax - 1);
  double *den = dvector(0, QGmax*QGmax - 1);
  double *num = dvector(0, QGmax*QGmax - 1);
  Element  *E = B->elmt;
  Element  *t = V[index]->flist[id];
  TFVIIaW = dparam("TFVIIaW");
	Linj = 1.0;

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

  switch (index){
    case 0: {
      Element *eT = V[1]->flist[id];
      eT->GetFace(eT->h_3d[0][0], face, wk);

      dsadd(qa*qb, K79M, wk, 1, den, 1);
      dsmul(qa*qb, k79*TFVIIaW*(Linj*Pec), wk, 1, num, 1);
      dvdiv(qa*qb, num, 1, den, 1, wk, 1);

      eT->InterpToFace1(face, wk, tmp);
      break;
    }
    case 1: {
      Element *eT = V[1]->flist[id];
      eT->GetFace(eT->h_3d[0][0], face, wk);

      dsadd(qa*qb, K79M, wk, 1, den, 1);
      dsmul(qa*qb, -k79*TFVIIaW*(Linj*Pec), wk, 1, num, 1);
      dvdiv(qa*qb, num, 1, den, 1, wk, 1);

      eT->InterpToFace1(face, wk, tmp);
      break;
    }
    case 6: {
      Element *eT = V[7]->flist[id];
      eT->GetFace(eT->h_3d[0][0], face, wk);

      dsadd(qa*qb, K710M, wk, 1, den, 1);
      dsmul(qa*qb, k710*TFVIIaW*(Linj*Pec), wk, 1, num, 1);
      dvdiv(qa*qb, num, 1, den, 1, wk, 1);

      eT->InterpToFace1(face, wk, tmp);
      break;
    }
    case 7: {
      Element *eT = V[7]->flist[id];
      eT->GetFace(eT->h_3d[0][0], face, wk);

      dsadd(qa*qb, K710M, wk, 1, den, 1);
      dsmul(qa*qb, -k710*TFVIIaW*(Linj*Pec), wk, 1, num, 1);
      dvdiv(qa*qb, num, 1, den, 1, wk, 1);

      eT->InterpToFace1(face, wk, tmp);
      break;
    }
    default:
      dzero(qa*qb, tmp, 1);
      break;
  }
  t->MakeFlux(B, face, tmp);

  free(wk); free(tmp); free(den); free(num);
	return;
}
#endif

#endif
#endif
