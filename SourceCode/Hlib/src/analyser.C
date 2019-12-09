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
#include <string.h>
#include "nektar.h"

/*
 * Run-time Analyzer ... called every time step
 */

static int  History       (Domain *omega, double time);
static void hisHeader     (Domain *omega);
static void addfields     (Element_List *V[]);
static void replacefields (Element_List *V[], int nsteps);
static void avgfields (Element_List *V[], Element_List *Vb[], int nsteps);
static int check_number=0;

int init = 1, verbose, iostep, hisstep, nsteps, timeavg, startavg, initavg=1,step_iavg;
double    TSAVG_mod;
double    TEAVG_mod;
static double last_time=0.0;
double    time_start;	

double    TSAVG;
double    TEAVG;
double    Tp;

void Analyser (Domain *omega, int step, double time)
{
  FILE     *fp[2];
  FILE     *fpb[2];
	
  int      nfields = omega->U->fhead->dim()+2;

  int i;
  Element_List  **V;
  Element_List  **Vb;
	
  char      fname[FILENAME_MAX];
  double    step_length;
  


  dparam_set("t", time);

  V = (Element_List**) malloc(nfields*sizeof(Element_List*));

  /* ..........  Field Files   ......... */
  
  V[0]   = omega->U;
  V[1]   = omega->V;
  V[2]   = omega->W;
  V[3]   = omega->P;
  V[4]   = omega->T;	
	
	// Calculating turbulent statistics	
	// This part has to be moved to a separate function
	// Hessam- Jan  2012 
	Vb = (Element_List**) malloc(nfields*sizeof(Element_List*));
	Vb[0]   = omega->Ub;
	Vb[1]   = omega->Vb;
	Vb[2]   = omega->Wb;
	Vb[3]   = omega->Pb;
	Vb[4]   = omega->Tb;
	
	
	
	
	Vb[0]->Set_state('t');
	Vb[1]->Set_state('t');
	Vb[2]->Set_state('t');
	Vb[3]->Set_state('t');
	Vb[4]->Set_state('t');
	
	
	
	

  if(init){
    verbose   = option("verbose");
    iostep    = iparam("IOSTEP");
    hisstep   = option("hisstep");
    if(!hisstep) hisstep = iparam("HISSTEP");
    nsteps    = iparam("NSTEPS");
    timeavg   = option("timeavg");
    TSAVG = dparam("TSAVG");
    TEAVG = dparam("TEAVG");
    Tp = dparam("Tp");
    time_start = time;  
    if (omega->his_list) hisHeader (omega);
    init = 0;
  }

// Start time and end time for average are modified so they begin exactly at the
// begining of a pulsation period 
if ((time-time_start)>= TSAVG && (time/Tp - floor(time/Tp)  )<1e-8  ) {
	if(initavg){
		step_iavg = step;
		initavg = 0;
		TSAVG_mod = time;
		TEAVG_mod = time + (TEAVG-TSAVG);
		ROOTONLY fprintf(stdout, "Time average started at TSAVG = %g  and ends at TEVAG = %g\n", TSAVG_mod,TEAVG_mod);
		   fflush(stdout);
	}
	
				
}
if (initavg ==0 && time>=TSAVG_mod && time<TEAVG_mod){
         
	avgfields (V, Vb, step-step_iavg+1);
}

  /* .......... General Output ......... */
  step_length = (clock()-last_time)/(double)CLOCKS_PER_SEC;
  last_time = clock();
  

  ROOTONLY fprintf(stdout, "Time step = %d, Time = %g Cpu-Time = %g\n", 
		   step, time, step_length);
  fflush(stdout);
  if (step == 0){                         /* Everything else is for step > 0 */
    if (omega->his_list)                  /* Do initial history point        */
      History (omega, time);
    return;
  }
  if ((step % hisstep) == 0){
    if (omega->his_list){
      History (omega, time);
      fflush(omega->his_file);
    }
    
    forces(omega,step,time);
    ROOTONLY 
      fflush(omega->fce_file);
 
    /* flush stdout at step as well */
    ROOTONLY 
      fflush(stdout);
    
    
    if(option("SurForce")){
      cnt_srf_elements(omega);
    }
    
    if(verbose){
      if(omega->soln)
	for(i=0;i<nfields;++i)
	  V[i]->Terror(omega->soln[i]);
      else
	V[0]->Terror("0.0");
    }
  }


 if(option("SurForce")){

 int a = cnt_srf_elements(omega);
 }
 


  if (step % iostep == 0 && step < nsteps) {         
    if (option ("checkpt")) {
      DO_PARALLEL{
	if(option("SLICES")&&check_number<100){
	  sprintf (fname, "%s_%d.chk.hdr.%d", omega->name, check_number,
		   pllinfo.procid);
	  fp[0] = fopen(fname,"w");
	  
	  sprintf (fname, "%s_%d.chk.dat.%d", omega->name, check_number,
		   pllinfo.procid);
	  fp[1] = fopen(fname,"w");
	  ++check_number;
	}
	else{
	  sprintf (fname, "%s.chk.hdr.%d", omega->name, pllinfo.procid);
	  fp[0] = fopen(fname,"w");
	
	  sprintf (fname, "%s.chk.dat.%d", omega->name, pllinfo.procid);
	  fp[1] = fopen(fname,"w");
	}

	Writefld (fp, omega->name, step, time, nfields, V);
	fclose(fp[0]);
	fclose(fp[1]);      
      }
      else{
	if(option("SLICES")&&check_number<100){
	  sprintf (fname, "%s_%d.chk",  omega->name,check_number);
	  ++check_number;
	}
	else
	  sprintf (fname, "%s.chk", omega->name);
	fp[1] = fp[0] = fopen(fname,"w");
		  
	Writefld (fp, omega->name, step, time, nfields, V);
		  
	fclose(fp[0]);
      }
    }
  }
  
	
	if (step % iostep == 0 && step < nsteps && time>=TSAVG_mod && time<TEAVG_mod) {         
		if (option ("checkpt")) {
			DO_PARALLEL{
				if(option("SLICES")&&check_number<100){
					sprintf (fname, "stat%s_%d.chk.hdr.%d", omega->name, check_number,
							 pllinfo.procid);
					fp[0] = fopen(fname,"w");
					
					sprintf (fname, "stat%s_%d.chk.dat.%d", omega->name, check_number,
							 pllinfo.procid);
					fp[1] = fopen(fname,"w");
					++check_number;
				}
				else{
// 				  
					sprintf (fname, "%s_stat.chk.hdr.%d", omega->name, pllinfo.procid);
					fp[0] = fopen(fname,"w");
					
					sprintf (fname, "%s_stat.chk.dat.%d", omega->name, pllinfo.procid);
					fp[1] = fopen(fname,"w");
				}
				    	
				Writefld (fp, omega->name, step, time, nfields, Vb);
				fclose(fp[0]);
				fclose(fp[1]);     
			}
			else{
				if(option("SLICES")&&check_number<100){
					sprintf (fname, "stat%s_%d.chk",  omega->name,check_number);
					++check_number;
				}
				else
					sprintf (fname, "stat%s.chk", omega->name);
				fp[1] = fp[0] = fopen(fname,"w");
				
				Writefld (fp, omega->name, step, time, nfields, Vb);
				
				fclose(fp[0]);
			}
		}
	}
	
	/*
	if ( step < nsteps) {         
		if (option ("checkpt")) {
			DO_PARALLEL{
				if(option("SLICES")&&check_number<100){
					sprintf (fname, "%s_%d.chk.hdr.%d", omega->name, check_number,
							 pllinfo.procid);
					fp[0] = fopen(fname,"w");
					
					sprintf (fname, "%s_%d.chk.dat.%d", omega->name, check_number,
							 pllinfo.procid);
					fp[1] = fopen(fname,"w");
					++check_number;
				}
				else{
					sprintf (fname, "%s.chk.hdr.%d", omega->name, pllinfo.procid);
					fp[0] = fopen(fname,"w");
					
					sprintf (fname, "%s.chk.dat.%d", omega->name, pllinfo.procid);
					fp[1] = fopen(fname,"w");
				}
				
				Writefld (fp, omega->name, step, time, nfields, V);
				fclose(fp[0]);
				fclose(fp[1]);      
			}
			else{
				if(option("SLICES")&&check_number<100){
					sprintf (fname, "%s_%d.chk",  omega->name,check_number);
					++check_number;
				}
				else
				sprintf (fname, "%s.stat.chk", omega->name);
				fpb[1] = fpb[0] = fopen(fname,"w");
				
				Writefld (fpb, omega->name, step, time, nfields, Vb);
				
				fclose(fpb[0]);	  
			}
		}
	}
	*/
  if(timeavg){
    addfields (V);
    if(step == nsteps) replacefields(V,nsteps);
  }

  V[0]->Set_state('t');
  V[1]->Set_state('t');
  V[2]->Set_state('t');
  V[3]->Set_state('t');




	
	
	free(V);
	free(Vb);
	
  return;		
}
    

/* ------------------------------------------------------------------------- *
 * History() -- Process history points                                       *
 *                                                                           *
 * This function processes the history point list and transforms data points *
 * if necessary to physical space.                                           *
 * ------------------------------------------------------------------------- */
static int gatherPts     (HisPoint *hp, Element_List **V, double *vbuf[HP_MAX]);

static int History (Domain *omega, double time){
  FILE     *fp = omega->his_file;
  HisPoint *hp = omega->his_list;
  Element_List  *V [MAXDIM+1];
  double   *vbuf[HP_MAX];
  register  int n, cnt;

  if (!hp) return 0;

  V [0]   = omega->U; 
  V [1]   = omega->V;
  V [2]   = omega->W;
  V [3]   = omega->P;

  gatherPts (hp, V, vbuf);
  cnt     = 0;
  
  do 
    { fprintf (fp, "%lf ", time);
      for (n = 0; n < strlen(hp->flags); n++)
	fprintf (fp, "%#13.6g ", vbuf[cnt][n]);
      fprintf (fp, ":%d\n", cnt+1);
      free (vbuf[cnt++]);    } 
  while 
    (hp = hp->next);
  
  return cnt;
}

/* Collect history points from the Elements */
static double *modecenter(Element *E, Edge *e);

static int gatherPts (HisPoint *hp, Element_List **V, double *vbuf[HP_MAX])
{
  register  int i, j, n, pos;
  int cnt;

  cnt = 4;

  for (i = 0; hp ; ++i, hp = hp->next) {
    vbuf[i] = dvector (0, (int) strlen (hp->flags)-1);
    for  (n = pos = 0; n < cnt; ++n) {
      if (strchr (hp->flags, V[n]->fhead->type)) {
	switch (hp->mode) {
	case TransVert:
	  if(V[n]->fhead->state == 't')
	    vbuf[i][pos++] = V[n]->flist[hp->id]->vert[hp->i].hj[0];
	  break;
	case TransEdge:
	  if(V[n]->fhead->state == 't'){
	    Element *E = V[n]->flist[hp->id];
	    Edge   *e    = E->edge+hp->i;
	    double *mode = modecenter(E,e);
	    
	    vbuf[i][pos] =0.5*(E->vert[E->ednum(hp->i,0)].hj[0]
			       +E->vert[E->ednum(hp->i,1)].hj[0]);
	    
	    for(j = 0; j < e->l; ++j) vbuf[i][pos] += e->hj[j]*mode[j];
	    pos++;
	  }
	  break;
	case Physical:
	  break;
	default:
	  error_msg (History -- unknown history point mode);
	    break;
	}
      }
    }
  }
  
  return i;
}

/* find the center of modes  and store */
static double *modecenter(Element *E, Edge *edg){
  static double *mode;
  static int Lmode;

  if(!(mode&&(edg->l <= Lmode))){
    int i,qa = E->qa;
    Mode *e = E->getbasis()->edge[0];
    
    if(!mode) free(mode);
    mode = dvector(0,edg->l);

    if(qa%2 == 0){ /* if even spacing interpolate to center point */
      double **im;
      getim(qa,qa+1,&im,a2a);
      
      for(i = 0; i < edg->l; ++i)
	mode[i] = ddot(qa,im[qa/2],1,e[i].a,1);
    }      
    else           /* else use center value which is always at center point */
      for(i = 0; i < edg->l; ++i)
	mode[i] = e[i].a[qa/2];
  }

  return mode;
}

/* Write the header for the history point file */

static void hisHeader (Domain *omega)
{
  FILE      *fp = omega->his_file;
  HisPoint  *hp = omega->his_list;
  Element_List   *U  = omega->U;
  Element *E;
  int        n  = 1;

  if (!fp) return;

  fputs ("# Nektar history point file\n"
	 "# \n"
	 "# History points:\n", fp);

  do 
    {
      E = U->flist[hp->id];
      if(hp->mode == TransVert){
	fprintf (fp, "#  %3d: x = %#6.2lf, y = %#6.2lf, z = %#6.2lf,", n++, 
		 E->vert[hp->i].x,
		 E->vert[hp->i].y,
		 E->vert[hp->i].z);
      }
      else if(hp->mode == TransEdge){
	fprintf (fp, "#  %3d: x = %#6.2lf, y = %#6.2lf, z = %#6.2lf, ", n++,
		 0.5*(E->vert[E->edvnum(hp->i,0)].x + 
		      E->vert[E->edvnum(hp->i,1)].x),
		 0.5*(E->vert[E->edvnum(hp->i,0)].y + 
		      E->vert[E->edvnum(hp->i,1)].y),
		 0.5*(E->vert[E->edvnum(hp->i,0)].z + 
		      E->vert[E->edvnum(hp->i,1)].z));
      }	

      fprintf (fp, "fields = %4s, [%2d %3d]\n",
	       hp->flags,  hp->i+1, hp->id+1);  
    }
  while
    (hp = hp->next);
  fputs ("#\n", fp);
  return;
}

static double **avg;

static void addfields (Element_List *V[]){
  register int i;
  Element  *E;
  double   *s;
  int       eDIM = V[0]->fhead->dim();
  int       nf = eDIM+1;

  if(!avg){
    int ntot = 0;
    /* count transform storage */
    for(E = V[0]->fhead; E; E = E->next)  ntot += E->Nmodes;
    
    avg = dmatrix(0,eDIM,0,ntot-1);
    dzero (nf*ntot,avg[0],1);
  }

  for(i = 0; i < nf; ++i)
    for(E = V[i]->fhead, s = avg[i]; E; E = E->next){
      dvadd (E->Nmodes,E->vert->hj,1,s,1,s,1);
      s += E->Nmodes;
    }
}

static void replacefields (Element_List *V[], int nsteps){
  register int i;
  Element  *E;
  double   *s,fac;
  int       eDIM = V[0]->fhead->dim();
  int       nf = eDIM+1;

  fac = 1.0/(double)nsteps;
  
  for(i = 0; i < nf; ++i)
    for(E = V[i]->fhead, s = avg[i]; E; E = E->next){
      dsmul(E->Nmodes,fac,s,1,E->vert->hj,1);
      s += E->Nmodes;
    }

  free_dmatrix(avg,0,0);
}

static void avgfields (Element_List *V[], Element_List *Vb[], int step){
	register int i;
	Element  *E;
	Element  *Eb;
	double   *s,fac;
	int       eDIM = V[0]->fhead->dim();
	int       nf = eDIM+2;
	int ntot;
	
	fac = (double)step - 1.0;
	ntot = Vb[0]->hjtot;
	for(i = 0; i < nf; ++i){
			dsvtvp (ntot, fac, Vb[i]->base_hj, 1, V[i]->base_hj, 1, Vb[i]->base_hj, 1);
	}
	
	fac = 1.0/(double)step;
	
	for(i = 0; i < nf; ++i){
			dsmul(ntot,fac,Vb[i]->base_hj,1,Vb[i]->base_hj,1);
		}
	
}
