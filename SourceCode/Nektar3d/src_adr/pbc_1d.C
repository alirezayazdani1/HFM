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
#include "nekstruct.h"
#include "pbc_1d.h"
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <veclib.h>


#ifdef PBC_1D

/* external function */
double BoundaryArea(Bndry *Bc, char *label);
double BoundaryAreaInlet(Bndry *Bc, char *label);
void M2Ptransform(int Nmodes, double om, my_dcmplx *MODES, int N, double Time_period, double *f);
void M2Ptransform(int Nmodes, double om, my_dcmplx* MODES , int N, double Time_period, int index_start, int index_end, double* f);
void FilterDFT(int N, int Nmodes, double  *f, double *Acos_sin, int FLAG_averag);
void Convolve(int N, double* f1, double* f2, double* ans, int time_index);
void ParallelConvolve(int N, double* f1, double* f2, double* ans, int it);
void ParallelConvolve(int N, double* f1, double* f2, double* ans, int index_start, int index_end, int it);
void ParallelConvolve_test(int N, double* f1, double* f2, double* ans, int index_start, int index_end, int it);

#define N_MAX_OUTLET 26 

int PBC1D::Create(Bndry *Ubc){

  int i,j;
  Bndry *Bc;
  int LABEL_FLAGS[N_MAX_OUTLET];
  memset(&LABEL_FLAGS[0],'\0',N_MAX_OUTLET*sizeof(int));

  standard_labels[0] = "bcA(x,y,z)\n";
  standard_labels[1] = "bcB(x,y,z)\n";
  standard_labels[2] = "bcC(x,y,z)\n";
  standard_labels[3] = "bcD(x,y,z)\n";
  standard_labels[4] = "bcE(x,y,z)\n";
  standard_labels[5] = "bcF(x,y,z)\n";
  standard_labels[6] = "bcG(x,y,z)\n";
  standard_labels[7] = "bcH(x,y,z)\n";
  standard_labels[8] = "bcI(x,y,z)\n";
  standard_labels[9] = "bcJ(x,y,z)\n";
  standard_labels[10] = "bcK(x,y,z)\n";
  standard_labels[11] = "bcL(x,y,z)\n";
  standard_labels[12] = "bcM(x,y,z)\n";
  standard_labels[13] = "bcN(x,y,z)\n";
  standard_labels[14] = "bcO(x,y,z)\n";
  standard_labels[15] = "bcP(x,y,z)\n";
  standard_labels[16] = "bcQ(x,y,z)\n";
  standard_labels[17] = "bcR(x,y,z)\n";
  standard_labels[18] = "bcS(x,y,z)\n";
  standard_labels[19] = "bcT(x,y,z)\n";
  standard_labels[20] = "bcU(x,y,z)\n";
  standard_labels[21] = "bcV(x,y,z)\n";
  standard_labels[22] = "bcW(x,y,z)\n";
  standard_labels[23] = "bcX(x,y,z)\n";
  standard_labels[24] = "bcY(x,y,z)\n";
  standard_labels[25] = "bcZ(x,y,z)\n";

  Nout = iparam("NOUTLETS");
  Ninl = iparam("NINLETS");
  

  /* STEP 1 check for consistency */

 /* create vector of integers "LABEL_FLAGS" 
    initialize with zero 
    loop over Bondaries look for ALL labels
    if label of outlet was found set corresponding flag to 1;
    sum results from all pertitions if parallel
    if Nout != sum (LABEL_FLAGS) -> setup is incosistent
       error message and exit.
    else
      deallocate  "LABEL_FLAGS" and continue setup
   */
  
   for (i = 0; i < N_MAX_OUTLET; ++i){

     for(Bc=Ubc;Bc;Bc=Bc->next){
       if ( (Bc->type == 'O' || Bc->type == 'o') && strcmp(Bc->blabel,standard_labels[i]) == 0 ){
          LABEL_FLAGS[i] = 1;
          break;
       }
     }
   }

   DO_PARALLEL{
      int *iwork;
      iwork = ivector(0,N_MAX_OUTLET-1);
      memset(iwork,'\0',N_MAX_OUTLET*sizeof(int));
      gisum (LABEL_FLAGS, N_MAX_OUTLET, iwork);
      free(iwork);
   }

   j = 0;
   for (i = 0; i < N_MAX_OUTLET; ++i){
     if (LABEL_FLAGS[i] != 0)
       j++;
   }
   if (j != Nout){
     fprintf(stderr, "Error in setup outlet BC. NOUTLETS in REA file = %d actual outlet count %d \n ",Nout,j);
     exit(-1);
   }

   for (i = 0; i < Nout; ++i){
     if (LABEL_FLAGS[i] == 0){
        fprintf(stderr, "Error in setup outlet BC. NOUTLETS in REA file = %d, outlet No. %d has no label \n",Nout,i);
        exit(-1);
     }
   }

   if (Nout == 0) 
    return 0;


   /* STEP 2 define parameters and allocate memory*/
  double dt = dparam("DELT");
  double t  = dparam("t");

  Nimpedance_modes  =  iparam("NIMPEDMODES");  
  Tperiod           =  dparam("DTIMEPERIOD");
 
  if (Tperiod == 0)
    Tperiod = dt*10.0;

  Nsteps_per_cycle = int ((Tperiod + dt*0.1)/dt)+1;

  j = numnodes();
  if ((Nsteps_per_cycle/j) > 100 && j > 1)
     parallel_convolution = 1;
  else
     parallel_convolution = 0;

  parallel_convolution = 0;

  /* for parallel convolution split "impedance" in Nproc. parts */

   if (parallel_convolution){

     Imp_length_local = (int) (Nsteps_per_cycle/j);
     i = Nsteps_per_cycle - Imp_length_local*j;
     if (mynode() < i)
       Imp_length_local++;

    /* to partition the global vector "impedance" find 
        first index corresponding to local processor */
  
     int *itemp;
     itemp = new int[j*2];
     memset(itemp,'\0',j*2*sizeof(int));
     DO_PARALLEL{
       itemp[mynode()] = Imp_length_local; 
       gisum(itemp,j,&itemp[j]);
     }

     Imp_index_start = 0;
     for (i = 1; i <= mynode(); i++) 
       Imp_index_start += itemp[i-1];

     free(itemp);
  }
  t = t-((int)((t+dt*0.1)/Tperiod))*Tperiod;
  time_step_in_cycle = int ((t+0.1*dt)/dt);

  FR_history = dmatrix(0,Nout-1,0,Nsteps_per_cycle-1);
  if (!parallel_convolution)
    impedance  = dmatrix(0,Nout-1,0,Nsteps_per_cycle-1);
  else
    impedance  = dmatrix(0,Nout-1,0,Imp_length_local);
  
  A_cos_sin  = dmatrix(0,Nout-1,0,Nimpedance_modes*2);
  impedance_modes = new my_dcmplx*[Nout];
  for (i = 0; i < Nout; ++i)
   impedance_modes[i] = new my_dcmplx[Nimpedance_modes*2+1];

  Pressure   = dvector(0,Nout+Ninl);
  Area       = dvector(0,Nout+Ninl);
  memset(Pressure,'\0',(Nout+Ninl)*sizeof(double));
  memset(Area,    '\0',(Nout+Ninl)*sizeof(double));

  for (i = 0; i < Nout; ++i){
    memset(FR_history[i],'\0',Nsteps_per_cycle*sizeof(double));
    if (!parallel_convolution)
      memset(impedance[i], '\0',Nsteps_per_cycle*sizeof(double));
    else
      memset(impedance[i], '\0',(Imp_length_local+1)*sizeof(double));
    memset(A_cos_sin[i], '\0',(Nimpedance_modes*2+1)*sizeof(double));
  }


  /* STEP 3 */
  /*  compute impedance */

  double rad,r_temp; 
  
  for (i = 0; i < Nout; ++i){
     rad =  BoundaryArea(Ubc, standard_labels[i]);

     DO_PARALLEL{
       gdsum(&rad,1,&r_temp);
     }
     Area[i] = rad;
     ROOTONLY
      printf("Area[outlet No. %d]  = %f \n",i,Area[i]);

     /*  radius must be passed in dimensional units! [cm] !!!! */ 
     rad = sqrt(rad/M_PI)*0.1;
     {
       SMALL_VESSEL root_vessel(rad);
       /* set omega  */
       ROOTONLY fprintf(stdout,"Ws = %2.16f \n",root_vessel.Ws);    
       int index = Nimpedance_modes+1;
       for (j = 1; j <= Nimpedance_modes; j++){
           impedance_modes[i][index] = root_vessel.GetZ0(j);
           impedance_modes[i][Nimpedance_modes-j].my_dcmplx_conj(impedance_modes[0][index]);
           ROOTONLY
             fprintf(stdout,"rank %d: impedance_modes[out = %d][mode = %d] = (%e,%e) \n",mynode(),
                             i,index,impedance_modes[i][index].real,impedance_modes[i][index].imag);
           index++;
       }
       impedance_modes[i][Nimpedance_modes] = root_vessel.GetZ0(0);
       ROOTONLY
         fprintf(stdout,"rank %d:impedance_modes[out = %d][mode = %d] = (%e,%e) \n",mynode(),
                               i,Nimpedance_modes,impedance_modes[i][Nimpedance_modes].real,
                                              impedance_modes[i][Nimpedance_modes].imag);

       ROOTONLY
         printf("SMALL_VESSEL - done, r = %f \n", rad );

      if (!parallel_convolution)
         M2Ptransform(Nimpedance_modes,omega_small_atree,impedance_modes[i],
                      Nsteps_per_cycle,2.0*M_PI/omega_small_atree,impedance[i]);
      else
         M2Ptransform(Nimpedance_modes,omega_small_atree,impedance_modes[i],
                      Nsteps_per_cycle,2.0*M_PI/omega_small_atree,
                      Imp_index_start,(Imp_index_start+Imp_length_local),impedance[i]);
    }
  }

  if (mynode() == (numnodes()-1))
    Imp_length_local--;

  /* compute area at inlet */
  for (i = 0; i < Ninl; ++i){
    rad =  BoundaryAreaInlet(Ubc, standard_labels[i]);
    r_temp = 0.0;
    DO_PARALLEL
     gdsum(&rad,1,&r_temp);
    Area[Nout+i] = rad;
    ROOTONLY
      printf("Area - inlet[%d] = %f \n",i,Area[Nout+i]);
  }

  /* compute number of faces with type='O' for each outlet in each partition */
  Nfaces_per_outlet = new int[Nout];
  Nfaces_per_inlet  = new int[Ninl];
  int face_index;

  for (i = 0; i < Nout; ++i){
     Nfaces_per_outlet[i] = 0;
     for(Bc=Ubc;Bc;Bc=Bc->next){
       if ( (Bc->type == 'O' || Bc->type == 'o') && strcmp(Bc->blabel,standard_labels[i]) == 0 )
          Nfaces_per_outlet[i]++;
     }
   } 
  

  for (i = 0; i < Ninl; ++i){
     Nfaces_per_inlet[i] = 0;
     for(Bc=Ubc;Bc;Bc=Bc->next){
       if ( (Bc->type == 'V' || Bc->type == 'v') && strcmp(Bc->blabel,standard_labels[i]) == 0 )
          Nfaces_per_inlet[i]++;
     }
   }

 
   /*  for each outlet get store ID's of faces with type == 'O' */
   ID_faces_per_outlet = new int*[Nout];
   for (i = 0; i < Nout; ++i)
     ID_faces_per_outlet[i] = new int[Nfaces_per_outlet[i]]; 

   for (i = 0; i < Nout; ++i){
     j = 0;
     face_index = 0;
     for(Bc=Ubc;Bc;Bc=Bc->next){
       if ( (Bc->type == 'O' || Bc->type == 'o') && strcmp(Bc->blabel,standard_labels[i]) == 0 ){
          ID_faces_per_outlet[i][face_index] = j;
          face_index++;
       }
       j++;
     }
   }

   ID_faces_per_inlet = new int*[Ninl];
   for (i = 0; i < Ninl; ++i)
     ID_faces_per_inlet[i] = new int[Nfaces_per_inlet[i]];

   for (i = 0; i < Ninl; ++i){
     j = 0;
     face_index = 0;
     for(Bc=Ubc;Bc;Bc=Bc->next){
       if ( (Bc->type == 'V' || Bc->type == 'v') && strcmp(Bc->blabel,standard_labels[i]) == 0 ){
          ID_faces_per_inlet[i][face_index] = j;
          face_index++;
       }
       j++;
     }
   }

   Nnodes_inlet  = new int[Ninl];
   Nnodes_outlet = new int[Nout];

#ifdef PARALLEL
   /* create communicator for outlets */
//   create_comm_BC(Nout, Nfaces_per_outlet);
   create_comm_BC_inlet_outlet(Ninl, Nfaces_per_inlet, Nnodes_inlet, Nout, Nfaces_per_outlet, Nnodes_outlet);
#endif

  return 0;
}

void PBC1D::SetGeofac(Bndry *Ubc, Bndry *Vbc, Bndry *Wbc){
  
   int index;
   Bndry *BcU, *BcV, *BcW; 
  
   for (BcU=Ubc,BcV=Vbc,BcW=Wbc;BcU;BcU = BcU->next,BcV = BcV->next,BcW = BcW->next){
     if ( BcU->type == 'o' || BcU->type == 'O' || BcU->type == 'v' || BcU->type == 'V' ){
               BcU->elmt->Surface_geofac(BcU);
               BcV->elmt->Surface_geofac(BcV);
               BcW->elmt->Surface_geofac(BcW);
     }      
   }
} 
  

void PBC1D::SetRC(char *name){

    /* set-up RCR boundary condition */
    R1 = new double[Nout*3];	// proximal resistance
		C1 = R1+Nout;		// capacitance
		R2 = C1+Nout;		// distal resistance
    memset(R1,'\0',Nout*3*sizeof(double));
    flowrate_RCR_old = new double[Nout];
    memset(flowrate_RCR_old,'\0',Nout*sizeof(double));
    /* check if RCRfile exist if yes read values of R1, R2 and C1,
       otherwise compute R1, R2 and C1, then create the RCRfile */

    int i;
    FILE *pRC_File;
    char fname_RC[BUFSIZ];
    sprintf (fname_RC, "%s.RCR", name );
    pRC_File = fopen(fname_RC,"r");
    if (pRC_File==NULL){
      ROOTONLY
        pRC_File = fopen(fname_RC,"w");

      for (i = 0; i < Nout; ++i){
        R2[i] = impedance_modes[i][0].real;// /sqrt(Area[i]/M_PI);
        R2[i] *= (0.1*0.1)/(ni_small_atree*density_small_atree); //scaling  L^2/mu
        C1[i]  = 0.2/R2[i];
				R1[i]  = 0.0;
        ROOTONLY
          fprintf(pRC_File,"%2.16f  %2.16f %2.16f \n",R1[i],C1[i],R2[i]);
      }
      ROOTONLY
         fclose(pRC_File);
    }
    else{
      for (i = 0; i < Nout; ++i)
        fscanf(pRC_File," %lf %lf %lf ",&R1[i],&C1[i],&R2[i]);
      fclose(pRC_File);
    }
    /* print summary */

    ROOTONLY{
      for (i = 0; i < Nout; ++i)
        fprintf(stdout,"PBC1D::setRCR -- outlet %d: R1 = %f  C1 = %f R2 = %f \n",i,R1[i],C1[i],R2[i]);
    }
}

void PBC1D::ResetRC(char *name){
    int i;
    FILE *pRC_File;
    char fname_RC[BUFSIZ];
    sprintf (fname_RC, "%s.RCR", name );
    pRC_File = fopen(fname_RC,"r");
    if (pRC_File==NULL){
      ROOTONLY
         fprintf(stdout,"PBC1D::resetRC -- can not open RCfile \n");
    }
    else{
      for (i = 0; i < Nout; ++i)
        fscanf(pRC_File," %lf %lf %lf ",&R1[i],&C1[i],&R2[i]);
      fclose(pRC_File);
      ROOTONLY{
        for (i = 0; i < Nout; ++i)
          fprintf(stdout,"PBC1D::resetRCR -- outlet = %d R1 = %f  C1 = %f R2 = %f \n",i,R1[i],C1[i],R2[i]);
      }
    }
}

void PBC1D::ReadFlowrateHistory(char *name){
 /* if history of flow rate exists - get it */
  FILE *pFile;
  char fname[BUFSIZ];

  /* open name.imp file for reading
     if file is not found return,
     flow rate history remains zero
     if file is present analyse it. */
  sprintf (fname, "%s.imp", name);
  pFile = fopen(fname,"r");
  if (pFile==NULL){
    ROOTONLY
      fprintf(stdout,
      "PBC1D::readFlowrateHistoryNew -- imp file was not found, setting flow rate history to 0.0 \n");
    return;
  }        
    
  /* analyse the data in name.imp file */
  int i,j;
  char buf[BUFSIZ];
  int Nsteps_modes; /* number of steps if space=='P' , number of modes if space=='M'*/ 
  int Noutlets_local;
  char space,Convolve_Fix;
  
  /* 1. check in if the flow rate history is saved in modal or physical space */
  fgets (buf, BUFSIZ, pFile);
  sscanf(buf,"%c",&space);

  /* 2. read the number of outlets */
  fgets (buf, BUFSIZ, pFile);
  sscanf(buf,"%d",&Noutlets_local);

  /* 3. read the number of steps if flow rate is given in physical space 
  or number of modes if it is given in modal space */
  fgets (buf, BUFSIZ, pFile);
  sscanf(buf,"%d",&Nsteps_modes);

  /* 4. read instruction on how to proceed, 
        Convolve_Fix =='F' compute pressure from given flow rate history
        and computed Impedance and fix it
        Convolve_Fix =='C' means keep updating flow rate every time step 
        and use convolution of flow rate and impedance to compute Pressure.*/
  fgets (buf, BUFSIZ, pFile);
  sscanf(buf,"%c",&Convolve_Fix);

 
  
  /* if space = modal we have enough information to proceed */
  
  if (Noutlets_local != Nout){ // inconsistent file
     ROOTONLY
       fprintf(stdout,
       "PBC1D::readFlowrateHistoryNew -- Number of outlets (%d) does not match number of outlet in *imp file (%d) \n",
       Nout,Noutlets_local);
     fclose(pFile);  
     return;  
  }
  
  char *p;
  rewind(pFile);	
  while (p = fgets (buf, BUFSIZ, pFile)){
    if (strstr (p, "FlowRateData")){
        break;
    }
  }
  if (space=='M'){	
    for (i = 0; i < Nout; i++){
    	fscanf(pFile,"%lf ",&A_cos_sin[i][0]);
    	for (j = 1; j <= Nsteps_modes; j++)
          fscanf(pFile,"%lf %lf",&A_cos_sin[i][2*j-1],&A_cos_sin[i][2*j]);
    }	
  }
  else{
  
     if (Nsteps_modes == Nsteps_per_cycle){
       for (i = 0; i < Nout; i++){
       	 for (j = 0; j < Nsteps_per_cycle; j++)	
    	    fscanf(pFile,"%lf", &FR_history[i][j]);
       }
     } 
     else{ // need to interpolate
       double dt = dparam("DT");  
       double dt_local = Tperiod/(Nsteps_modes-1.0);
       double tj, inv_dt_local_p2,coef_a,coef_b,coef_c;
       int J0,J1,J2;
       inv_dt_local_p2 = 1.0/(dt_local*dt_local);
       double *FR_local;
       FR_local = new double[Nsteps_modes];
       
       for (i = 0; i < Nout; i++){
       	
       	 for (j = 0; j < Nsteps_modes; j++)	
    	    fscanf(pFile,"%lf", &FR_local[j]);
       
       	 for (j = 0; j < Nsteps_per_cycle; j++){	
            tj = j*dt;
#if defined (__blrts__)            
            J1 = (int) round(tj/dt_local);
#else
            J1 = (int) (tj/dt_local);
#endif 
            if ( (J1 > 0) && J1 < (Nsteps_modes-1) ){
               J0 = J1-1;
               J2 = J1+1;
            }
            else{
                if (J1 <= 0){
                  J0 = 0;
                  J1 = 1;
                  J2 = 2;   
                }
                else{ 
                  J2 = Nsteps_modes-1;
                  J1 = J2-1;
                  J0 = J2-2;
                }
            }
            /*  compute interpolation coefficients use second order Lagrange */
            coef_a = (tj-dt_local*J1)*(tj-dt_local*J2)*(0.5*inv_dt_local_p2);
            coef_b = (tj-dt_local*J0)*(tj-dt_local*J2)*(-inv_dt_local_p2);
            coef_c = (tj-dt_local*J1)*(tj-dt_local*J0)*(0.5*inv_dt_local_p2);  
            FR_history[i][j] = FR_local[J0]*coef_a+FR_local[J1]*coef_b+FR_local[J2]*coef_c;
         } // end of for (j = 0; ...
       }    // end of for (i = 0; ...
       delete[] FR_local;
     }
   }

   rewind(pFile);
   j = 0;
   while (p = fgets (buf, BUFSIZ, pFile)){
     if (strstr (p, "FlowRateModal")){
       fscanf(pFile,"%d", &j); // total number of fourier modes per outlet        
       for (i = 0; i < Nout; i++){
         for (j = 0; j < (Nimpedance_modes*2+1); j++)
           fscanf(pFile,"%lf", &A_cos_sin[i][j]);
       }
       j = 1; //flag telling that data in modal space is also provided
       break;
     }
   }

   fclose(pFile);

     /* if data in modal space is provided - skip the next section */
    double *FR_history_temp;
    if (j == 1) 
      goto skip_Ftransform;
    
    /* if time > Tcycle do fourier transform of flow-rate history, obtain fourier coef.  */
    /* if Flow-rate is given in Physical space and  time > Tcycle 
        obtain fourier coefficient */
    /* if Flow-rate is given in modal space transform it into physical space.  */
    FR_history_temp = dvector(0,Nsteps_per_cycle-2);    

    if ((dparam("t") > Tperiod) && (space=='P')){
      for (i = 0; i < Nout; i++){
        memcpy(FR_history_temp,FR_history[i],(Nsteps_per_cycle-1)*sizeof(double));
        FilterDFT(Nsteps_per_cycle-1, Nimpedance_modes, FR_history_temp,A_cos_sin[i],0);
      }
    }
    free(FR_history_temp);

    skip_Ftransform:

    if ( space=='M'){
      double dt = dparam("DT");
      double arg;
      int k;
      
      for (i = 0; i < Nout; i++){
        for (j = 0; j < Nsteps_per_cycle; j++){
          arg = omega_small_atree*j*dt;
          FR_history[i][j] = A_cos_sin[i][0];
          for (k=1; k <= Nsteps_modes; k++)
          FR_history[i][j] += A_cos_sin[i][k*2-1]*cos(arg*k) + A_cos_sin[i][k*2]*sin(arg*k);
        }
      }
    }
}

void PBC1D::SaveFlowrateHistory(char *name){

  int i,j;
  FILE *pFile;
  char fname[BUFSIZ];

  sprintf (fname, "%s.imp", name);
  pFile = fopen(fname,"w");

  fprintf(pFile," P /* P - physical, M - modal space */ \n");
  fprintf(pFile,"%d /* number of outlets */ \n",Nout);
  fprintf(pFile,"%d /* Number of steps   */  \n",Nsteps_per_cycle);  
  fprintf(pFile,"C  /* instruction: C - P(t) = conv. (F,Z) , F - fix P(t) */ \n");
  fprintf(pFile,"FlowRateData \n");
  for (i = 0; i < Nout; i++){
    for (j = 0; j < Nsteps_per_cycle; j++)
      fprintf(pFile," %.10f \n",FR_history[i][j]);
  }
  fprintf(pFile,"FlowRateModal \n");
  fprintf(pFile,"%d  /* number of fourier modes */ \n",Nimpedance_modes*2+1);
  for (i = 0; i < Nout; i++){
    for (j = 0; j < (Nimpedance_modes*2+1); j++)
      fprintf(pFile," %.10f \n",A_cos_sin[i][j]);
  }
  fclose(pFile);
}


double PBC1D::GetPval(Bndry *Pbc){

   int i;
   for (i = 0; i < Nout; ++i){
     if (strcmp(Pbc->blabel,standard_labels[i]) == 0)
        return Pressure[i];
   }  
   /* if no match for standard_labels was found return zero */
   return 0.0;
}


void PBC1D::UpdateTimestepCycle(){
   double t = dparam("t");
   double dt = dparam("DELT");
   /* parameter "t" = time at the end of the cuurent time step,
      that is t^(n+1), thus we need to substruct "dt" from "t" to get 
      the time_step_in_cycle correctly */

   t -= dt; 

   t = t-((int)((t+dt*0.1)/Tperiod))*Tperiod;
     time_step_in_cycle = int ((t+0.1*dt)/dt);
}

void PBC1D::UpdateFRHistory(double *flowrate){

  /* INPUT: values of flowd-rate at outlets at current time 
     function updates arrays FR_history where flow-rate history 
     over a cycle is stored 
     if (t == Tperiod-dt) filter FR_history using Fourier transform 
     cut all frecuencies higher then  Nimpedance_modes 
     store values of Fourier coefficients.
   
     starting from the end of second cycle average values of Fourier
     coefficients with those from previous cycle
  
*/

  int i,FLAG_filter;
  for (i = 0; i < Nout; ++i)
    FR_history[i][time_step_in_cycle] = flowrate[i];


   /* predict flow-rate for the rest of the first time period */
  if   (dparam("t") < (Tperiod-dparam("DELT"))){
    double factor;
    int j;

    for (i = 0; i < Nout; ++i){
      factor = FR_history[i][time_step_in_cycle];
      for (j = time_step_in_cycle; j < Nsteps_per_cycle; j++)
        FR_history[i][j] = factor*cos(M_PI*0.5*(j-time_step_in_cycle)/(Nsteps_per_cycle-time_step_in_cycle));
    }
  }

  if (dparam("t") < (1.5*Tperiod))
     FLAG_filter = 0;
  else
     FLAG_filter = 1;
 
  if  (time_step_in_cycle == (Nsteps_per_cycle-2) ) {
    fprintf(stdout,"PBC1D::update_FR_history, FLAG_filter = %d time_step_in_cycle = %d Nsteps_per_cycle = %d \n",FLAG_filter,time_step_in_cycle,Nsteps_per_cycle);
    for (i = 0; i < Nout; ++i)
      FilterDFT(Nsteps_per_cycle-1, Nimpedance_modes, FR_history[i],A_cos_sin[i],FLAG_filter);
  }
}

void PBC1D::ComputePressureSteady(double *flowrate){

  int i;
  double Po=0.0,alpha;
  double inv_dt = 1.0/dparam("DELT");
  double DPSCAL = dparam("DPSCAL");
  double Rp, Rd, alpha1, alpha2;

	for (i = 0; i < Nout; ++i){	//Alireza: 2-element Windkessel model
		Rp = R1[i]; Rd = R2[i];
		if (!Rp){	//Alireza: 2-element Windkessel model
  		alpha = Rd*C1[i]*inv_dt;
    	Pressure[i] = DPSCAL*(Rd*flowrate[i]+Po+Pressure[i]*alpha)/(1.0+alpha);
    	flowrate_RCR_old[i] = flowrate[i];
  	}
		else{	//Alireza: 3-element Windkessel model
			alpha1 = Rp*C1[i]*inv_dt; alpha2 = Rd*C1[i]*inv_dt;
    	Pressure[i] = DPSCAL*((Rp+Rd)*flowrate[i]+Rd*alpha1*(flowrate[i]-flowrate_RCR_old[i])+Po+Pressure[i]*alpha2)/(1.0+alpha2);
    	flowrate_RCR_old[i] = flowrate[i];
		}
	}
}

void PBC1D::ComputePressureRCNonsteady(double *flowrate){

  double Po=0.0,alpha;
  double inv_dt = 1.0/dparam("DELT");
  int i;
  double DPSCAL = dparam("DPSCAL");
  static int INIT_FLAG=0;
  static double **Asin_coef, **Acos_coef;
  static double *Ao;
  int Number_rc_modes = 24;
  if (INIT_FLAG == 0){

     Ao = dvector(0,Nout-1);
     Asin_coef = dmatrix(0,Nout-1,0,Number_rc_modes-1);
     Acos_coef = dmatrix(0,Nout-1,0,Number_rc_modes-1);
     memset(Ao,'\0',Nout*sizeof(double));

     for (i = 0; i < Nout; ++i){
			 Ao[i] = 1.0;
       memset(Asin_coef[i],'\0',Number_rc_modes*sizeof(double));
       memset(Acos_coef[i],'\0',Number_rc_modes*sizeof(double));
     }
Asin_coef[1][0] = -0.000000;
Asin_coef[1][1] = 0.445290;
Asin_coef[1][2] = -0.179185;
Asin_coef[1][3] = -0.062988;
Asin_coef[1][4] = -0.150998;
Asin_coef[1][5] = 0.077661;
Asin_coef[1][6] = 0.019965;
Asin_coef[1][7] = 0.047670;
Asin_coef[1][8] = -0.005210;
Asin_coef[1][9] = -0.012413;
Asin_coef[1][10] = -0.000894;
Asin_coef[1][11] = -0.000005;
Asin_coef[1][12] = 0.002729;
Asin_coef[1][13] = -0.003319;
Asin_coef[1][14] = -0.000996;
Asin_coef[1][15] = 0.001082;
Asin_coef[1][16] = 0.002930;
Asin_coef[1][17] = 0.001716;
Asin_coef[1][18] = 0.001216;
Asin_coef[1][19] = 0.000138;
Asin_coef[1][20] = 0.000119;
Asin_coef[1][21] = 0.000714;
Asin_coef[1][22] = 0.000304;
Asin_coef[1][23] = -0.000105;
Acos_coef[1][0] = 1.189586;
Acos_coef[1][1] = -0.227300;
Acos_coef[1][2] = -0.291123;
Acos_coef[1][3] = -0.005615;
Acos_coef[1][4] = 0.072508;
Acos_coef[1][5] = 0.100112;
Acos_coef[1][6] = -0.003292;
Acos_coef[1][7] = -0.000505;
Acos_coef[1][8] = -0.041716;
Acos_coef[1][9] = 0.006509;
Acos_coef[1][10] = -0.004455;
Acos_coef[1][11] = 0.002528;
Acos_coef[1][12] = -0.001558;
Acos_coef[1][13] = 0.001292;
Acos_coef[1][14] = 0.004483;
Acos_coef[1][15] = 0.003780;
Acos_coef[1][16] = 0.002304;
Acos_coef[1][17] = -0.001212;
Acos_coef[1][18] = -0.000923;
Acos_coef[1][19] = -0.000259;
Acos_coef[1][20] = -0.000300;
Acos_coef[1][21] = -0.000402;
Acos_coef[1][22] = -0.000234;
Acos_coef[1][23] = -0.000296;

Asin_coef[2][0] = -0.000000;
Asin_coef[2][1] = 2.484296;
Asin_coef[2][2] = -0.511727;
Asin_coef[2][3] = -0.458104;
Asin_coef[2][4] = -0.571988;
Asin_coef[2][5] = 0.136195;
Asin_coef[2][6] = 0.020224;
Asin_coef[2][7] = 0.103077;
Asin_coef[2][8] = 0.010429;
Asin_coef[2][9] = 0.015921;
Asin_coef[2][10] = 0.009467;
Asin_coef[2][11] = 0.016223;
Asin_coef[2][12] = -0.019136;
Asin_coef[2][13] = 0.003484;
Asin_coef[2][14] = 0.013837;
Asin_coef[2][15] = 0.001671;
Asin_coef[2][16] = -0.005214;
Asin_coef[2][17] = -0.006719;
Asin_coef[2][18] = 0.005969;
Asin_coef[2][19] = 0.004835;
Asin_coef[2][20] = 0.002678;
Asin_coef[2][21] = 0.003294;
Asin_coef[2][22] = -0.001092;
Asin_coef[2][23] = -0.001204;
Acos_coef[2][0] = 4.604210;
Acos_coef[2][1] = -0.804558;
Acos_coef[2][2] = -1.565847;
Acos_coef[2][3] = -0.143919;
Acos_coef[2][4] = 0.216213;
Acos_coef[2][5] = 0.317577;
Acos_coef[2][6] = 0.003618;
Acos_coef[2][7] = 0.034981;
Acos_coef[2][8] = -0.068688;
Acos_coef[2][9] = 0.021200;
Acos_coef[2][10] = -0.009799;
Acos_coef[2][11] = -0.011636;
Acos_coef[2][12] = -0.006179;
Acos_coef[2][13] = 0.020968;
Acos_coef[2][14] = 0.002307;
Acos_coef[2][15] = -0.003286;
Acos_coef[2][16] = -0.004565;
Acos_coef[2][17] = 0.008501;
Acos_coef[2][18] = 0.009874;
Acos_coef[2][19] = -0.000037;
Acos_coef[2][20] = -0.003922;
Acos_coef[2][21] = -0.004068;
Acos_coef[2][22] = 0.000377;
Acos_coef[2][23] = -0.000229;

Asin_coef[3][0] = -0.000000;
Asin_coef[3][1] = 0.778630;
Asin_coef[3][2] = -0.033367;
Asin_coef[3][3] = -0.111566;
Asin_coef[3][4] = -0.181885;
Asin_coef[3][5] = -0.042586;
Asin_coef[3][6] = -0.036422;
Asin_coef[3][7] = 0.049674;
Asin_coef[3][8] = 0.019296;
Asin_coef[3][9] = 0.025934;
Asin_coef[3][10] = 0.011611;
Asin_coef[3][11] = 0.000424;
Asin_coef[3][12] = -0.016480;
Asin_coef[3][13] = -0.006750;
Asin_coef[3][14] = -0.004467;
Asin_coef[3][15] = -0.003787;
Asin_coef[3][16] = 0.006445;
Asin_coef[3][17] = 0.000483;
Asin_coef[3][18] = 0.000772;
Asin_coef[3][19] = 0.000062;
Asin_coef[3][20] = -0.000353;
Asin_coef[3][21] = -0.000619;
Asin_coef[3][22] = -0.000765;
Asin_coef[3][23] = -0.000827;
Acos_coef[3][0] = 1.422126;
Acos_coef[3][1] = -0.112834;
Acos_coef[3][2] = -0.481964;
Acos_coef[3][3] = -0.118227;
Acos_coef[3][4] = -0.005797;
Acos_coef[3][5] = 0.070345;
Acos_coef[3][6] = 0.065351;
Acos_coef[3][7] = 0.044728;
Acos_coef[3][8] = 0.000544;
Acos_coef[3][9] = 0.005489;
Acos_coef[3][10] = -0.024514;
Acos_coef[3][11] = -0.021025;
Acos_coef[3][12] = -0.009657;
Acos_coef[3][13] = 0.013400;
Acos_coef[3][14] = 0.000489;
Acos_coef[3][15] = 0.011419;
Acos_coef[3][16] = 0.004271;
Acos_coef[3][17] = -0.003366;
Acos_coef[3][18] = -0.001143;
Acos_coef[3][19] = -0.000579;
Acos_coef[3][20] = -0.001219;
Acos_coef[3][21] = -0.000547;
Acos_coef[3][22] = 0.000781;
Acos_coef[3][23] = -0.000128;

Asin_coef[4][0] = -0.000000;
Asin_coef[4][1] = 1.516236;
Asin_coef[4][2] = -0.171087;
Asin_coef[4][3] = -0.262555;
Asin_coef[4][4] = -0.327033;
Asin_coef[4][5] = -0.024935;
Asin_coef[4][6] = 0.009933;
Asin_coef[4][7] = 0.098320;
Asin_coef[4][8] = 0.048728;
Asin_coef[4][9] = 0.021884;
Asin_coef[4][10] = -0.000834;
Asin_coef[4][11] = -0.020888;
Asin_coef[4][12] = -0.009980;
Asin_coef[4][13] = -0.000298;
Asin_coef[4][14] = 0.000478;
Asin_coef[4][15] = 0.005916;
Asin_coef[4][16] = 0.002569;
Asin_coef[4][17] = -0.001567;
Asin_coef[4][18] = 0.004371;
Asin_coef[4][19] = 0.001023;
Asin_coef[4][20] = -0.000759;
Asin_coef[4][21] = 0.001808;
Asin_coef[4][22] = 0.001090;
Asin_coef[4][23] = 0.000033;
Acos_coef[4][0] = 2.518640;
Acos_coef[4][1] = -0.328228;
Acos_coef[4][2] = -0.930788;
Acos_coef[4][3] = -0.135116;
Acos_coef[4][4] = 0.041390;
Acos_coef[4][5] = 0.204423;
Acos_coef[4][6] = 0.089588;
Acos_coef[4][7] = 0.046089;
Acos_coef[4][8] = -0.034630;
Acos_coef[4][9] = -0.012215;
Acos_coef[4][10] = -0.033807;
Acos_coef[4][11] = -0.011331;
Acos_coef[4][12] = 0.009254;
Acos_coef[4][13] = 0.010264;
Acos_coef[4][14] = 0.007906;
Acos_coef[4][15] = 0.004663;
Acos_coef[4][16] = -0.000819;
Acos_coef[4][17] = -0.001270;
Acos_coef[4][18] = 0.001076;
Acos_coef[4][19] = 0.000024;
Acos_coef[4][20] = -0.000364;
Acos_coef[4][21] = 0.000162;
Acos_coef[4][22] = 0.000367;
Acos_coef[4][23] = -0.000083;

     INIT_FLAG = 1;
  }

     /* fix R1 compute R2,  R2/R1 ratio is given in fourier space, then R2 = R1 * ratio */
		 /* RC model only; distal resistance considered */
  for (i = 1; i < Nout; ++i)
    R2[i] = R2[0]*Acos_coef[i][0]*Ao[i];

  double t = dparam("t"),sin_alpha,cos_alpha;
  double omega_t = 2.0*M_PI*t/Tperiod;
  int j;
  
  for (i = 1; i < Number_rc_modes; ++i){
    alpha = omega_t*(double) i;
    sin_alpha = sin(alpha);
    cos_alpha = cos(alpha);
    for (j = 1; j < Nout; ++j)
      R2[j] += R2[0]*(Asin_coef[j][i]*sin_alpha+Acos_coef[j][i]*cos_alpha);
  }

  for (i = 0; i < Nout; ++i){
     alpha = R2[i]*C1[i]*inv_dt;
     Pressure[i] = DPSCAL*(R2[i]*flowrate[i]+Po+Pressure[i]*alpha)/(1.0+alpha);
     flowrate_RCR_old[i] = flowrate[i];
  }
}

void PBC1D::ComputePressure(){
/* compute pressure using convolution 
   presure is computed in [cm,sec,gr], so scale it to nondim units */

   int i,time_step_in_cycle_advanced;
   double DPSCAL = dparam("DPSCAL");

   if (time_step_in_cycle == (Nsteps_per_cycle-2)) 
       time_step_in_cycle_advanced = 0;
   else
       time_step_in_cycle_advanced = time_step_in_cycle+1;


   for (i = 0; i < Nout; ++i){
      if (!parallel_convolution)
         ParallelConvolve(Nsteps_per_cycle-1, FR_history[i], impedance[i], &Pressure[i], time_step_in_cycle_advanced);
      else{
       ParallelConvolve(Nsteps_per_cycle-1, FR_history[i], impedance[i], &Pressure[i], Imp_index_start, (Imp_index_start+Imp_length_local), time_step_in_cycle_advanced);

       DO_PARALLEL{
         double temp = 0.0;
         gdsum(&Pressure[i],1,&temp);
       }
     }
     Pressure[i] *= DPSCAL*0.1/(Tperiod*density_small_atree); 
  }
}

void PBC1D::ComputeCoronaryPressure(double *flowrate){

    int i, j;
    static int INITCOR_FLAG = 0;
    static double *r1, *r2, *r3, *c1, *c2;
    static double *beta1, *beta2, *beta3, *phi0, *phi1, *Pc1, *Qc2, 
									*Pc2, *Pc1_old, *Pc2_old, *Pim, *Pim_old;
    double DPSCAL = dparam("DPSCAL");
  
    double Po=0.0,alpha;
    double inv_dt = 1.0/dparam("DELT");
    double Rp, Rd, alpha1, alpha2;

    if (INITCOR_FLAG == 0){
			beta1 = dvector(0, NoutCor-1);
			beta2 = dvector(0, NoutCor-1);
			beta3 = dvector(0, NoutCor-1);
			phi0 = dvector(0, NoutCor-1);
			phi1 = dvector(0, NoutCor-1);

			r1 = dvector(0, NoutCor-1);
			r2 = dvector(0, NoutCor-1);
			r3 = dvector(0, NoutCor-1);
			c1 = dvector(0, NoutCor-1);
			c2 = dvector(0, NoutCor-1);
	
			Pc1 = dvector(0, NoutCor-1);
			Pc2 = dvector(0, NoutCor-1);
			Pc1_old = dvector(0, NoutCor-1);
			Pc2_old = dvector(0, NoutCor-1);
			Qc2 = dvector(0, NoutCor-1);
			Pim = dvector(0, NoutCor-1);
			Pim_old = dvector(0, NoutCor-1);
	
			for (i = 0; i < NoutCor; ++i){
				r1[i] = Ra[i];
				r2[i] = Ra_mic[i];
				r3[i] = Rv[i];
				c1[i] = Ca[i];
				c2[i] = Cim[i];
				beta1[i] = r2[i]/r3[i];
				beta2[i] = r2[i]*c2[i]*inv_dt;
				beta3[i] = c1[i]*inv_dt;
				phi0[i] = r2[i]*(1.0+beta1[i]+beta2[i]);
				phi1[i] = beta3[i]+1.0/r2[i]-1.0/phi0[i];
				Pc1[i] = 0.0;
				Pc2[i] = 0.0;
				Pc1_old[i] = 0.0;
				Pc2_old[i] = 0.0;
				Qc2[i] = 0.0;
				Pim[i] = 0.0;
				Pim_old[i] = 0.0;
			}
			INITCOR_FLAG = 1;
    }

 	for (i = 0; i < NoutRCR; ++i){
    Rp = R1[i]; Rd = R2[i];
    if (!Rp){ //Alireza: 2-element Windkessel model
      alpha = Rd*C1[i]*inv_dt;
      Pressure[i] = DPSCAL*(Rd*flowrate[i]+Po+Pressure[i]*alpha)/(1.0+alpha);
      flowrate_RCR_old[i] = flowrate[i];
    }
    else{ //Alireza: 3-element Windkessel model
      alpha1 = Rp*C1[i]*inv_dt; alpha2 = Rd*C1[i]*inv_dt;
      Pressure[i] = DPSCAL*((Rp+Rd)*flowrate[i]+Rd*alpha1*(flowrate[i]-flowrate_RCR_old[i])+Po+Pressure[i]*alpha2)/(1.0+alpha2);
      flowrate_RCR_old[i] = flowrate[i];
    }
	}

  double t = dparam("t"), sin_alpha, cos_alpha;
  double omega_t = Pim_wnum*t;

  for (j = 0; j < NoutCor; ++j)
  	Pim[j] = Pim_scal[j]*Pim_ccos[0];
  for (i = 1; i < Pim_nmodes; ++i){
    alpha = omega_t*(double) i;
    sin_alpha = sin(alpha);
    cos_alpha = cos(alpha);
    for (j = 0; j < NoutCor; ++j)
      Pim[j] += Pim_scal[j]*(Pim_csin[i]*sin_alpha + Pim_ccos[i]*cos_alpha);
  }
	for (i = 0; i < NoutCor; ++i){
		int ii = i + NoutRCR;
		Pressure[ii] = DPSCAL*(r1[i]*flowrate[ii] + 
									(1.0/phi1[i])*(flowrate[ii]+beta3[i]*Pc1_old[i]+(1.0/phi0[i])*
									(beta1[i]*Po+beta2[i]*(Pim[i]-Pim_old[i]+Pc2_old[i]))));
		Pc1[i] = Pressure[ii] - flowrate[ii]*r1[i];
		Qc2[i] = flowrate[ii] - c1[i]*(Pc1[i]-Pc1_old[i])*inv_dt;
		Pc2[i] = Pc1[i] - Qc2[i]*r2[i];
		Pc1_old[i] = Pc1[i];
		Pc2_old[i] = Pc2[i];
		Pim_old[i] = Pim[i];
		flowrate_RCR_old[ii] = flowrate[ii];
//		if (i==0) ROOTONLY fprintf(stdout,"%e %e %e %e %e %e %e\n",r1[i],c1[i],r2[i],c2[i],r3[i],Pim_scal[i],Pim[i]);
	}
}
#endif
