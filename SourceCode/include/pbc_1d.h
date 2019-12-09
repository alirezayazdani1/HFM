#ifndef PBC_1D_H
#define PBC_1D_H

#define Lrr_small_atree  50.0
#define Rmin_small_atree (0.02)
#define k1_small_atree   (2.0E+7) 
#define k2_small_atree   (-22.53)
#define k3_small_atree   (8.65E+5)
#define ni_small_atree   0.03525
#define omega_small_atree (2.0*M_PI/0.699148936170213)
#define density_small_atree 1.0

class my_dcmplx{

 public:
    double real,imag;

    void set_val(double re, double im);
    void my_dcmplx_axpy(double a, my_dcmplx x, my_dcmplx y);
    void my_dcmplx_axpy(double a, my_dcmplx x, double b, double y);
    void my_dcmplx_dot(my_dcmplx x, my_dcmplx y);
    void my_dcmplx_inv_dot(my_dcmplx x, my_dcmplx y);
    void my_dcmplx_sin(my_dcmplx x);
    void my_dcmplx_cos(my_dcmplx x);
    void my_dcmplx_sqrt(my_dcmplx x);
    void my_dcmplx_conj(my_dcmplx x);

};


class PBC1D {

  public:
  
  int Nout;               /* number of outlets                                    */ 
  int Ninl;               /* number of inlets                                     */
  int Nimpedance_modes;   /* number of impedance modes                            */
  int Nsteps_per_cycle;   /* number of time steps per cycle                       */
  int time_step_in_cycle; /* index of time step in currend cycle                  */ 
  int Imp_index_start;    /* index of element of global array corresponding 
                             to the first element in local partition              */ 
  int Imp_length_local;   /* number of elements of global array to be stored
                             in local partition                                   */
  int parallel_convolution; /* = 1 for parallel convolution, = 0 otherwise         */
 
  double Tperiod;         /* length of cycle - time period                        */
  double **FR_history;    /* Flow-Rate history; dimension[Nout][Nsteps_per_cycle] */
  double **A_cos_sin;     /* Fourier coef. of flow-rate history                   */
  double **impedance;     /* impedance if phys. space;
                                               dimension[Nout][Nsteps_per_cycle]  */
  double *Pressure;       /* values of pressure computed by convolution                                
                                                   MPIdimension[Nout]                */ 
  double *Area;           /* area of all outlet and inlet boundaries              */  

  char *standard_labels[26]; /* standatd labels for inlets and outlets            */

  int *Nfaces_per_outlet;    /* number of faces at each outlet - local variable   */
  int **ID_faces_per_outlet; /* local indices of faces at each outlet             */ 

  int *Nfaces_per_inlet;    /* number of faces at each outlet - local variable   */
  int **ID_faces_per_inlet; /* local indices of faces at each outlet             */

  my_dcmplx **impedance_modes; /* values of impedance in fourier space
                                  for each outlet                                 */

  double *R1,*C1,*R2,*flowrate_RCR_old;

  int *Nnodes_inlet, *Nnodes_outlet;

	int NoutRCR, NoutCor, Pim_nmodes;
	double Pim_wnum;
	double *Ra, *Ca, *Ra_mic, *Cim, *Rv, *Pim, *Pim_scal, *Pim_ccos, *Pim_csin;

  int  Create(Bndry *Ubc);
	void ReadLPN(FILE* fp, Domain* omega); 
  void SetGeofac(Bndry *Ubc, Bndry *Vbc, Bndry *Wbc);
  void SetRC(char *name);
  void ResetRC(char *name);
  void ReadFlowrateHistory(char *name);
  void SaveFlowrateHistory(char *name);
  double GetPval(Bndry *Ubc);
  void UpdateTimestepCycle();
  void UpdateFRHistory(double *flowrate); 
  void ComputePressure(); 
  void ComputePressureSteady(double *flowrate);
  void ComputePressureRCNonsteady(double *flowrate);
	void ComputeCoronaryPressure(double *flowrate);
};


class SMALL_VESSEL{

 public:
  double  radius;    /* radius of vessel        */
  double  Length;    /* length of vessel        */
  my_dcmplx* Zo;       /* Impedance at the inlet  */
  my_dcmplx* Zl;       /* Impedance at the outlet */
  char type;         /* type = T if this is a terminal vessel
                   = C if this is connecting vessel */
  double Compliance; /* Complines */
  double Ws;         /* Womersley number */
  my_dcmplx* Fj;

  SMALL_VESSEL* child_1;
  SMALL_VESSEL* child_2;

  SMALL_VESSEL(double ro);

  my_dcmplx GetZL(int mode);
  my_dcmplx GetZ0(int mode);

};

#endif
