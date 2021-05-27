#ifndef H_NEKTAR
#define H_NEKTAR

/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $Source: /homedir/cvs/Nektar/include/nektar.h,v $  
 * $Revision: 1.9 $
 * $Date: 2006/08/15 19:31:06 $    
 * $Author: ssherw $  
 * $State: Exp $   
 *---------------------------------------------------------------------------*/

#include "hotel.h"
#include "stokes_solve.h"

#ifdef PBC_1D
#include "pbc_1d.h"
#endif

using namespace MatSolve;

#ifdef PARALLEL
#include <mpi.h>
#define exit(a) {MPI_Abort(MPI_COMM_WORLD, -1); exit(a-1);}
extern "C"
{
#include "gs.h"
}
#endif

#define MAXDIM 3

/* general include files */
#include <math.h>
#include "veclib.h"

/* parameters */
#define HP_MAX  128   /* Maximum number of history points */

typedef enum {
  Splitting,                    /* Splitting scheme                        */
  StokesSlv,                    /* Stokes solver                           */
  SubStep,                      /* Substepping scheme                      */
  SubStep_StokesSlv,            /* Stokes solver                           */ 
  Adjoint_StokesSlv,
  TurbCurr                      /* Splitting scheme with turbulent current */ 
} SLVTYPE;


typedef enum {                  /* ......... ACTION Flags .......... */
  Rotational,                   /* N(U) = U x curl U                 */
  Convective,                   /* N(U) = U . grad U                 */
  Stokes,                       /* N(U) = 0  [drive force only]      */
  Alternating,                  /* N(U^n,T) = U.grad T               */
                                /* N(U^(n+1),T) = div (U T)          */
  Linearised,                   /* Linearised convective term        */
  StokesS,                      /* Steady Stokes Solve               */
  Oseen,                        /* Oseen flow                        */
  OseenForce,
  StokesForce,                  /* Stokes forcing terms              */
  LambForce,                    /* N(U) = u^2 Curl w - j             */
  Pressure,                     /* div (U'/dt)                       */
  Viscous,                      /* (U' - dt grad P) Re / dt          */
  Prep,                         /* Run the PREP phase of MakeF()     */
  Post,                         /* Run the POST phase of MakeF()     */
  Transformed,                  /* Transformed  space                */
  Physical,                     /* Phyiscal     space                */
  TransVert,                    /* Take history points from vertices */
  TransEdge,                    /* Take history points from mid edge */
  Poisson   = 0,
  Helmholtz = 1,
  Laplace   = 2
} ACTION;

typedef struct hpnt {           /* ......... HISTORY POINT ......... */
  int         id         ;      /* ID number of the element          */
  int         i, j, k    ;      /* Location in the mesh ... (i,j,k)  */
  char        flags             /* The fields to echo.               */
            [_MAX_FIELDS];      /*   (up to maxfields)               */
  ACTION      mode       ;      /* Physical or Fourier space         */
  struct hpnt *next      ;      /* Pointer to the next point         */
} HisPoint;

// Ce107
typedef struct intepts {
  int     npts;
  Coord    X;
  double **ui;
} Intepts;


typedef struct gf_ {                  /* ....... Green's Function ........ */
  int           order                 ; /* Time-order of the current field   */
  Bndry        *Gbc                   ; /* Special boundary condition array  */
  Element_List *basis                 ; /* Basis velocity (U or W)           */
  Element_List *Gv[MAXDIM][_MAX_ORDER]; /* Green's function velocities       */
  Element_List *Gp        [_MAX_ORDER]; /* Green's function pressure         */
  double        Fg        [_MAX_ORDER]; /* Green's function "force"          */
} GreensF;


/* local structure for global number */
typedef struct gmapping {
  int nvs;
  int nvg; 
  int nes;
  int neg; 
  int *vgid;
  int *egid;
  int nfs;
  int nfg; 
  int *fgid;

  int nsolve;
  int nglobal;
  int nv_psolve;
} Gmap;

/* time dependent function inputs -- ie. womersley flow */
typedef struct time_str {  /* string input      */
  char *TimeDepend;
} Time_str;

/* info for a function of the form:
   t = (ccos[i]*cos(wnum[i]*t) + csin[i]*sin(wnum[i]*i*t)) */
   
typedef struct time_fexp {  /* fourier expansion */
  int nmodes;
  double wnum;
  double *ccos;
  double *csin;
} Time_fexp;

typedef enum{
  ENCODED,                /* fourier encoded values */
  RAW,                    /* actual measure,ents    */
  PRESSURE,
  MASSFLOW
} TwomType;

typedef struct time_wexp {  /* Womersley expansion fourier expansion */
  int nmodes;
  int nraw;
  double scal0;
  double radius;
  double axispt[MAXDIM];
  double axisnm[MAXDIM];
  double wnum;
  double *ccos;
  double *csin;
  double *raw;
  TwomType type;
  TwomType form;
} Time_wexp;

typedef union tfuninfo { 
  Time_str  string [MAXDIM];   
  Time_fexp fexpand[MAXDIM];   
  Time_wexp wexpand;   
} TfunInfo;              

typedef enum{
  String,                 /* input string            */
  Fourier,                /* input fourier expansion */
  Womersley               /* input womersley expansion */
} TfunType;

typedef struct tfuct {
  TfunType    type;             /* type of time function input */
  TfunInfo    info;             /* info for function           */
}Tfunct;


#ifdef BMOTION2D
typedef enum motiontype {
  Prescribed,             // no body motion 
  Heave,                 // heave motion only
  InLine,                // Inline motion only
  Rotation,              // Rotation only
  HeaveRotation          // heave and rotation
} MotionType;

typedef struct motion2d {
  MotionType motiontype; /* enumeration list with motion type           */
  double stheta;       /* inertial moment of the section               */
  double Itheta;       /* static moment of the section                 */
  double mass;         /* mass of the section                          */
  double cc[2];        /* damping coefficients for each degree of      */
                       /* freedom, cc[0] for y and cc[1] for theta     */
  double ck[2];        /* stiffness coefficients                       */
  double cl[2];        /* scaling constants for force and moment       */
  double dp[2];        /* displacement at previous time level          */
  double dc[2];        /* displacement at current time level           */
  double vp[2];        /* first derivatives at previous time level     */
  double vc[2];        /* first derivatives at current time level      */
  double ap[2];        /* second derivatives at previous time level    */
  double ac[2];        /* second derivatives at current time level     */

  double f[3];        /* for storing forces at the previous step.      */
  double aot;          /* Angle of attach in HM mode                   */
  double x[3];
  double alpha[3];
  double Ured_y;
  double Ured_alpha;
  double zeta_y;
  double zeta_alpha;
  double mass_ratio;
  double inertial_mom;
} Motion2d;
#endif

#ifdef FCM
  class FCM_Domain;
  class FCM_Boundary_Info;
#endif //FCM

/* Solution Domain contains global information about solution domain */
typedef struct dom {
  int      step;
  char     *name;           /* Name of run                       */
  char     **soln;
  char     **sources;
  double   dt;

  FILE     *fld_file;       /* field file pointer                */
  FILE     *dat_file;       /* field file pointer                */
  FILE     *his_file;       /* history file pointer              */    
  FILE     *fce_file;       /* force file                        */
  FILE     *flo_file;       /* flow rate outlets/inlet           */
  FILE     *pre_file;       /* pressure at outlets/inlet        */
  FILE     *mom_file;       /* momentum flux at outlets/inlet    */
  HisPoint *his_list;       /* link list for history points      */
  Tfunct   *Tfun;           /* time dependent into conditions    */

  Element_List  *U, *V, *W, *P;  /* Velocity and Pressure fields      */
#ifdef ADR
  Element_List  **T;	             /* Temperature */
	double *Pr;
#ifdef CCLF
	double  *ind_func;							/* Indicator Function */
	int  		*number;								/* PLT number */
	double  *volume, **timeAct;			/* Activation Time */
#endif
#endif
	double **Vavg;
  Element_List  *Uf;             /* --------------------------------- */
  Element_List  *Vf;             /*        Multi-step storage         */
  Element_List  *Wf;             /* --------------------------------- */
  Element_List  *Pf;
#ifdef ADR
  Element_List  **Tf;	
#endif
  Element_List  *Lfx,*Lfy;       /* forcing terms to lamb vector      */
	
  MatSolve::StokesMatrix *StkSlv;

  double **u;                   /* Field storage */           
  double **v;                   /* Field storage */           
  double **w;                   /* Field storage */
#ifdef ADR
  double ***t;                   /* Field storage */	
#endif

  double **uf;                   /* Non-linear forcing storage */           
  double **vf;                   /* Non-linear forcing storage */           
  double **wf;                   /* Non-linear forcing storage */  
#ifdef ADR
  double ***tf;                   /* Non-linear forcing storage */  	
#endif
  
  double **lfx;                 /* lamb vector storage */
  double **lfy;                 
    

  double **us;                  /* multistep storage */
  double **vs;
  double **ws;
  double **ps;
#ifdef ADR
  double ***ts;
#endif
	

  double **ul;                  /* Lagrangian velocity for Op Int Spl. */
  double **vl;
  double **wl;

  double **mu, **mv, **mw, **mx, **my, **mz;
  
  Bndry    *Ubc,*Vbc,*Wbc;      /* Boundary  conditons               */
#ifdef ADR
  Bndry    **Tbc;      /* Boundary  conditons               */	
#endif
  Bndry    *Pbc;
#ifdef ALE
  Bndry    *dUdt, *dVdt;
#endif
  Bsystem  *Usys;               /* Velocity Helmholtz matrix system  */
  Bsystem  *Vsys;               /* Velocity Helmholtz matrix system  */
  Bsystem  *Wsys;               /* Velocity Helmholtz matrix system  */
#ifdef ADR
  Bsystem  **Tsys;               /* Temperature Helmholtz matrix system  */	
#endif
  Bsystem  *Pressure_sys;       /* pressure Poisson   matrix system  */

  double   **ForceFuncs;        /* storage for variable forcing      */
  char     **ForceStrings;      /* string definition of forcing      */

  //  Metric     *kinvis;

  // ALE structures
#ifdef ALE
  double **ztri;
  Element_List *Ptri;
  Element_List *Ptri_f;
  Bndry *Ptri_bc;
  Bsystem *Ptri_sys;
  Corner *corn;

  double  *velocity;
  double *position;
  Element_List **MeshX;
  Element_List **MeshV;
  Element_List *MeshVf;
  Bndry        **MeshBCs;
  Bsystem      *Mesh_sys;

  Element_List *Psegs;
  Element_List *PsegsF;
  Bsystem      *Psegs_sys;

  Element_List *Msegs;
  Element_List *MsegsUF;
  Element_List *MsegsVF;
  Bsystem      *Msegs_sys;

  int          *update_list;
#endif

#ifdef BMOTION2D
  FILE     *bdd_file, *bda_file;       /* motion of body file              */
  FILE     *bgy_file;                  /* energy data                      */
  
  Motion2d   *motion2d;                /* Body motion info                 */
#endif

#ifdef MAP
  // Ce107
  FILE     *int_file;       /* interpolation file pointer        */
  Map      *mapx,  *mapy;   /* Mapping in x and y                */
  MStatStr *mstat;          /* Moving Statistics information     */
  Intepts  *int_list;       /* link list for interpolation points*/
#endif

#ifdef PBC_1D
  PBC1D  *IMPBC;    /* Impedance boundary condition      */
#endif

	Element_List **delV;
	Element_List *gamdot;
#ifdef FCM
	Element_List *PartNum;
  FCM_Domain* FCM_Omega;
  FCM_Boundary_Info* FCM_BI;
#endif	//FCM
} Domain;

/* function in drive.c  */
void MakeF   (Domain *omega, ACTION act, SLVTYPE slvtype);
void solve(Element_List *U, Element_List *Uf,Bndry *Ubc,Bsystem *Ubsys,SolveType Stype,int step);
int Navier_Stokes(Domain *Omega, double t0, double tN);
void SVVExplicit(Element_List *U, Element_List *Uf);
double setQfactor(int ind, int limit, int low, int high);

/* functions in prepost.c */
void      parse_args (int argc,  char **argv);
void      PostProcess(Domain *omega, int, double);
Domain   *PreProcess (int argc, char **argv);
void      set_vertex_links(Element_List *UL);
void LocalNumScheme  (Element_List *E, Bsystem *Bsys, Gmap *gmap);
Gmap *GlobalNumScheme(Element_List *E, Bndry *Ebc);
void free_gmap(Gmap *gmap);
void free_Global_info(Element_List *Mesh, Bndry *Meshbc, 
       		     Gmap *gmap, int lnel);

Bsystem *gen_bsystem(Element_List *UL, Gmap *gmap);
 
/* functions in io.c */
void      ReadParams     (FILE *rea);
void      ReadPscals     (FILE *rea, Domain *);
void      ReadLogics     (FILE *rea);
Element_List  *ReadMesh       (FILE *rea,char *);
void      ReadKinvis     (Domain *);
void      ReadICs        (FILE *, Domain *);
void      ReadDF         (FILE *, Domain *);
void      ReadDF         (FILE *fp, int nforces, ...);
void      summary        (void);
void      ReadSetLink    (FILE *fp, Element_List *U);
void      ReadSetLink    (FILE *fp, Element *U);
Bndry    *ReadBCs        (FILE *fp, Element *U);
Bndry *ReadMeshBCs (FILE *fp, Element_List *Mesh);
Bndry *ReadMeshBCsTemp (FILE *fp, Element_List *Mesh);
Bndry    *bsort          (Bndry *, int );
void      read_connect   (FILE *name, Element_List *);
void      ReadOrderFile  (char *name,Element_List *E);
void      ReadHisData    (FILE *fp, Domain *omega);
void      ReadSoln       (FILE* fp, Domain* omega);
void      ReadDFunc      (FILE *fp, Domain *Omega);
void      ReadWave       (FILE *fp, double **wave, Element_List *U);
void      ReadTimeDepend (FILE *fp, Domain *omega);

/* structure specific to bwoptim and recurSC */
typedef struct facet{
  int  id;
  struct facet *next;
} Facet;

typedef struct fctcon{
  int ncon;
  Facet *f;
} Fctcon;


/* function in bwoptim.c */
void bandwidthopt (Element *E, Bsystem *Bsys, char trip);
void MinOrdering   (int nsols,  Fctcon *ptcon, int *newmap);
void addfct(Fctcon *con, int *pts, int n);
void free_Fctcon   (Fctcon *con, int n);

/* functions in recurrSC.c */
void Recursive_SC_decom(Element *E, Bsystem *B);

/* functions in convective.c */
void VdgradV (Domain *omega);
void CdgradV (Domain *omega);

/* functions in DN.C */
void DN(Domain *omega);

/* functions in rotational.c */
void VxOmega (Domain *omega);

/* functions in divergenceVv.c */
void DivVv (Domain *omega);

/* functions in lambforce.C */
void LambForcing(Domain *omega);

/* functions in pressure.c */
Bndry *BuildPBCs (Element_List *P, Bndry *temp);
void   SetPBCs   (Domain *omega);
#ifdef ADR
void   SetTBCs   (Domain *omega);
#endif
void Set_Global_Pressure(Element_List *Mesh, Bndry *Meshbcs);
void Replace_Numbering(Element_List *UL, Element_List *GUL);
void Replace_Local_Numbering(Element_List *UL, Element_List *GUL);
void set_delta(int Je);

/* function in stokes.c      */
void StokesBC (Domain *omega);

/* functions in analyser.c   */
void Analyser (Domain *omega, int step, double time);
//void Statistics (Domain *omega, int step, double time);

/* functions in forces */
void  forces (Domain *omega, int step, double time);
void  Forces (Domain *omega, double time, double *F,int writoutput);

/* functions in sections */
int cnt_srf_elements(Domain *omega);
Element_List *setup_surflist(Domain *omega, Bndry *Ubc, char type);


// Functions in ALE
Bndry *BuildMeshBCs(Element_List *M, Bndry *Ubc);
void Update_Mesh(Domain *Omega, int Je, double dt);
//void   Update_Mesh (Domain *Omega);
void   setup_ALE   (Domain *Omega);
void   setup_ALE   (Domain *Omega, Element_List *, Bndry *);
void Update_Mesh_Velocity(Domain *Omega);
void Set_Mesh(Domain *Omega);
void set_ALE_ICs(Domain *omega);
void set_soliton_ICs(Element_List *U, Element_List *V);
void set_soliton_BCs(Bndry *Ubc, Bndry *Vbc, char ch);
void update_paddle(Element_List *, Bndry *);

// Functions in smoother
Element_List *setup_seglist(Bndry *Ubc, char type);
Bsystem      *setup_segbsystem(Element_List *seg_list);
void          fill_seglist(Element_List *seg_list, Bndry *Ubc);
void          fill_bcs(Element_List *seg_list, Bndry *Ubc);
void          update_seg_vertices(Element_List *, Bndry *);
void update_surface(Element_List *EL, Element_List *EL_F, Bsystem *Bsys, Bndry *BCs);
//void smooth_surface(Element_List *EL, Element_List *EL_F,
//		    Bsystem *Bsys, Bndry *BCs);
void smooth_surface(Domain *, Element_List *segs, Bsystem *Bsys, Bndry *BCs);
void test_surface(Bndry *Ubc, char type);

// functions in magnetic
void Set_Global_Magnetic(Element_List *Mesh, Bndry *Meshbcs);
void ReadAppliedMag(FILE* fp, Domain* omega);

double cfl_checker     (Domain *omega, double dt);
double full_cfl_checker(Domain *omega, double dt, int *eid_max);

// functions in mpinektar.C
void init_comm(int*, char ***);
void exit_comm();
void exchange_sides(int Nfields, Element_List **Us);
void SendRecvRep(void *buf, int len, int proc);

/* functions in mlevel.C */
void Mlevel_SC_decom(Element_List *E, Bsystem *B);

/* functions in womersley.C */
void WomInit(Domain *Omega);
void SetWomSol (Domain *Omega, double time, int save);
void SaveWomSol(Domain *Omega);

void zbesj_(double *ZR, double *ZI, double *FNU, int *KODE, int *N,
		     double *CYR, double *CYI, int *NZ, int *IERR);
void WomError    (Domain *omega, double time);
void SetWomField (Domain *omega, double *u, double *v, double *w, double time);

/* functions in wannier.C */
#ifdef HEIMINEZ
void Heim_define_ICs (Element_List  *V);
void Heim_reset_ICs  (Element_List **V);
void Heim_error      (Element_List **V);
#endif
/* functions in structsolve.C */
void CalcMeshVels(Domain *omega);


#ifdef MAP
typedef struct mppng {            /* ......... Mapping ............... */
  int       NZ                  ; /* Number of z-planes                */
  double    time                ; /* Current time                      */
  double   *d                   ; /* Displacement                      */
  double   *z                   ; /*   (z-deriv)                       */
  double   *zz                  ; /*   (zz-deriv)                      */
  double   *t                   ; /* Velocity                          */
  double   *tt                  ; /* Acceleration                      */
  double   *tz                  ; /*   (tz-deriv)                      */
  double   *tzz                 ; /*   (tzz-deriv)                     */
  double   *f                   ; /* Force                             */
} Map;

typedef struct mstatstr {         /* Moving Statistics information     */
  double *x;                      /* vector holding the x-coordinate   */
  double *y;                      /* vector holding the y-coordinate   */
  int    *sgnx;                   /* vector holding the sign of dx/dt  */
  int    *sgny;                   /* vector holding the sign of dy/dt  */
  int    *nstat;                  /* vector holding the # of samples   */
  int     n;                      /* number of sampling (x,y) points   */
} MStatStr;
#endif

/* ce107 changes begin */
int       backup         (char *path1);
void      ReadIntData    (FILE *fp, Domain *omega);
void      ReadSoln       (FILE *fp, Domain *omega);
void      ReadMStatPoints(FILE* fp, Domain* omega);
void      averagefields  (Domain *omega, int nsteps, double time);
double    VolInt         (Element_List *U, double shift);
double    VolInt         (Element_List *U);

double            L2norm (Element_List *V);
void       average_u2_avg(Domain *omega, int step, double time);
void             init_avg(Domain *omega);

/* functions in interp */
void interp(Domain *omega);

/* functions in nektar.C */
void AssembleLocal(Element_List *U, Bsystem *B);

/* functions in dgalerkin.c */
void set_elmt_sides    (Element_List *E);
void set_sides_geofac  (Element_List *EL);
void Jtransbwd_Orth    (Element_List *EL, Element_List *ELf);
void EJtransbwd_Orth   (Element *U, double *in, double *out);
void InnerProduct_Orth (Element_List *EL, Element_List *ELf);
int  *Tri_nmap         (int l, int con);
int  *Quad_nmap        (int l, int con);

/* subcycle.C */
void SubCycle(Domain *Omega, int Torder);
void Upwind_edge_advection(Element_List *U, Element_List *V, 
				  Element_List *UF, Element_List *VF);
void Add_surface_contrib(Element_List *Us, Element_List *Uf);
void Fill_edges(Element_List *U, Element_List *Uf, Bndry *Ubc,double t);


#ifdef BMOTION2D
/* functions in bodymotion.c */
void set_motion (Domain *omega, char *bddfile, char *bdafile);
void Bm2daddfce (Domain *omega);
void ResetBc    (Domain *omega);
void IntegrateStruct(Domain *omega, double time, int step, double dt);
#endif

#ifdef ALE
Corner *Find_corners(Element_List *EL,Domain *omega);
void sort_gids(Element_List *tri_list, Bndry *tri_bc);
void make_diff_coeff(Bsystem *Bsys, Element_List *P, Corner *corn);
void SetHisData(Domain *omega);
void setup_2d(Domain *omega, Corner *corn, Bndry *Meshbc);
Element_List *setup_trilist(Bndry *Ubc, char type);
Bndry *setup_tribndry(Element_List *tri_list, Corner *corn);
Bsystem *setup_tribsystem(Element_List *tri_list, Bndry *tri_bc);
void      ReadPosition(Domain *omega);
Gmap *SurfaceNumScheme(Element_List *UL, Bndry *Ebc);
Bndry *bsort_zeta(Bndry *bndry_list, int nbcs);
#endif
// MSB: Added for PSE------------------------------------
// Needs to be included after Domain is defined        //
#include "stokes_solve_F.h"                              //

#ifdef ADR
#if defined (CCHF) || (CCLF)
void AddCCSource(Domain *, int);
void AddBrinkman(Domain *);
void CalcTbc(Domain *, Bndry *, int, int);
enum model {
	Fogelson,
	Anand
};

const double tscale = 0.01*10.0;				// problem time scale (sec)
const double lscale = 1.0e-6;			// problem length scale (m)
const double cscale = 1.0e-9;			// problem length scale (m)

const double mu = 3.0e-1/tscale, stddev = sqrt(0.2e-1)/tscale;		// agonist release Gaussian parameters (50 msec release time constant of ADP)

#ifdef CCHF
// Fogelson's Model
/*
enum SPECIES {
	II, IIa, V, Va, VII, VIIa, VIII, VIIIa, IX, IXa, X, Xa, 
	V_IIa, VII_IIa, VII_Xa, VIII_IIa, APC, TFPI, TFPIa, ADP //, TxA2
};
const int agonists[] = { 2, 1, 19 };					// number of agonists and their indices
*/

const double N2 = 2000., N5 = 3000., N8 = 450., N9 = 250., N9ast = 250., N10 = 2700.;
// setting the reaction rates ( 1/sec, 1/(nMsec) )
const double kz7see10f = 5e6*cscale*tscale, kz7see10r = 1.0*tscale, kz7see10cat = 5.0*tscale, kz7see2f = 3.92e5*cscale*tscale, kz7see2r = 1.0*tscale, kz7see2cat = 6.1e-2*tscale,
						 kz10e7sef = 8.95e6*cscale*tscale, kz10e7ser = 1.0*tscale, kz10e7secat = 1.15*tscale, kz9e7sef = 8.95e6*cscale*tscale, kz9e7ser = 1.0*tscale, kz9e7secat = 1.15*tscale,
						 k7on = 5e7*cscale*tscale, k7off = 5e-3*tscale, kz7e10f = 5e6*cscale*tscale, kz7e10r = 1.0*tscale, kz7e10cat = 5.0*tscale,
						 kz7e2f = 3.92e5*cscale*tscale, kz7e2r = 1.0*tscale, kz7e2cat = 6.1e-2*tscale, kz5e2f = 1.73e7*cscale*tscale, kz5e2r = 1.0*tscale, kz5e2cat = 0.23*tscale,
						 kz8e2f = 2.64e7*cscale*tscale, kz8e2r = 1.0*tscale, kz8e2cat = 0.9*tscale, k9on = 1e7*cscale*tscale, k9off = 2.5e-2*tscale, k10on = 1e7*cscale*tscale, k10off = 2.5e-2*tscale,
						 k5on = 5.7e7*cscale*tscale, k5off = 0.17*tscale, k8on = 5e7*cscale*tscale, k8off = 0.17*tscale, k2on = 1e7*cscale*tscale, k2off = 5.9*tscale,
						 kz5me10mf = 1e8*cscale*tscale, kz5me10mr = 1.0*tscale, kz5me10mcat = 4.6e-2*tscale, kz5me2mf = 1.73e7*cscale*tscale, kz5me2mr = 1.0*tscale, kz5me2mcat = 0.23*tscale,
						 kz8me10mf = 5.1e7*cscale*tscale, kz8me10mr = 1.0*tscale, kz8me10mcat = 2.3e-2*tscale, kz8me2mf = 2.64e7*cscale*tscale, kz8me2mr = 1.0*tscale, kz8me2mcat = 0.9*tscale,
						 kz10mtenf = 1.31e8*cscale*tscale, kz10mtenr = 1.0*tscale, kz10mtencat = 20.*tscale, kz2mprof = 1.03e8*cscale*tscale, kz2mpror = 1.0*tscale, kz2mprocat = 30.*tscale,
						 ktenf = 1e8*cscale*tscale, ktenr = 0.01*tscale, kprof = 1e8*cscale*tscale, kpror = 0.01*tscale, k9in = 0.1*tscale, k10in = 0.1*tscale, k2in = 0.2*tscale,
						 kapce8mf = 1.2e8*cscale*tscale, kapce8mr = 1.0*tscale, kapce8mcat = 0.5*tscale, kapce5mf = 1.2e8*cscale*tscale, kapce5mr = 1.0*tscale, kapce5mcat = 0.5*tscale,
						 ktfpiae10f = 1.6e7*cscale*tscale, ktfpiae10r = 3.3e-4*tscale, ktfpiae7sef = 1e7*cscale*tscale, ktfpiae7ser = 1.1e-3*tscale;

// Anand's Model
enum SPECIES {
	IXa, Ix, VIIIa, VIII, Va, V, Xa, X, IIa, II, Ia, I, XIa, XI,
	ATIII, TFPI, APC, PC, L1AT, ADP //, TxA2
};
const int agonists[] = { 2, 8, 19 };					// number of agonists and their indices

// setting the reaction rates ( 1/sec, 1/(nMsec), nM )
const double k9 = (11.0/60.0)*tscale, K9M = 160.0, h9 = (0.0162/60.0)*tscale, k8 = (194.4/60.0)*tscale, K8M = 1.12e5, h8 = (0.222/60.0)*tscale, 
					hC8 = (10.2/60.0)*tscale, HC8M = 14.6, k5 = (27.0/60.0)*tscale, K5M = 140.5, h5 = (0.17/60.0)*tscale, hC5 = (10.2/60.0)*tscale, 
					HC5M = 14.6, k10 = (2391.0/60.0)*tscale, K10M = 160.0, h10 = (0.347/60.0)*tscale, hTFPI = (0.48/60.0)*tscale, k2 = (1344.0/60.0)*tscale, 
					K2M = 1060.0, h2 = (0.714/60.0)*tscale, k1 = (3540.0/60.0)*tscale, K1M = 3160.0, k11 = (0.0078/60.0)*tscale, K11M = 50.0,
        	h11A3 = (1.6e-3/60.0)*tscale, h11L1 = (1.3e-5/60.0)*tscale, kPC = (39.0/60.0)*tscale, KPCM = 3190.0, hPC = (6.6e-7/60.0)*tscale, 
					KdZ = 0.56, KdW = 0.1, k79 = (32.4/60.0)*tscale, K79M = 24.0, k710 = (103.0/60.0)*tscale, K710M = 240.0;

//	TODO: these numbers must be finalized !!!!
const double weight[] = { 1.0, 1.0 };				// weights for activation by each specie
const double thresh[] = { 1.0, 1000.0 };	// activation threshold for each specie (nM)
const double relCon[] = { 0.0, 0.03/(stddev*sqrt(2.0*M_PI)) };		// release content (nM)
#endif

#ifdef CCLF
enum SPECIES {
	PLTmu, PLTma, IIa, II, ADP, TxA2, Ia, PLTba
};

const int 	 Ndis = 100;
const double Vscale = lscale*lscale*lscale;
const int 	 agonists[] = { 3, 2, 4, 5 };					// number of agonists and their indices
const double weight[] = { 1.0, 1.0, 1.0 };				// weights for activation by each specie
const double thresh[] = { 1.0, 1000.0, 600.0 };	// activation threshold for each specie (nM)
const double relCon[] = {0.0, 170.0/(stddev*sqrt(2.0*M_PI)), 0.0};		// release content (nmol)

const double kadh = 2.0e2*tscale, kcohxPLTmax = 1.0e4*tscale, kadp = 0.34*tscale, kiia = 0.5*tscale, ktxa2 = 1.0*tscale, k1 = (2.82*1.0e1/60.0)*tscale, ksurf = 5.0e-6*tscale, kiiAP = 1.0*tscale, kin = 0.025*tscale; // TODO: all physical units; must be converted !!!

const double PLTmax = 0.12*1.0e12*lscale*lscale, ksi = 0.05, g0 = 1.0/0.9988, etat = 0.05, etastar = pow(0.1, 3), phi0B = 0.5, tact = 1.0*1.0e-1/tscale;
#endif

#endif
#endif

#endif
