/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $Source: /homedir/cvs/Nektar/Nektar3d/src/prepost.C,v $  
 * $Revision: 1.3 $
 * $Date: 2006/05/08 10:01:39 $    
 * $Author: ssherw $  
 * $State: Exp $   
 *---------------------------------------------------------------------------*/
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <veclib.h>
#include "nektar.h"
#include "pbc_1d.h"

/* local functions */
static void     Summary     (void);
int             bandwidth   (Element *E, Bsystem *Bsys);
static int      suminterior (Element *E);

struct ssn_tag {
  char name[FILENAME_MAX];    /* name of session without postfix   */
  char fld [FILENAME_MAX];    /* Output (field) file               */
  char his [FILENAME_MAX];    /* History point file                */
  char rea [FILENAME_MAX];    /* Input file (from pre-processor)   */
#ifdef ADR
  char adr [FILENAME_MAX];    /* Input file (for temp. B.C.)	     */
#endif
  char fce [FILENAME_MAX];    /* Forces file                       */
  char flo [FILENAME_MAX];    /* Flow Rate file                    */
  char pre [FILENAME_MAX];    /* Pressure at outlets file          */
  char mom [FILENAME_MAX];    /* momentum flux at outlets/inlet file  */
} session;

/* function from flowrate.C */
double PressureMean(Domain *omega, char *label);

static struct {              /* Default options for Nekscal         */
  char *name;
  int   val;
  char *descrip;
} nektar_ops[] = {
  "order"     ,0, "-n #   ... run with # of modes ",
  "binary"    ,1, "-a     ... write output in ascii format",
  "verbose"   ,0, "-v     ... plot divergence error",
  "variable"  ,0, "-var   ... run with variable order (uses [].ord file)",
  "checkpt"   ,0, "-chk   ... dump checkpoint fields at \'iostep\' steps",
  "SurForce"  ,0, "-sfce  ... dump sectional forces",
  "oldhybrid" ,0, "-old   ... read old hybrid format field files",
  "NoVertBlock",0, "-nvb   ... No vertex block in LE Precon",
  "MRHS_NRHS", 0,"-mrhs # ... Multiple Right Hand side solver using # vectors",
  "iterative" ,0, "-i     ... iterative solve",
  "Initcond"  ,0, "-I     ... set up initial condition",
  "recursive" ,0, "-r #   ... recursive static condesation with # recursions",
  "timeavg"   ,0, "-T     ... evaluate time average field",
  "tvarying"  ,0, "-V     ... time varying calculation",
	"paraview"  ,0, "-vtk   ... dump vtu files at \'vtksteps\' steps",
	"SVVE"      ,0, "-svve  ... use explicit SVV",
	"SVVI"      ,0, "-svvi  ... use implicit SVV",
   0, 0, 0
};

void parse_args(int argc, char *argv[])
{
  int     i;
  char    c;

#ifdef DEBUG
  init_debug();
  // mallopt(M_DEBUG,0);
#endif

  /* initialize the parser and install the defaults */  
  manager_init();

  for (i = 0; nektar_ops[i].name; i++) 
    option_set(nektar_ops[i].name, nektar_ops[i].val);
  
  if(argc == 1) goto showUsage;

  while (--argc && (*++argv)[0] == '-') {
    if (strcmp (*argv+1, "chk") == 0) {
      option_set("checkpt",1);
      continue;
    }
    else if (strcmp (*argv+1, "sfce") == 0) {
       option_set("SurForce",1);
       continue;
    }
    else if (strcmp (*argv+1, "var") == 0) {
      option_set("variable",1);
      continue;
    }
    else if (strcmp (*argv+1, "old") == 0) {
      option_set("oldhybrid",1);
      continue;
    }
    else if (strcmp (*argv+1, "nvb") == 0) {
      option_set("NoVertBlock",1);
      continue;
    }
    else if (strcmp (*argv+1, "vtk") == 0) {
      *argv++;
      iparam_set("VTKSTEPS",atoi(*argv));
      argc--;
      continue;
    }
    else if (strcmp (*argv+1, "svve") == 0) {
      option_set("SVVE",1);
      continue;
    }
    else if (strcmp (*argv+1, "svvi") == 0) {
      option_set("SVVI",1);
      continue;
    }
    else if (strcmp (*argv+1, "mrhs") == 0) {
      int n; 
      argv[0]+=4;
      if (*++argv[0]) 
	n = atoi (*argv);
      else {
	n = atoi (*++argv);
	argc--;
      }
      option_set("MRHS_NRHS",n);
      (*argv)[1] = '\0'; 
      continue;
    }
    while (c = *++argv[0])
      switch (c) {
      case 'a':
	option_set("binary",0);
	break;
      case 'b':
	fprintf(stderr,"Option b is now automatic\n");
	break;
      case 'c':
	{ int n; 
	  if (*++argv[0]) 
	    n = atoi (*argv);
	  else {
	    n = atoi (*++argv);
	    argc--;
	  }
	  option_set("TORDER",n);
	  (*argv)[1] = '\0'; 
	break;
	}
      case 'i':
	option_set("iterative",1);
	break;
      case 'I': 
	option_set("Initcond",1);
	break;
      case 'm':
	option_set("mixediter",1);
	break;
      case 'n':
	{ int n; 
	  if (*++argv[0]) 
	    n = atoi (*argv);
	  else {
	    n = atoi (*++argv);
	    argc--;
	  }
	  option_set("NORDER.REQ",n);
	  (*argv)[1] = '\0'; 
	}
	break;
      case 'r':
	{ int n; 
	  if (*++argv[0]) 
	    n = atoi (*argv);
	  else {
	    n = atoi (*++argv);
	    argc--;
	  }
	  option_set("recursive",n);
	  (*argv)[1] = '\0'; 
	}
	break;
      case 'S':
	{
	  option_set("SLICES",1);
	  break;
	}
      case 't':
	{ double n; 
	  if (*++argv[0]) 
	    n = atof (*argv);
	  else {
	    n = atof (*++argv);
	    argc--;
	  }
	  dparam_set("THETA",n);
	  (*argv)[1] = '\0'; 
	}
	break;
      case 'T':
	option_set("timeavg",1);
	break;
      case 'v':
	option_set("verbose",2);
	break;
      case 'V':
	option_set("tvarying",1);
	break;
      default:
	goto showUsage;
      }
  }
  return;

 showUsage:
  ROOTONLY{
    fputs("usage: nektar [options] file[.rea]\n\n"
	  "options:\n", stderr);
    for (i = 0; nektar_ops[i].name; i++)
      fprintf(stderr, "%s\n", nektar_ops[i].descrip);
  }
  exit(-1);
}
  
void LocalNumScheme  (Element_List *E, Bsystem *Bsys, Gmap *gmap);
Gmap *GlobalNumScheme(Element_List *E, Bndry *Ebc);
void free_gmap(Gmap *gmap);
void free_Mesh_Bcs(Bndry *MeshBcs);
void Reflect_Global_Velocity(Element_List *Mesh, Bndry *Meshbcs, int dir);

Element_List *Mesh;
Element_List *MeshTemp;

Domain *PreProcess(int argc, char **argv)
{
  FILE     *rea_file;
  Element  *E;
  Bndry    *Meshbc;
  Domain   *omega;
  int       k,Je,lev;
  double    Re,dt;
  Gmap     *gmap;
  Gmap     *gmap_p;
  int       procid = mynode();
  SLVTYPE   slvtype;
#ifdef ADR
	FILE		 *adr_file;
  Bndry    *MeshbcTemp;	
  Gmap     *gmap_Temp;	
#endif

  /* Create the new domain and open the input files */
  parse_args(argc, argv);

  option_set("GSLEVEL", min(pllinfo.nprocs,8));
  option_set("FAMOFF" , 1);

  omega = (Domain*) calloc(1,sizeof(Domain));

  sprintf(session.name,"%s"     , strtok(argv[argc-1],"."));
  sprintf(session.rea, "%s.rea" , argv[argc-1]);
#ifdef ADR
  sprintf(session.adr, "%s.adr" , argv[argc-1]);
#endif
  sprintf(session.fld, "%s.fld" , argv[argc-1]);
  sprintf(session.his, "%s.his" , argv[argc-1]);

#ifdef DEBUG
  DO_PARALLEL{
    char *buf = (char*) calloc(BUFSIZ,sizeof(char));
    sprintf(buf, "%s.dbx.%d", argv[argc-1], procid);
    debug_out = fopen(buf, "w");
  }
#endif

  ROOTONLY{
    sprintf(session.fce, "%s.fce" , argv[argc-1]);
    sprintf(session.flo, "%s.flo" , argv[argc-1]);
    sprintf(session.pre, "%s.pre" , argv[argc-1]);
    sprintf(session.mom, "%s.mom" , argv[argc-1]);
  }

  omega->name     = argv[argc-1];
  rea_file        = fopen(session.rea,"r");
#ifdef ADR
  adr_file        = fopen(session.adr,"r");
#endif

  if (rea_file == (FILE*) NULL)
    error_msg(PreProcess(): Could not open input file(s));

  ROOTONLY{
    omega->fce_file = fopen(session.fce,"a");
    omega->flo_file = fopen(session.flo,"a");
    omega->pre_file = fopen(session.pre,"a");
    omega->mom_file = fopen(session.mom,"a");
  }

  /* Read the input parameters */
  ReadParams  (rea_file);
  ReadPscals  (rea_file, omega);
  ReadLogics  (rea_file);

  Je = iparam("INTYPE");
  slvtype = (SLVTYPE) iparam("SLVTYPE");

#if 0
    set_LZero(1);
    iparam_set("LQUAD",iparam("MODES")+1);
    iparam_set("MQUAD",iparam("MODES")+1);
    iparam_set("NQUAD",iparam("MODES")+1);
#endif
  if(slvtype == SubStep){
    set_LZero(1);
    iparam_set("LQUAD",iparam("MODES")+1);
    iparam_set("MQUAD",iparam("MODES")+1);
    iparam_set("NQUAD",iparam("MODES")+1);
    option_set("OpSplit",1);
  }
  
  if(option("tvarying"))
    dparam_set("t",0);

  /* Build the mesh */
  Mesh   = ReadMesh(rea_file, session.name);
  Meshbc = ReadMeshBCs(rea_file, Mesh);   
  gmap   = GlobalNumScheme (Mesh, Meshbc);
	
#ifdef ADR
  MeshbcTemp = ReadMeshBCsTemp(adr_file, Mesh); 
//  gmap_Temp  = GlobalNumScheme (Mesh, MeshbcTemp);
#endif
	
  ReadTimeDepend (rea_file, omega);

  /* set up velocity system using U */
  if(option("SurForce"))
   omega->U    = LocalMesh(Mesh,Meshbc,session.name);
  else
   omega->U    = LocalMesh(Mesh,session.name);
  
  if(slvtype == SubStep){
    init_ortho_basis();
    for(E = omega->U->fhead;E; E = E->next) E->dgL = E->lmax;
    set_elmt_sides   (omega->U);
    set_sides_geofac (omega->U);
  }

  DO_PARALLEL{ // recall global numbering to put partition vertices first
    free_gmap(gmap);
    Reflect_Global_Velocity    (Mesh, Meshbc, 0);
    gmap  = GlobalNumScheme    (Mesh, Meshbc);
    Replace_Numbering          (omega->U, Mesh);
  }

  omega->Ubc  = ReadBCs (rea_file,omega->U->fhead);
  /* also sets up edge links system */
  omega->Usys = gen_bsystem(omega->U,gmap);
  free_gmap(gmap);
  
  /* set up velocity mat structure */
  ROOTONLY Summary ();

  /* set up V field */
  omega->V   = omega->U->gen_aux_field ('v');
  if(slvtype == SubStep) set_elmt_sides(omega->V);
  omega->Vbc = ReadBCs       (rea_file,omega->V->fhead);   
  if(option("REFLECT1")||option("REFLECT0")){
    Reflect_Global_Velocity    (Mesh, Meshbc, 1);
    gmap  = GlobalNumScheme    (Mesh, Meshbc);
    Replace_Numbering          (omega->V, Mesh);
    omega->Vsys = gen_bsystem  (omega->V, gmap);
    free_gmap(gmap);
  }
  else
    omega->Vsys = omega->Usys;
  
  /* set up W field */
  omega->W   = omega->U->gen_aux_field ('w');
  if(slvtype == SubStep) set_elmt_sides(omega->W);
  omega->Wbc = ReadBCs(rea_file,omega->W->fhead);   
  if(option("REFLECT2")||
     (option("REFLECT0")&&option("REFLECT1")&&(!option("REFLECT2")))){
    Reflect_Global_Velocity   (Mesh, Meshbc, 2);
    gmap  = GlobalNumScheme   (Mesh, Meshbc);
    Replace_Numbering         (omega->W, Mesh);
    omega->Wsys = gen_bsystem (omega->W, gmap);
    free_gmap(gmap);
  }
  else if (!option("REFLECT0"))
    omega->Wsys = omega->Usys;
  else if (!option("REFLECT1")){
    Replace_Local_Numbering(omega->W,omega->V);
    omega->Wsys = omega->Vsys;
  }

  /* set up pressure system */
  Set_Global_Pressure    (Mesh, Meshbc);
  gmap = GlobalNumScheme (Mesh, Meshbc);
  omega->P            = omega->U->gen_aux_field ('p');    
  Replace_Numbering                 (omega->P, Mesh);
  omega->Pbc          = BuildPBCs   (omega->P, omega->Ubc);
  omega->Pressure_sys = gen_bsystem (omega->P, gmap);
  free_gmap(gmap);

#ifdef ADR
	int nspecs = iparam("NSPEC");
	omega->T = (Element_List **) malloc(nspecs*sizeof(Element_List *));
	omega->Tf = (Element_List **) malloc(nspecs*sizeof(Element_List *));
	omega->Tbc = (Bndry **) malloc(nspecs*sizeof(Bndry *));
	omega->Tsys = (Bsystem **) malloc(nspecs*sizeof(Bsystem *));

	Reflect_Global_Velocity(Mesh, MeshbcTemp, 0);
  gmap_Temp				= GlobalNumScheme(Mesh, MeshbcTemp);
	for (int i = 0; i < nspecs; i++){
		omega->T[i]			= omega->U->gen_aux_field('t');
		Replace_Numbering(omega->T[i], Mesh);
	  if(slvtype == SubStep) set_elmt_sides(omega->T[i]);
  	omega->Tbc[i]		= ReadBCs(adr_file, omega->T[i]->fhead);   
  	omega->Tsys[i]	= gen_bsystem(omega->T[i], gmap_Temp);
	}
  free_gmap(gmap_Temp);

#if 0
	for (int i = 0; i < nspecs; i++){
	  /* set up temperature system */
  	omega->T[i]    = omega->U->gen_aux_field ('t');
		DO_PARALLEL{
			free_gmap(gmap_Temp);    
			Reflect_Global_Velocity    (Mesh, MeshbcTemp, 0);
			gmap_Temp  = GlobalNumScheme    (Mesh, MeshbcTemp);
			Replace_Numbering          (omega->T[i], Mesh);
		}
		if(slvtype == SubStep) set_elmt_sides(omega->T[i]);
		omega->Tbc[i]  = ReadBCs (adr_file, omega->T[i]->fhead);
		omega->Tsys[i] = gen_bsystem(omega->T[i], gmap_Temp);
	}
	free_gmap(gmap_Temp);
#endif
#endif

  /* set up forcing function storage */
  omega->Pf = omega->P->gen_aux_field('p');
  omega->Uf = omega->U->gen_aux_field('u');
  omega->Vf = omega->V->gen_aux_field('v');
  omega->Wf = omega->W->gen_aux_field('w');
#ifdef ADR
	for (int i = 0; i < nspecs; i++)
	  omega->Tf[i] = omega->T[i]->gen_aux_field('t');	
#endif
  
  if(slvtype == SubStep){
    set_elmt_sides(omega->Uf);
    set_elmt_sides(omega->Vf);
    set_elmt_sides(omega->Wf);
#ifdef ADR
		for (int i = 0; i < nspecs; i++)
			set_elmt_sides(omega->Tf[i]);
#endif
  }

  dt = dparam("DELT");
  omega->soln = NULL;
  ReadSoln(rea_file, omega);

  Re = 1.0/dparam("KINVIS");
  //  ReadKinvis (omega);  

  omega->Pressure_sys->lambda = (Metric*)calloc(omega->P->nel, sizeof(Metric));

  ROOTONLY{
    fprintf(stdout,"Generating pressure system [."); 
    fflush(stdout);
  }
  GenMat (omega->P,omega->Pbc,omega->Pressure_sys,
	  omega->Pressure_sys->lambda,Helm);
  ROOTONLY{
    fprintf(stdout,"]\n");
    fflush(stdout);
  }
    
  omega->Usys->lambda = (Metric*) calloc(omega->U->nel, sizeof(Metric));
  
  for(k=0;k<omega->U->nel;++k)
    omega->Usys->lambda[k].d = Re*getgamma(1)/dt;

  omega->Vsys->lambda = omega->Usys->lambda;
  omega->Wsys->lambda = omega->Usys->lambda;

#ifdef ADR
	for (int i = 0; i < nspecs; i++){
  	omega->Tsys[i]->lambda = (Metric*) calloc(omega->T[i]->nel, sizeof(Metric));

  	for(k=0;k<omega->T[i]->nel;++k)
			omega->Tsys[i]->lambda[k].d = Re*omega->Pr[i]*getgamma(1)/dt;  
	}
#endif

  ROOTONLY {
    fprintf(stdout,"Generating velocity system [."); 
    fflush(stdout);
  }

  GenMat (omega->U,omega->Ubc,omega->Usys,omega->Usys->lambda,Helm);
  
  ROOTONLY {
    fprintf(stdout,"]\n"); 
    fflush(stdout);
  }

  if(option("REFLECT1")||option("REFLECT0")){
    ROOTONLY {
      fprintf(stdout,"Generating V-velocity system [."); 
      fflush(stdout);
    }
    if(omega->Vsys->smeth == iterative){
      ROOTONLY fprintf(stdout,"..re-using u-sys..");
      omega->Vsys->Gmat = omega->Usys->Gmat; // use previous declaration
    }
    GenMat (omega->V,omega->Vbc,omega->Vsys,omega->Vsys->lambda,Helm);
    ROOTONLY {
      fprintf(stdout,"]\n"); 
      fflush(stdout);
    }
  }

  if(option("REFLECT2")||
     (option("REFLECT0")&&option("REFLECT1")&&(!option("REFLECT2")))){
    ROOTONLY {
      fprintf(stdout,"Generating W-velocity system [."); 
      fflush(stdout);
    }
    if(omega->Wsys->smeth == iterative){
      ROOTONLY fprintf(stdout,"..re-using u-sys..");
      omega->Wsys->Gmat = omega->Usys->Gmat; // use previous declaration
    }
    GenMat (omega->W,omega->Wbc,omega->Wsys,omega->Wsys->lambda,Helm);
    ROOTONLY {
      fprintf(stdout,"]\n"); 
      fflush(stdout);
    }
  }

#ifdef ADR
  ROOTONLY {
		fprintf(stdout,"Generating temperature system [."); 
		fflush(stdout);
	}
	for (int i = 0; i < nspecs; i++)	
		GenMat (omega->T[i],omega->Tbc[i],omega->Tsys[i],omega->Tsys[i]->lambda,Helm);
	
	ROOTONLY {
		fprintf(stdout,"]\n"); 
		fflush(stdout);
	}
#endif	

#ifndef PARALLEL
  if(!option("tvarying")){
    DirBCs(omega->U,omega->Ubc,omega->Usys,Helm);
    DirBCs(omega->V,omega->Vbc,omega->Vsys,Helm);
    DirBCs(omega->W,omega->Wbc,omega->Wsys,Helm);
	  
    omega->U->zerofield();
    omega->V->zerofield();
    omega->W->zerofield();	 
#ifdef ADR
		for (int i = 0; i < nspecs; i++){
			DirBCs(omega->T[i],omega->Tbc[i],omega->Tsys[i],Helm);
			omega->T[i]->zerofield();
		}
#endif
  }
#endif

  if(omega->Tfun && omega->Tfun->type == Womersley) WomInit(omega);
	omega->sources = NULL;

  ReadICs     (rea_file, omega);
  ReadDF      (rea_file, omega);
  ReadHisData (rea_file, omega);
  ReadDFunc   (rea_file, omega);

  DO_PARALLEL
    sprintf(session.his,"%s.%d",session.his,pllinfo.procid);

	if (omega->his_list)
  	omega->his_file = fopen(session.his,"w");
  ROOTONLY omega->fce_file = fopen(session.fce,"a");
  ROOTONLY omega->flo_file = fopen(session.flo,"a");
  ROOTONLY omega->pre_file = fopen(session.pre,"a");
  ROOTONLY omega->mom_file = fopen(session.mom,"a");

  // set up initial BCs - needs to be called after readics 
  if(omega->Tfun && omega->Tfun->type == Womersley)
    SetWomSol(omega,dparam("STARTIME"),0);
  else if(option("tvarying") == 1){
    Bndry *Bc;
    dparam_set("t",dparam("STARTIME"));
    for(Bc=omega->Ubc;Bc;Bc=Bc->next)
      Bc->elmt->update_bndry(Bc,1);
    for(Bc=omega->Vbc;Bc;Bc=Bc->next)
      Bc->elmt->update_bndry(Bc,1);
    for(Bc=omega->Wbc;Bc;Bc=Bc->next)
      Bc->elmt->update_bndry(Bc,1);
#ifdef ADR
  	for (int i = 0; i < nspecs; i++)
    	for(Bc=omega->Tbc[i];Bc;Bc=Bc->next)
      	Bc->elmt->update_bndry(Bc,1);
#endif
  }
  // free edges, vertices
  delete        (Mesh);
  free_Mesh_Bcs (Meshbc);
#ifdef ADR
  free_Mesh_Bcs (MeshbcTemp);
#endif
  
  // Set up multistep storage 
  if(slvtype == SubStep){
    omega->ul = dmatrix(0,Je-1,0,omega->U->htot-1);
    dzero(Je*omega->U->htot, *omega->ul, 1);
    omega->vl = dmatrix(0,Je-1,0,omega->V->htot-1);
    dzero(Je*omega->V->htot, *omega->vl, 1);
    omega->wl = dmatrix(0,Je-1,0,omega->W->htot-1);
    dzero(Je*omega->V->htot, *omega->wl, 1);
    lev = Je+1;
  }
  else
    lev = Je;

  omega->u  = dmatrix(0,lev-1,0,omega->U->htot-1);
  dzero(lev*omega->U->htot, *omega->u, 1);
  omega->uf = dmatrix(0,lev-1,0,omega->U->htot-1);
  dzero(lev*omega->U->htot, *omega->uf, 1);

  omega->v  = dmatrix(0,lev-1,0,omega->U->htot-1);
  dzero(lev*omega->U->htot, *omega->v, 1);
  omega->vf = dmatrix(0,lev-1,0,omega->U->htot-1);
  dzero(lev*omega->U->htot, *omega->vf, 1);

  omega->w  = dmatrix(0,lev-1,0,omega->U->htot-1);
  dzero(lev*omega->U->htot, *omega->w, 1);
  omega->wf = dmatrix(0,lev-1,0,omega->U->htot-1);
  dzero(lev*omega->U->htot, *omega->wf, 1);

#ifdef ADR
	omega->t = (double ***) malloc(nspecs*sizeof(double **));
	omega->tf = (double ***) malloc(nspecs*sizeof(double **));
	for (int i = 0; i < nspecs; i++){
  	omega->t[i]  = dmatrix(0,lev-1,0,omega->T[i]->htot-1);
	  dzero(lev*omega->T[i]->htot, *omega->t[i], 1);
  	omega->tf[i] = dmatrix(0,lev-1,0,omega->T[i]->htot-1);
  	dzero(lev*omega->T[i]->htot, *omega->tf[i], 1);	
	}
#ifdef CCLF
	omega->ind_func = dvector(0,omega->U->htot-1);
	omega->number = ivector(0,omega->U->nel-1);
	izero(omega->U->nel, omega->number, 1);
	omega->volume = dvector(0,omega->U->nel-1);
	dzero(omega->U->nel, omega->volume, 1);
	omega->timeAct = dmatrix(0,omega->U->nel-1,0,Ndis-1);
#endif
#endif

 /* this needs to be set after LocalMesh so that pllinfo is defined */
#if 0		//Alireza: No need to write the field files at the end of simulation
  DO_PARALLEL{
    sprintf(session.fld, "%s.fld.hdr.%d", omega->name,pllinfo.procid);
    omega->fld_file = fopen(session.fld,"w");
    sprintf(session.fld, "%s.fld.dat.%d", omega->name,pllinfo.procid);
    omega->dat_file = fopen(session.fld,"w");
  }
  else{
    sprintf(session.fld, "%s.fld", omega->name);
    omega->fld_file = fopen(session.fld,"w"); 
    omega->dat_file = omega->fld_file;
 }    
#endif

#ifdef PBC_1D
  omega->IMPBC = new PBC1D[1];
  omega->IMPBC[0].Create(omega->Ubc);
  omega->IMPBC[0].ReadFlowrateHistory(omega->name);
  omega->IMPBC[0].SetGeofac(omega->Ubc, omega->Vbc, omega->Wbc);

  gsync();
  ROOTONLY
    fprintf(stdout, "mega->IMPBC[0].setGeofac - done\n");

  if (omega->IMPBC[0].Nimpedance_modes == 0){
		if (iparam("ILPN"))
			omega->IMPBC[0].ReadLPN(rea_file, omega);
		else
	    omega->IMPBC[0].SetRC(omega->name);

    int ii;
    double *Pressure_temp;
    Pressure_temp = dvector(0,omega->IMPBC[0].Nout);
    for (ii = 0; ii < omega->IMPBC[0].Nout; ++ii)
      omega->IMPBC[0].Pressure[ii] = PressureMean(omega,omega->IMPBC[0].standard_labels[ii])/omega->IMPBC[0].Area[ii];
    gdsum(omega->IMPBC[0].Pressure,omega->IMPBC[0].Nout,Pressure_temp);

    for (ii = 0; ii < omega->IMPBC[0].Nout; ++ii)
      omega->IMPBC[0].flowrate_RCR_old[ii] = 0.0;
    gdsum(omega->IMPBC[0].flowrate_RCR_old,omega->IMPBC[0].Nout,Pressure_temp);

    free(Pressure_temp);
  }
#endif

  fclose(rea_file);
#ifdef ADR
  fclose(adr_file);
#endif

  return omega;
}

void PostProcess(Domain *omega, int step, double time){
  Element_List **V;
  FILE  *fp[2];
  int nfields;
  nfields = omega->U->fhead->dim()+1;

#if 0		//Alireza: No need to write the field files at the end of simulation
  fp[0] = omega->fld_file;
  fp[1] = omega->dat_file;


  V = (Element_List**) malloc(nfields*sizeof(Element_List*));
  
  V[0] = omega->U;
  V[1] = omega->V;
  V[2] = omega->W;
  V[3] = omega->P;
  

  Writefld(fp, omega->name, step, time, nfields, V);
    
  DO_PARALLEL{
    fclose(fp[0]);
    fclose(fp[1]);
  }
  else
    fclose(fp[0]);
#endif

	 if ( (omega->his_list) )
  	fclose(omega->his_file);
   ROOTONLY  fclose(omega->fce_file);
   ROOTONLY  fclose(omega->flo_file);
   ROOTONLY  fclose(omega->pre_file);
   ROOTONLY  fclose(omega->mom_file);
}

/* ----------------------------------------------------------------------- *
 * Summary() - Print a summary of the input data                           *
 *                                                                         *
 * This function collects the echo of the input data which used to be      *
 * scattered throughout this file.                                         *
 *                                                                         *
 * ----------------------------------------------------------------------- */

static void Summary (void)
{
  printf ("Input File          : %s\n",session.name);
  printf (" Concat files       : %g\n",dparam("concat"));
  printf ("Re (U=%lg, L=%lg)       : %g\n",dparam("Re_Uinf"),dparam("Re_Len"), dparam("Re_Uinf")*dparam("Re_Len")/dparam("KINVIS"));
#ifdef OSSCYL
  printf ("KC   number         : %g\n",dparam("KC"));
  printf ("beta number         : %g\n",1./(dparam("KC")*dparam("KINVIS")));
#endif
  printf ("Time step           : %g\n",dparam("DELT"));
  printf ("Integration order   : %d\n",iparam("INTYPE"));
    if(!option("variable"))
     printf("Number of modes     : %d\n",iparam("MODES"));
  else
     printf("Number of modes     : variable\n");
  printf ("Number of elements  : %d\n",iparam("ELEMENTS"));
  printf ("Number of Families  : %d\n",iparam("FAMILIES"));
  DO_PARALLEL{
    printf ("Number of Processors: %d\n",pllinfo.nprocs);
    printf ("Gs level            : %d\n",option("GSLEVEL"));
  }
  fputs ("Equation type       : ", stdout);
  switch (iparam("EQTYPE")) {
  case Rotational:
    puts("Navier-Stokes (rotational)");
    break;
  case Convective:
    puts("Navier-Stokes (convective)");
    break;
  case Stokes:
    puts("Stokes flow");
    break;
  default:
    puts("undefined");
    break;
  }

  printf("Integration time    : %g, or %d steps\n", 
	 dparam("FINTIME"), iparam("NSTEPS"));
  printf("I/O time for saves  : %g, or %d steps",
	 dparam("IOTIME"),  iparam("IOSTEP"));
  
  if  (option("checkpt"))
    fputs (" [checkpoint]\n", stdout);
  else { 
    putchar ('\n');
    if(iparam("NSTEPS"))
      if (iparam("NSTEPS") / iparam("IOSTEP") > 10)
	fputs ("Summary: " "You have more than 10 dumps..."
		      "did you want to use -chk?\n", stderr);
  } 


  printf("NS Scheme Type      : ");
  switch(iparam("SLVTYPE")){
  case Splitting:
    printf("Splitting scheme \n");
    break;
  case SubStep:
    printf("DG-Advection Substepping  \n");
    break;
  default:
    printf("Unknown Solver \n");
  }

  if(option("iterative")){
    int n;
    if(n = option("MRHS_NRHS"))
      printf("Solver Type         : Iterative  (mrhs = %d) ",n);
    else
      printf("Solver Type         : Iterative ");
    
    switch((int) dparam("PRECON")){
    case Pre_Diag:
      printf("(Precon = diagonal)\n");
      break;
    case Pre_Block:
      printf("(Precon = block)\n");
      break;
    case Pre_None:
      printf("(Precon = none)\n");
      break;
    case Pre_LEnergy:
      if(option("NoVertBlock"))
	printf("(Precon = low energy  block (No Vert Block) )\n");
      else
	printf("(Precon = low energy  block)\n");
      break;
    }
  }
  else{
    printf("Solver Type         : Direct ");
    
    if(option("recursive"))
      printf("Mutlilevel \n");
    else
      printf("RCM \n");
  }

  return;
}

/* This is a function to sort out the global numbering scheme for      *
 * vertices, edges and faces. It also sets up the list of cumulative   *
 * indices for edges and vertices.                                     */

static void set_bmap       (Element *, Bsystem *);

Bsystem *gen_bsystem(Element_List *UL, Gmap *gmap){
  Bsystem  *Bsys;
  Element  *E=UL->fhead;

  Bsys       = (Bsystem *)calloc(1,sizeof(Bsystem));

  DO_PARALLEL /* force iterative solver if in parallel */
    option_set("iterative",1);

  if(option("iterative")){
    Bsys->smeth = iterative;
    Bsys->Precon = (enum PreType) (int) dparam("PRECON");
  }

  /* basis numbering scheme for vertices, edges and faces */
  LocalNumScheme(UL,Bsys,gmap);

  if(option("recursive"))
    Mlevel_SC_decom(UL,Bsys);
  else  if(Bsys->smeth == direct){
    bandwidthopt(E,Bsys,'a');
    ROOTONLY{
      fprintf(stdout,"rcm bandwidth (%c)   : %d [%d (%d)] \n",E->type,
	      bandwidth(E,Bsys),Bsys->nsolve,suminterior(E));
      fflush(stdout);
    }
  }

  set_bmap(E,Bsys);

  Bsys->families = iparam("FAMILIES");
  return Bsys;
}

static void  setGid(Element_List *E);

Gmap *GlobalNumScheme(Element_List *UL, Bndry *Ebc){
  Element  *U = UL->fhead;
  register int i;
  int       nvg,nvs,neg,nes,nfg,nfs,scnt,ncnt;
  int      *gsolve,*gbmap,l,l1;
  Bndry    *Ubc; 
  Vert     *v;
  Edge     *e;
  Face     *f;
  Gmap     *gmap;
  Element *E;

  gmap = (Gmap *)malloc(sizeof(Gmap));

  setGid (UL);  /* setup a global numbering scheme i.e. without boundaries*/

  /* This part of the routine re-orders the vertices, edges and then
     the faces so that the knowns are listed first. For the edges and
     the faces it also initialises a cummalative list stored in Bsys. */
  
  /*--------------------*/
  /* Vertex re-ordering */
  /*--------------------*/

  /* find maximum number of global vertices; */
  for(E=U,nvg=0; E; E = E->next)
    for(i = 0; i < E->Nverts; ++i)
      nvg = max(E->vert[i].gid, nvg);
  ++nvg;

  gsolve = ivector(0,nvg-1);
  gbmap  = ivector(0,nvg-1);

  ifill(nvg, 1, gsolve,1);
  // Assemble vertex solve mask to sort out multiplicity issues 
  for(E=U; E; E = E->next)
    for(i = 0; i < E->Nverts; ++i) // note vert[i].solve can have mag of 1 or 2
      gsolve[E->vert[i].gid] = (gsolve[E->vert[i].gid]*E->vert[i].solve)? 
	max(gsolve[E->vert[i].gid],E->vert[i].solve):0;
       //gsolve[E->vert[i].gid] &= E->vert[i].solve;
  
  // copy back mask 
  for(E=U; E; E = E->next){
    v = E->vert;
    for(i = 0; i < E->Nverts; ++i)
     v[i].solve = gsolve[E->vert[i].gid];
  }
  
  scnt = 0; ncnt = 0;

  DO_PARALLEL{
    /* reset vertex solve id ordering so that points on parallel 
       patches are first  */
    for(i = 0; i < nvg; ++i)
      if(gsolve[i] == 2)
	gbmap[i] =  scnt++;
      else if(gsolve[i] == 0)
	gbmap[i] = ncnt++;
    
    gmap->nv_psolve = scnt;

    for(i = 0; i < nvg; ++i)
      if(gsolve[i] == 1)
	gbmap[i] =  scnt++;
  }
  else{
    for(i = 0; i < nvg; ++i)
      gbmap[i] = gsolve[i]? scnt++:ncnt++;
  }
  nvs = scnt;

  /* place unknowns at the end of the pile */
  for(i = 0; i < nvg; i++) 
    gbmap[i] += gsolve[i]?  0:nvs;

  /* replace vertices numbering into vertex structures */
  for(E=U; E; E = E->next)
    for(i = 0; i < E->Nverts; ++i)
      E->vert[i].gid = gbmap[E->vert[i].gid];


  /* optimise vertex ordering for iterative solver */
  if(option("iterative")&&0){ // turn off form moment since going to use static
    Bsystem B;                // condensation on vertex solve 

    /* set up dummy B system */
    B.nel = UL->nel;
    B.nv_solve = nvs;
    bandwidthopt(UL->fhead,&B,'v');
  }

  free(gsolve); free(gbmap);

  /*------------------*/
  /* Edge re-ordering */
  /*------------------*/

  /* find maximum number of global edges; */
  for(E = U,neg=0; E; E = E->next)
    for(i = 0; i < E->Nedges; ++i)
      neg = max(E->edge[i].gid,neg);
  ++neg;

  gsolve = ivector(0,neg-1);
  gbmap  = ivector(0,neg-1);

  /* form edge gsolve and gbmap */
  ifill(neg, 1, gsolve,1);
  for(Ubc = Ebc; Ubc; Ubc = Ubc->next)
    if((Ubc->type == 'V')||(Ubc->type == 'W')||(Ubc->type == 'o')){
      E = Ubc->elmt;
      l = Ubc->face;
      if(E->dim() == 2)
	gsolve[E->edge[l].gid] = 0;
      else 
	for(i = 0; i < E->Nfverts(l); ++i)
	  gsolve[E->edge[E->ednum(l,i)].gid] = 0;
    }
  
  scnt = 0;  ncnt = 0;
  for(i = 0; i < neg; ++i)
    gbmap[i] = gsolve[i]? (scnt++):(ncnt++);
  nes   = scnt;

  /* place unknowns at the end of the pile */
  for(i = 0; i < neg; ++i)
    gbmap[i] += gsolve[i]?  0:nes;

  /* replace sort gid's */
  for(E = U; E; E = E->next)
    for(i = 0; i < E->Nedges; ++i)
      E->edge[i].gid = gbmap[E->edge[i].gid];
  
  if(U->dim() == 3){
    /*------------------*/
    /* Face re-ordering */
    /*------------------*/
    
    /* find maximum number of global faces; */
    for(E = U, nfg = 0; E; E = E->next)
      for(i = 0; i < E->Nfaces; ++i)
	nfg = (E->face[i].gid > nfg)? E->face[i].gid:nfg;
    ++nfg;
    
    gsolve = ivector(0,nfg-1);
    gbmap  = ivector(0,nfg-1);
    
    /* form faces part of gsolve, gbmap */
    ifill(nfg, 1, gsolve,1);
    for(Ubc = Ebc; Ubc; Ubc = Ubc->next)
      if((Ubc->type == 'V')||(Ubc->type == 'W')||(Ubc->type == 'o'))
	gsolve[Ubc->elmt->face[Ubc->face].gid] = 0;
    
    scnt = 0; ncnt = 0;
    for(i = 0; i < nfg; ++i)
      gbmap[i] = gsolve[i]? (scnt++):(ncnt++);
    nfs = scnt;
    
    /* place unknowns at the end of the pile */
    for(i = 0; i < nfg; ++i)
      gbmap[i] += gsolve[i]?  0:nfs;
    
    /* replace sorted gid's */
    for(E=U;E;E=E->next){
      f=E->face;
      for(i = 0; i < E->Nfaces; ++i)
	f[i].gid = gbmap[f[i].gid];
    }
    
    free(gsolve); free(gbmap);
  }


 /* finally set up a global id information  */
  gmap->nvs = nvs;
  gmap->nvg = nvg;
  gmap->nes = nes;
  gmap->neg = neg;
  gmap->nfs = nfs;
  gmap->nfg = nfg;

  gmap->vgid = ivector(0,nvg-1);
  gmap->egid = ivector(0,neg-1);
  gmap->fgid = ivector(0,nfg-1);

  int *elen  = ivector(0,neg-1);
  int *flen  = ivector(0,nfg-1);

  /* fill list with size */
  for(E=U;E;E=E->next){
    for(i = 0; i < E->Nedges; ++i)
      elen[E->edge[i].gid] = E->edge[i].l;

    for(i = 0; i < E->Nfaces; ++i)
      flen[E->face[i].gid] = (E->Nfverts(i) == 3) ?
	(E->face[i].l)*(E->face[i].l+1)/2 :
      	(E->face[i].l)*(E->face[i].l);
  }
  
  /* set up global numbering scheme */
  int cnt = 0;
  for(i = 0; i < nvs; ++i)
    gmap->vgid[i] = cnt++;

  for(i = 0; i < nes; cnt += elen[i],++i)
    gmap->egid[i] = cnt;

  for(i = 0; i < nfs;  cnt +=flen[i],++i)
    gmap->fgid[i] = cnt;

  gmap->nsolve = cnt;
  
  for(i = nvs; i < nvg; ++i)
    gmap->vgid[i] = cnt++;

  for(i = nes; i < neg; cnt += elen[i],++i)
    gmap->egid[i] = cnt;

  for(i = nfs; i < nfg;  cnt +=flen[i],++i)
    gmap->fgid[i] = cnt;

  gmap->nglobal = cnt;

  free(elen); 
  free(flen);
  
  return gmap;

}

void free_gmap(Gmap *gmap){
  free(gmap->vgid);
  free(gmap->egid);
  free(gmap->fgid);
  free(gmap);  
}

void LocalNumScheme(Element_List *EL, Bsystem *Bsys, Gmap *gmap){
  register  int i,k;
  int       nvk,nvs,nek,nes,ncnt,cnt;
  int       *vmap,*emap,*elen;
  int       nel = EL->nel; /* local elements */
  int       nfk,nfs;
  int      *fmap,*flen;
  Element   *E;

  /* This part of the routine re-orders the vertices, edges and then
     the faces so that the knowns are listed first. For the edges and
     the faces it also initialises a cummalative list stored in Bsys. */

  /*--------------------*/
  /* Vertex re-ordering */
  /*--------------------*/

  /* find the number of solved and global vertices in local mesh; */

  vmap  = ivector(0,gmap->nvg-1);
  ifill(gmap->nvg,-1,vmap,1);

  for(E=EL->fhead;E;E=E->next)
    for(i = 0; i < E->Nverts; ++i)
      vmap[E->vert[i].gid] = 1;

  /* count the number of locally solved vertices and set up mapping */
  nvs = 0;
  for(i = 0; i < gmap->nvs; ++i)
    if(vmap[i]+1) vmap[i] = nvs++;
  
  /* count the number of locally known vertices and set up mapping */
  /* Note: haven't added nvs onto vmap here since need to add nsolve
     at end of routine. */
  nvk = 0;
  for(i = gmap->nvs; i < gmap->nvg; ++i)
    if(vmap[i]+1) vmap[i] = nvk++;
  
  /*--------------------*/
  /* Edge re-ordering   */
  /*--------------------*/

  /* find the number of solved and known edges in local mesh; */
  emap = ivector(0,gmap->neg-1);
  elen = ivector(0,gmap->neg-1);
  ifill(gmap->neg,-1,emap,1);
  ifill(gmap->neg,-1,elen,1);
  
  for(E=EL->fhead;E;E=E->next)
    for(i = 0; i < E->Nedges; ++i){
      emap[E->edge[i].gid] = 1;
      elen[E->edge[i].gid] = E->edge[i].l;
    }

  /* count the number of locally solved edges and make up cummalative
     edge list */
  nes  = 0;
  for(i = 0; i < gmap->nes; ++i)
    if(emap[i]+1)  emap[i]  = nes++;

  /* count the number of locally known edges and make up cummalative 
     edge list*/
  nek = 0;
  for(i = gmap->nes; i < gmap->neg; ++i)
    if(emap[i]+1) emap[i]  = nes + nek++;
  
  /* set up cumulative edge list */
  Bsys->edge = ivector(0,nes+nek);
  ncnt = 0;
  for(i = 0, cnt = 0; i < gmap->neg; ++i){
    if(elen[i]+1){
      Bsys->edge[cnt++] = ncnt;
      ncnt += elen[i];
    }
  }
  Bsys->edge[cnt] = ncnt; /* this puts total edge length in final value */

  /*--------------------*/
  /* Face re-ordering   */
  /*--------------------*/

  /* find the number of solved and known faces  in local mesh; */

  fmap  = ivector(0,gmap->nfg-1);
  flen  = ivector(0,gmap->nfg-1);
  ifill(gmap->nfg,-1,fmap,1);
  ifill(gmap->nfg,-1,flen,1);
  
  for(E=EL->fhead;E;E=E->next)
    for(i = 0; i < E->Nfaces; ++i){
      fmap[E->face[i].gid] = 1;
      flen[E->face[i].gid] = (E->Nfverts(i) == 3) ? 
	E->face[i].l*(E->face[i].l+1)/2 : E->face[i].l*E->face[i].l;
    }
  
  /* count the number of locally solved faces and set up cummalative
     face list */
  nfs  = 0;
  for(i = 0; i < gmap->nfs; ++i)
    if(fmap[i]+1) fmap[i] = nfs++;
  
    
  /* count the number of locally known vertices and set up mapping */
  nfk = 0;
  for(i = gmap->nfs; i < gmap->nfg; ++i)
    if(fmap[i]+1) fmap[i] =  nfs + nfk++;

  /* set up cumulative face list */
  Bsys->face = ivector(0,nfs+nfk);
  ncnt = 0;
  for(i = 0,cnt= 0; i < gmap->nfg; ++i){
    if(flen[i]+1){
      Bsys->face[cnt++] = ncnt;
      ncnt += flen[i];
    }
  }
  Bsys->face[cnt] = ncnt; /* this puts total edge length in final value */
    
  /* work out total number of local solve values */
  Bsys->nv_solve  = nvs;
  Bsys->ne_solve  = nes;
  Bsys->nf_solve  = nfs;

  Bsys->nel       = EL->nel;
  Bsys->nsolve    = nvs + Bsys->edge[nes] + Bsys->face[nfs];
  Bsys->nglobal   = nvs + nvk + Bsys->edge[nes+nek] + Bsys->face[nfs+nfk];


  /* We want a local ordering scheme where the solved vertices are
     first followed by the solved edges and then solved faces. We then
     wish to order the known vertices followed by the known edges and
     finally the known face which is done below */
  
  /* Add nsolve to known vertices */
  for(i = gmap->nvs; i < gmap->nvg; ++i)
    if(vmap[i]+1) vmap[i] += Bsys->nsolve;
  
  /* Sort out values in cumalative list of face and edges. If faces
     are used these need to be done first because of requried
     information stored in edge list */

  int nstot  = Bsys->nsolve;

  int nvstot = nvs;
  int nestot = Bsys->edge[nes];
  int nfstot = Bsys->face[nfs];

  int nvktot = nvk;
  int nektot = Bsys->edge[nes+nek]-Bsys->edge[nes];

  for(i = 0; i < nes; ++i)
    Bsys->edge[i] += nvstot;
  for(i = 0; i < nfs; ++i)
    Bsys->face[i] += nvstot + nestot;
  
  for(i = 0; i < nek+1; ++i)
    Bsys->edge[nes+i] += nvstot + nfstot + nvktot ;
  for(i = 0; i < nfk+1; ++i)
    Bsys->face[nfs+i] += nvstot + nestot + nvktot + nektot;

  /* Scatter new local values back to vertex edge and face id's */
  for(E=EL->fhead;E;E=E->next){
    for(i = 0; i < E->Nverts; ++i)
      E->vert[i].gid = vmap[E->vert[i].gid];
    
    for(i = 0; i < E->Nedges; ++i)
      E->edge[i].gid = emap[E->edge[i].gid];
    
    for(i = 0; i < E->Nfaces; ++i)
      E->face[i].gid = fmap[E->face[i].gid];
  }  


  /* finally if parallel option then there is enough information in
     vmap,emap,fmap and gmap to define a global-local mapping for the
     boundary degrees of freedom */
#ifdef PARALLEL
    int j;
    Pllmap *p;
    /* set up mapping for solved points */
    p = Bsys->pll = (Pllmap*)malloc(sizeof(Pllmap));
    p->nv_solve   = gmap->nvs;
    p->nsolve     = gmap->nsolve;
    p->nglobal    = gmap->nglobal;
    p->nv_gpsolve = gmap->nv_psolve;
    p->nv_lpsolve = 0;
    
    /* set up solve map */
    p->solvemap = ivector(0,Bsys->nsolve);
    
    for(i = 0; i < gmap->nvs; ++i)
      if(vmap[i]+1){
	p->solvemap[vmap[i]] = gmap->vgid[i];
	if(gmap->vgid[i] < p->nv_gpsolve)
	  p->nv_lpsolve++;
      }
    
    for(i = p->nv_lpsolve; i < nvs; ++i)
      if(gmap->vgid[p->solvemap[i]] < p->nv_gpsolve)
	fprintf(stderr,	"Error in parallel vertex numbering: "
		"gid:%d, i:%d, nv_gpsolve: %d\n",gmap->vgid[p->solvemap[i]],
		i,p->nv_gpsolve);
    
    for(i = 0; i < gmap->nes; ++i)
      if(emap[i]+1)
	for(j = 0; j < elen[i]; ++j)
	  p->solvemap[Bsys->edge[emap[i]] + j] = gmap->egid[i] + j;
    
    for(i = 0; i < gmap->nfs; ++i)
      if(fmap[i]+1)
	for(j=0; j < flen[i]; ++j)
	  p->solvemap[Bsys->face[fmap[i]]+j] = gmap->fgid[i] + j;
    
    /* extra mappings for ddot sums in Bsolve_CG */
    p->solve = gs_init(p->solvemap,Bsys->nsolve,option("GSLEVEL"));
    
    if(Bsys->nsolve){
      p->mult = dvector(0,Bsys->nsolve-1);
      dfill(Bsys->nsolve,1.,p->mult,1);
      gs_gop(p->solve, p->mult, "+");
      dvrecp(Bsys->nsolve,p->mult,1,p->mult,1);
    }
    
    if(Bsys->nglobal - Bsys->nsolve){
      int nloc, nglo;
      
      nloc = Bsys->nsolve;
      nglo = gmap->nsolve;
      
      p->knownmap = ivector(0,Bsys->nglobal - Bsys->nsolve);
      
      for(i = gmap->nvs; i < gmap->nvg; ++i)
	if(vmap[i]+1)
	  p->knownmap[vmap[i]-nloc] = gmap->vgid[i]-nglo;
      
      for(i = gmap->nes; i < gmap->neg; ++i)
	if(emap[i]+1)
	  for(j=0; j < elen[i]; ++j)
	    p->knownmap[Bsys->edge[emap[i]]-nloc+j] = gmap->egid[i]+j -nglo;
      
      for(i = gmap->nfs; i < gmap->nfg; ++i)
	if(fmap[i]+1)
	  for(j=0; j < flen[i]; ++j)
	    p->knownmap[Bsys->face[fmap[i]]-nloc+j] = gmap->fgid[i]+j-nglo;
      
    }
    else
      p->knownmap = ivector(0,0); 
    
    /* this must be called outside of if statement since even if there
       are no points all processors must send gs_init */
    
    p->known=gs_init(p->knownmap,Bsys->nglobal-Bsys->nsolve,option("GSLEVEL"));
#endif
  
  free(vmap); free(emap);  free(elen);  free(fmap);  free(flen);
}

typedef struct vertnum {
  int    id; 
  struct vertnum *base;
  struct vertnum *link;
} Vertnum;


static void setGid(Element_List *UL){
  Element  *U = UL->fhead;
  register int i,j,k;
  int      nvg, face, eid, edgeid, faceid, vertmax;
  const    int nel = UL->nel;
  Edge     *e,*ed;
  Face     *f;
  Vertnum  *V,*vb,*v,*vc;

  Element *E;

  /* set vector of consequative numbers */
  /* set up vertex list */
  vertmax = Max_Nverts;

  V = (Vertnum *) calloc(vertmax*nel,sizeof(Vertnum));
  
  if(U->dim() == 2){
    for(E = U; E; E = E->next)
      for(i = 0; i < E->Nedges; ++i){
	for(j = 0; j < E->dim(); ++j){

	  v = V+E->id*vertmax + E->fnum(i,j);

	  if(E->edge[i].base){
	    if(E->edge[i].link){
	      eid  = E->edge[i].link->eid;
	      face = E->edge[i].link->id;
	    }
	    else{
	      eid  = E->edge[i].base->eid;
	      face = E->edge[i].base->id;
	    }

	    vb = V[eid*vertmax + UL->flist[eid]->fnum1(face,j)].base;
	    
	    if(eid < E->id){  /* connect to lower element */
	      if(!v->base) v->base = v;
	      
	      /* search through all points and assign to same base */
	      for(;vb->link;vb = vb->link);
	      if(vb->base != v->base) vb->link = v->base;
	      for(v = v->base;v; v = v->link) v->base = vb->base;
	    }
	    else if(!v->base) v->base = v;
	  }
	  else if(!v->base) v->base = v;
	}
      }
  }

  if(U->dim() == 3){
    double x, y, z, x1, y1, z1, cx, cy, cz, TOL = 1E-5;
    int vn, vn1, flag, nfv;
    Element *F;
    
    for(E = U; E; E = E->next)
      for(i = 0; i < E->Nfaces; ++i){
	for(j = 0; j < E->Nfverts(i); ++j){
	  
	  vn = E->vnum(i,j);

	  v = V+E->id*vertmax + vn;
	  
	  if(E->face[i].link){
	    eid  = E->face[i].link->eid;
	    face = E->face[i].link->id;
	    F    = UL->flist[eid];
	    
	    nfv = F->Nfverts(face);

	    x = E->vert[vn].x, y = E->vert[vn].y, z = E->vert[vn].z;
	    
	    cx = 0.0; cy = 0.0; cz = 0.0;
	    for(k=0;k<nfv;++k){
	      cx += F->vert[F->vnum(face,k)].x - E->vert[E->vnum(i,k)].x;
	      cy += F->vert[F->vnum(face,k)].y - E->vert[E->vnum(i,k)].y;
	      cz += F->vert[F->vnum(face,k)].z - E->vert[E->vnum(i,k)].z;
	    }
	    cx /= 1.*nfv;      cy /= 1.*nfv;      cz /= 1.*nfv;

	    // loop through vertices on neighbour face
	    // break out when match is made
	    flag = 1;
	    for(k=0;k < nfv;++k){
	      vn1 = F->vnum(face,k);
	      x1  = F->vert[vn1].x-cx; 
	      y1  = F->vert[vn1].y-cy;
	      z1  = F->vert[vn1].z-cz;
	      if(sqrt((x1-x)*(x1-x)+ (y1-y)*(y1-y) + (z1-z)*(z1-z)) < TOL){
		flag = 0;
		break;
	      }
	    }
	    
	    if(flag) 
	      fprintf(stderr, "Error in SetGid  Elmt:%d Face:%d Vertex:%d\n",
		      E->id, i, j);
	    
	    vb = V[eid*vertmax + vn1].base;
	    
	    /* connect to lower element */
	    
	    if(eid <= E->id){  
	      if(!v->base) v->base = v;
	      
	      /* search through all points and assign to same base */
	      if(vb){
		for(;vb->link;vb = vb->link);
		if(vb->base != v->base) vb->link = v->base;
		for(v = v->base;v; v = v->link) v->base = vb->base;
	      }
	    }
	    else if(!v->base) v->base = v;
	  }
	  else if(!v->base) v->base = v;
	}
      }
  }
  
  /* number vertices consequatively */
  for(E = U, nvg = 1; E; E = E->next)
    for(i = 0; i < E->Nverts; ++i)
      if(!V[E->id*vertmax+i].base->id) 
	V[E->id*vertmax+i].id = V[E->id*vertmax+i].base->id = nvg++;
      else                       
	V[E->id*vertmax+i].id = V[E->id*vertmax+i].base->id;
  nvg--;

  for(E = U; E; E = E->next)
    for(i = 0; i < E->Nverts; ++i)
      E->vert[i].gid = V[E->id*vertmax+i].id-1;
  
  /* set gid's to -1 */
  for(E = U; E; E = E->next){
    for(i = 0; i < E->Nedges; ++i) E->edge[i].gid = -1;
    for(i = 0; i < E->Nfaces; ++i) E->face[i].gid = -1;
  }
  
  /* at present just number edge and faces consequatively */
  faceid = 0; edgeid = 0;
  for(E = U; E; E = E->next){
    e = E->edge;
    for(i = 0; i < E->Nedges; ++i){
      if(e[i].gid==-1)
	if(e[i].base){
	  for(ed = e[i].base; ed; ed = ed->link)
	    ed->gid = edgeid;
	  ++edgeid;
	}
	else
	  e[i].gid = edgeid++;
    }
    f = E->face;
    for(i = 0; i < E->Nfaces; ++i)
      if(f[i].gid==-1)
	if(f[i].link){
	  f[i].gid       = faceid;
	  f[i].link->gid = faceid++;
	}
	else
	  f[i].gid       = faceid++;
  }
  
  free(V);
  return;
}

static int suminterior(Element *E){
  int sum=0;
  for(;E;E = E->next)
    sum += E->Nmodes - E->Nbmodes;

  return sum;
}

static void set_bmap(Element *U, Bsystem *B){
  register int i,j,k,n;
  int   l;
  const int nel = B->nel;
  int  **bmap;
  Element *E;

  /* declare memory */
  bmap = (int **) malloc(nel*sizeof(int *));
  for(E=U,l=0;E;E=E->next) l += E->Nbmodes;

  bmap[0] = ivector(0,l-1);
  for(i = 0, E=U; i < nel-1; ++i, E=E->next)
    bmap[i+1] = bmap[i] + E->Nbmodes;

  /* fill with bmaps */
  for(E=U;E;E=E->next){
    for(j = 0; j < E->Nverts; ++j)
      bmap[E->id][j] = E->vert[j].gid;
    
    for(j = 0,n = E->Nverts; j < E->Nedges; ++j,n+=l){
      l = E->edge[j].l;
      for(k = 0; k < l; ++k)
        bmap[E->id][n+k] = B->edge[E->edge[j].gid] + k;
    }

    if(E->dim() == 3)
      for(k = 0; k < E->Nfaces; ++k){
	l = E->face[k].l;

	if(E->Nfverts(k) == 3){ // triangle face
	  l = l*(l+1)/2;
	  for(i = 0; i < l; ++i)
	    bmap[E->id][n+i] = B->face[E->face[k].gid] + i;
	  n += l;
	}
	else{  // square face
	  if(E->face[k].con < 4){
	    for(i=0;i<l;++i)
	      for(j=0;j<l;++j)
		bmap[E->id][n+j+i*l] = B->face[E->face[k].gid]+j+i*l;
	  }
	  else{ // transpose modes
	    for(i=0;i<l;++i)
	      for(j=0;j<l;++j)
		bmap[E->id][n+i+j*l] = B->face[E->face[k].gid]+j+i*l;
	    //	    E->face[k].con -= 4;
	  }
	  n += l*l;
	}
      }
  }

  B->bmap = bmap;
}

void free_Mesh_Bcs(Bndry *MeshBcs){
  Bndry *Bc;

  for(Bc=MeshBcs;Bc;Bc = Bc->next)
    if(Bc->type == 'Z')
      free(Bc->bvert);

  free(Bc);
}

/* this function resets the global mesh and boundary conditions so that
   they are correctly set for the velocity */

void Reflect_Global_Velocity(Element_List *Mesh, Bndry *Meshbcs, int dir){
  register int i,j,id;
  Element *U;
  Bndry   *B;
  int flag = 0;

  /* reset solve mask */
  for(U = Mesh->fhead; U; U = U->next)
    for(i = 0; i < U->Nverts; ++i)
      U->vert[i].solve = 1;

  DO_PARALLEL{ /* set interface solve masks to 2 */
    if((Mesh->fhead->dim() == 3)&&(pllinfo.partition)){
      for(U = Mesh->fhead; U; U = U->next)
	for(i = 0; i < U->Nfaces; ++i)
	  if(U->face[i].link)
	    if(pllinfo.partition[U->id] 
	       != pllinfo.partition[U->face[i].link->eid])
	      for(j = 0; j < U->Nfverts(i); ++j){
		id = U->vnum(i,j);
		U->vert[id].solve *= 2;
		//clamp for mult calls
		U->vert[id].solve = min(U->vert[id].solve,2); 
	      }
    }
  }

  /* Loop through boundary conditions and set Dirichlet boundaries */
  for(B = Meshbcs; B; B = B->next){
    switch (B->usrtype) {
    case 'V': case 'v':
    case 'W': 
      B->type = 'V';
      U = B->elmt;
      for(i = 0; i < U->Nfverts(B->face); ++i) 
	U->vert[U->vnum(B->face,i)].solve = 0;
      break;
    case 'Z':
      U = B->elmt;
      if(B->bvert[0] == dir){
	B->type = 'V';
	for(i = 0; i < U->Nfverts(B->face); ++i) 
	  U->vert[U->vnum(B->face,i)].solve = 0;
      }
      else
	B->type = 'F';
      break;
    }
  }
}




