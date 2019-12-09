/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $Source: 
 * $Revision:
 * $Date: 
 * $Author: piv@cfm.brown.edu
 * $State: 
 *---------------------------------------------------------------------------*/
#ifndef _FCM_H_
#define _FCM_H_

#include "nektar.h"
#include "nekstruct.h"
#include "mvector.h"
#include <string.h>

typedef enum{SPHERE, TRACER}FCM_PARTICLE;

class Vertx{
 public:
  double x,y,z;
  Vertx(){};
  Vertx(const double x_,const double y_, const double z_){
    x=x_;y=y_;z=z_;
  }	
  ~Vertx(){};  

  Vertx(const Vertx& l){
    fprintf(stderr,"Vop1 is being called\n");
    x=l.x;y=l.y;z=l.z;
  };
  
  
  Vertx& operator =(const Vertx& l){
  fprintf(stderr,"Vop2 is being called\n");
    // assignment operator:
    if (&l == this)
      return *this;
    (*this).x=l.x;
    (*this).y=l.y;
    (*this).z=l.z;
    return *this;
  }
  bool operator ==(const Vertx &c) const {
    if (!((*this).x == c.x))
      return 0;
    if (!((*this).y == c.y))
      return 0;
    if (!((*this).z == c.z))
      return 0;
    return 1;
  }
};

typedef Mvector<int> Ivector;
typedef Mvector<double> Dvector;


//---******--- PLATELET MODEL ---******---//

class Bio_Model{
 private:
 public:
  FCM_Domain* FCM_Omega;
  Bio_Model(FCM_Domain* FCM_Omega_){FCM_Omega=FCM_Omega_;};
  ~Bio_Model(){};
  virtual void Bio_Force(const int pid, double &f0, double &f1, double &f2)
    //Add bio force to the particle force monopole
    {fprintf(stderr,"Virt func is being called!!\n"); return;};
  virtual void Add_Particle(const int biostate)
    //Add particle to the module
    {fprintf(stderr,"Virt func is being called!!\n"); return;};
  virtual void Remove_Particle(const int pid)
    //Remove particle to the module
    {fprintf(stderr,"Virt func is being called!!\n"); return;};
#if defined (ADR) && (CCHF)
  virtual void Set_Particle_Biostate()
    //Specific to activation mechanism, called each time step
    {fprintf(stderr,"Virt func is being called!!\n"); return;};
#endif
  virtual void Bio_Work()
    //Specific to Bio Model work, called each time step
    {fprintf(stderr,"Virt func is being called!!\n"); return;};
  virtual void Particle_Particle_Interaction(const int pid1, const int pid2)
    //Called if 2 particles interact
    {fprintf(stderr,"Virt func is being called!!\n"); return;};
  virtual void Particle_Wall_Interaction(const int pid, const int bid)
    //Called if particle wall interact
    {fprintf(stderr,"Virt func is being called!!\n"); return;};
  virtual int Get_Particle_State(const int pid)
    //Return particle current state
    {fprintf(stderr,"Virt func is being called!!\n"); return -999;};
  virtual int Free_Particle(const int pid)
    //Return particle current state
    {fprintf(stderr,"Virt func is being called!!\n"); return -999;};
  virtual void Write_Particle_FieldFile(FILE* fout,int step,double time)
    {fprintf(stderr,"Virt func is being called!!\n"); return;};
  virtual void Read_Particle_FieldFile(FILE* fin,int &step,double &time)
    {fprintf(stderr,"Virt func is being called!!\n"); return;};
 virtual int FBAR_Particle(const int pid)
    //Return particle current state
    {fprintf(stderr,"Virt func is being called!!\n"); return -999;};
  };//class Bio_Model

class Dummy_Model : public Bio_Model{
 private:
 public:
  Dummy_Model(FCM_Domain* FCM_Omega_) : Bio_Model(FCM_Omega_){};
  ~Dummy_Model(){};
  void Bio_Force(const int pid, double &f0, double &f1, double &f2){return;};
  void Add_Particle(const int biostate){return;};
  void Remove_Particle(const int pid){return;};
  void Bio_Work(){return;};
  void Particle_Particle_Interaction(const int pid1, const int pid2){return;};
	void Particle_Wall_Interaction(const int pid, const int bid){return;};
  int Get_Particle_State(const int pid){return -999;};
  int Free_Particle(const int pid){return -999;};
  void Write_Particle_FieldFile(FILE* fout,int step,double time){
    fprintf(fout,"Dummy\n"); return;
  };
  void Read_Particle_FieldFile(FILE* fin,int &step,double &time){
    return;
  };
 int FBAR_Particle(const int pid)
    //Return particle current state
    {return -999;};

};//class Dummy_Model


typedef Mvector<Vertx*> Vvector;

//--- PDR MODEL ---//

typedef enum{PASSIVE,PREACTIVE,ACTIVE,ADHERED}PDR_BSTATES;

class Particle_Bonds{
 private:
 public:
  Ivector PBP_list;
  Vvector PBW_list;
  Particle_Bonds(){};
  ~Particle_Bonds(){
    fprintf(stderr,"~Particle_Bonds()\n");
    //TODO -- fix bonds when removing particles!!
  };
  Particle_Bonds(const Particle_Bonds& l){
    fprintf(stderr,"Op1 is being called\n");
    PBP_list = l.PBP_list;
    PBW_list = l.PBW_list;
  };
  
  
  Particle_Bonds& operator =(const Particle_Bonds& l){
    // assignment operator:
    fprintf(stderr,"Assgn op = is being called\n");
    if (&l == this)
      return *this;
    (*this).PBP_list = l.PBP_list;
    (*this).PBW_list = l.PBW_list;
    return *this;
  };
  bool operator ==(const Particle_Bonds &c) const {
    fprintf(stderr,"Bool == is being called!\n");
    if (!((*this).PBP_list == c.PBP_list))
      return 0;
    if (!((*this).PBW_list == c.PBW_list))
      return 0;
    return 1;
  };
};
typedef Mvector<Particle_Bonds*> PBvector;

class PDR_Model : public Bio_Model{
 private:
 public:
  double PPBFcoef,PWBFcoef;//Part-part and -wall bond coeff
  double BioATime,BioPATime_min,BioPATime_max;
  //need pointers here!
  Ivector SC_list;//Particle current bioState list
  Ivector SN_list;//Particle new bioState list
#ifdef PARALLEL
  Ivector WK_list;//temp store for parallel runs
#endif //PARALLEL
  Dvector ST_list;//When changed to current state
  Dvector ST_Act_list;	// Alireza: when changed to activation state
  PBvector PB_list;//list of bonds structures
  
  PDR_Model(FCM_Domain* FCM_Omega_);
  ~PDR_Model(){};
  void Bio_Force(const int pid, double &f0, double &f1, double &f2);
  void Particle_Particle_Bond_Force(const int pid, double &f0, double &f1, double &f2);
  void Particle_Wall_Bond_Force(const int pid, double &f0, double &f1, double &f2);
  void Add_Particle(const int biostate);
  void Remove_Particle(const int pid);
#if defined (ADR) && (CCHF)
	void Set_Particle_Biostate();
#endif
  void Bio_Work();
  void Particle_Particle_Interaction(const int pid1, const int pid2);
	void Particle_Wall_Interaction(const int pid, const int bid);
  int Get_Particle_State(const int pid){return SC_list[pid];};
  int Free_Particle(const int pid);
  void Write_Particle_FieldFile(FILE* fout,int step,double time);
  void Read_Particle_FieldFile(FILE* fin,int &step,double &time);
  int FBAR_Particle(const int pid);
};//class PDR_Model

//---******--- BACKGROUND MESH ---******---//

class  BGElement{
 private:
 public:
  Ivector EID_list;//List of Nektar elements assos with this
  Ivector PID_list;//List of particles assos with this
  BGElement(){};
  ~BGElement(){};
};


class BGMesh{
 private:
  FCM_Domain* FCM_Omega;
  double xmin,xmax,ymin,ymax,zmin,zmax;
  double elsize,invelsize;//size of the elements in the BGmesh
  int xnel,ynel,znel;
 public:
  BGElement**** mesh;
  BGMesh(FCM_Domain *FCM_Omega_){
    FCM_Omega=FCM_Omega_;
    Set_Elsize();
    Set_Bbox();
    Set_Mesh();
    Distribute_Elements();
    //Distribute_Particles();
  }
  ~BGMesh(){};//need to fill in
  void Set_Elsize();
  void Set_Bbox();
  void Set_Mesh();
  void Distribute_Elements();
  void Add_Element(int i,int j,int k,int eid);
  void Distribute_Particles();
  inline int Get_Location(const double x,const double y,const double z, int &i, int &j, int &k){
    //computes i,j,k, given part coords; returns 0 if not in mesh, 1 otherwise
    i=-(int)floor((xmin-x)*invelsize);
    j=-(int)floor((ymin-y)*invelsize);
    k=-(int)floor((zmin-z)*invelsize);
    if((i<1)||(i>xnel)) return 0;
    if((j<1)||(j>ynel)) return 0;
    if((k<1)||(k>znel)) return 0;
    --i;--j;--k;
    if(mesh[i][j][k]==NULL) return 0;
    return 1;
  };
  inline int Valid_Index(const int i,const int j,const int k){
    if((i<0)||(i>=xnel)) return 0;
    if((j<0)||(j>=ynel)) return 0;
    if((k<0)||(k>=znel)) return 0;
    return 1;
  };
};

//---******--- PARTICLES ---******---//

class Particle{
 protected:
  int id;
  int uid;//unique id (need for visualization)
  FCM_PARTICLE type;
 public:
  double lf0,lf1,lf2;
  int inertia_flag;
  double mass,volume; //mass and volume of the particle
  double radius;//maximum radius
  double venvelope,fenvelope,msigma,fsigma;
  Particle(double mass_,double radius_);
  ~Particle(){//need to fill in destructor
    fprintf(stderr, "Particle destructor is being called!!\n");
  }
  void Set_id(const int id_){id=id_; return;};
  inline int Get_id(){return id;};
  void Set_uid(const int uid_){uid=uid_; return;};
  inline int Get_uid(){return uid;};

  inline FCM_PARTICLE Get_type(){return type;};
  virtual void Velocity_Element_Integral(FCM_Domain* FCM_Omega ,const int eid,
					 double &u_,double &v_,double &w_,double &s_);
#if defined (ADR) && (CCHF)
  virtual void Element_Integral(FCM_Domain* FCM_Omega ,const int eid,
					 double &u_,double &v_,double &w_,double &s_,double &actval);
  virtual void Add_Concen_Source(FCM_Domain* FCM_Omega, const int eid)
    {fprintf(stderr,"Virt func is being called!!\n"); return;};
#endif
  virtual void Add_Monopole(FCM_Domain* FCM_Omega, const int eid)
    {fprintf(stderr,"Virt func is being called!!\n"); return;};
  virtual void Add_Particle_Number(FCM_Domain* FCM_Omega, const int eid)
    {fprintf(stderr,"Virt func is being called!!\n"); return;};
  virtual void External_Force(FCM_Domain* FCM_Omega,
			      double &f0, double &f1, double &f2)
    {fprintf(stderr,"Virt func is being called!!\n"); return;};
  void Inertia_Force(FCM_Domain* FCM_Omega,
			      double &f0, double &f1, double &f2);
  virtual void Particle_Force(FCM_Domain* FCM_Omega,
			      double &f0, double &f1, double &f2){return;};  
  int Particle_Particle_Force(FCM_Domain* FCM_Omega, int pid2,
			      double &tx, double &ty, double &tz,
			      double &dx, double &dy, double &dz,
			      double &fm){return -1;};
  int Particle_Particle_Dissipation(FCM_Domain* FCM_Omega, int pid2,
                              double &dx, double &dy, double &dz,
                              double &fm){return -1;};
};//class Particle

class Sphere : public Particle{
 private:
 public:
  Sphere(double mass_,double radius_) :
    Particle(mass_,radius_){type = SPHERE;};
  ~Sphere(){};//fill in
  void Add_Monopole(FCM_Domain* FCM_Omega, const int eid);
  void Add_Particle_Number(FCM_Domain* FCM_Omega, const int eid);
#if defined (ADR) && (CCHF)
	void Add_Concen_Source(FCM_Domain* FCM_Omega, const int eid);
#endif
  void External_Force(FCM_Domain* FCM_Omega,
		      double &f0, double &f1, double &f2);
  void Particle_Force(FCM_Domain* FCM_Omega,double &f0, double &f1, double &f2);
  int Particle_Particle_Force(FCM_Domain* FCM_Omega, int pid2,
			      double &tx, double &ty, double &tz,
			      double &dx, double &dy, double &dz,
			      double &fm);
  int Particle_Particle_Dissipation(FCM_Domain* FCM_Omega, int pid2,
                              double &dx, double &dy, double &dz,
                              double &fm);
};//class Sphere

class Tracer : public Particle{
 private:
 public:
  Tracer(double mass_,double radius_) :
    Particle(mass_,radius_){type = TRACER;};
  ~Tracer(){};//fill in 
  void Add_Monopole(FCM_Domain* FCM_Omega, const int eid){return;};
  void Add_Particle_Number(FCM_Domain* FCM_Omega, const int eid){return;};
#if defined (ADR) && (CCHF)
	void Add_Concen_Source(FCM_Domain* FCM_Omega, const int eid){return;};
#endif
  void External_Force(FCM_Domain* FCM_Omega,
		      double &f0, double &f1, double &f2){return;};
  void Particle_Force(FCM_Domain* FCM_Omega,
		      double &f0, double &f1, double &f2){return;};
  int Particle_Particle_Force(FCM_Domain* FCM_Omega, int pid2,
			      double &tx, double &ty, double &tz,
			      double &dx, double &dy, double &dz,
			      double &fm){return -1;};
  int Particle_Particle_Dissipation(FCM_Domain* FCM_Omega, int pid2,
                              double &dx, double &dy, double &dz,
                              double &fm){return -1;};
  void Wall_Force(FCM_Domain* FCM_Omega, 
		  double &f0, double &f1, double &f2){return;};
};//class Tracer




//---******--- NEKTAR ELEMENT GEOM INFO ---******---//

class Element_Geometry_Info{
 private:
 public:
  double *cx,*cy,*cz;//coordinates of element centers
  double *rad;//the radius of enc. sphere
  Element_Geometry_Info(Element_List *U){
    Element* E;
    int i,nel,eid;
    double r2;
    nel=U->nel;
    cx=dvector(0,nel-1); dzero(nel,cx,1);
    cy=dvector(0,nel-1); dzero(nel,cy,1);
    cz=dvector(0,nel-1); dzero(nel,cz,1);
    rad=dvector(0,nel-1); dzero(nel,rad,1);
    //fiill in centers of elements
    for(eid=0;eid<nel;++eid){
      E=U->flist[eid];
      for(i=0;i<E->Nverts;++i){
	cx[eid]+=E->vert[i].x;
	cy[eid]+=E->vert[i].y;
	cz[eid]+=E->vert[i].z;  
      }
      cx[eid]=cx[eid]/(double)E->Nverts;
      cy[eid]=cy[eid]/(double)E->Nverts;
      cz[eid]=cz[eid]/(double)E->Nverts;
    }
    
    //find the radius of env.sphere
    for(eid=0;eid<nel;++eid){
      E=U->flist[eid];
      r2=0.0;
      for(i=0;i<E->Nverts;++i){
	r2=(cx[eid]-E->vert[i].x)*(cx[eid]-E->vert[i].x)
	  +(cy[eid]-E->vert[i].y)*(cy[eid]-E->vert[i].y)
	  +(cz[eid]-E->vert[i].z)*(cz[eid]-E->vert[i].z);
	rad[eid]=max(rad[eid],r2);
      }
      rad[eid]=sqrt(rad[eid]);
    }
  };
  ~Element_Geometry_Info(){};//need to fill in destructor
};//class Element_Geometry_Info

//---******--- FCM BOUNDARY ---******---//

class Boundary{
 private:
 public:
  int curved;//flag for curvature
  int npts;//number of points
  int nverts;//number of verts
  double cx,cy,cz;//coordinates of center
  double rad;//radius of env.sphere
  char type;//type of BC
  int eid, fid;
  Vertx vert[4];
  Dorp nx,ny,nz;//components of surface normal
  Coord X;//coordinates of surface points
  Boundary(int npts_, int curved_){
    npts=npts_;    curved=curved_;
    X.x=dvector(0,QGmax*QGmax*QGmax-1); 
    X.y=dvector(0,QGmax*QGmax*QGmax-1); 
    X.z=dvector(0,QGmax*QGmax*QGmax-1);
    if(curved){
      nx.p=dvector(0,npts-1); //need QGmax??
      ny.p=dvector(0,npts-1); nz.p=dvector(0,npts-1);
    }
    else{
      nx.p=NULL; ny.p=NULL; nz.p=NULL;
    }
  }
  ~Boundary(){};
};//class Boundary

typedef Mvector<Boundary*> Boundary_List;

class FCM_Boundary_Info{
 private:
 public:
  //periodic data
  double xpermax,xpermin,ypermax,ypermin,zpermax,zpermin;
  int xperflag,yperflag,zperflag;
  int outflowflag;

  char inflowtype;
  int inflowflag;
  int* inflow_list;
  double* inflow_weight;

  Boundary_List B_list;
  FCM_Boundary_Info(Element_List* Mesh, Bndry *Meshbc);
  ~FCM_Boundary_Info(){};
  void Check_Periodic_Boundary(FCM_Domain* FCM_Omega);
  void Check_Outflow_Boundary(FCM_Domain* FCM_Omega);
  void Set_Random_Inflow();
  void Add_Random_Particles(FCM_Domain* FCM_Omega);
  void Add_Random_Particle(FCM_Domain* FCM_Omega);
  int Get_Random_Inflow_Element();
  void Get_Random_Inflow_Coordinates(int eid,double &xp,double &yp,double &zp);
};//class FCM_Boundary_Info

typedef Mvector<Particle*> Particle_List;


//---****** --- FCM_Domain --- ******---//

class FCM_Domain{
 private:
 public:
  double dt;//timestep
  double time;//current nektar time
  int nel;//#of elements in Omega
  Element_List *U,*V,*W,*P;//pointers to nektar Elemet_lists
  //Temporary storage
  double *VEI_M,*VEI_b;
  Coord VEI_X;
  double PPFcoef,PWFcoef;//Particle-particle and part-wall force coeff

  double EFX,EFY,EFZ;//external force data
  Particle_List P_list; //particle list 
  Dvector X_list; Dvector Y_list; Dvector Z_list; //particle coordinate vector
  //SSTI->
  Dvector X1_list; Dvector Y1_list; Dvector Z1_list; //particle coordinate vector
  Dvector X2_list; Dvector Y2_list; Dvector Z2_list; //particle coordinate vector
  Dvector X3_list; Dvector Y3_list; Dvector Z3_list; //particle coordinate vector
  //SSTI=
  Dvector U_list; Dvector V_list; Dvector W_list; //particle velocity vector
  Dvector U1_list; Dvector V1_list; Dvector W1_list; //previous velocity vector
  Dvector U2_list; Dvector V2_list; Dvector W2_list; //previous velocity vector
	Dvector shear_list; // particle shear rate vector
#if defined (ADR) && (CCHF)
	Dvector C_list; // paricle concentration vector
#endif
#ifdef PARALLEL
  Dvector WK_list;//temp store for parallel runs
#endif //PARALLEL

  Domain *Omega;//Nektar domain
  int numpart; //total number of particles
  int numfixedpart; // HeLi total number of fixed particles
  Element_Geometry_Info* Element_GI;
  FCM_Boundary_Info* FCM_BI; //fcm_boundary infomation

  BGMesh* BG_Mesh;
  Ivector BG_list;

  //Bio_Model* BioM;
  PDR_Model* BioM;	//Alireza: need to make sure this is correct

  int uid_cnt; //counter for unique part ids(need for viz)
  
  
  FCM_Domain(Domain* Omega_);
  ~FCM_Domain(){};
  void Set_Current_Time(const double time_){time=time_;}; 
  double Get_Current_Time(){return time;}; 
  void Update_Particle_Coordinates();
  void Compute_Particle_Velocities(int step);
  void Add_FCM_Force();//update nektar forcing function
  int Particle_Particle_Force(int pid1, int pid2,
			      double &tx, double &ty, double &tz,
			      double &dx, double &dy, double &dz,
			      double &fm);
  void Particle_Wall_Force(const int pid,double &f0, double &f1, double &f2);

  void Add_Particle(const char ptype_, const int biostate_,
		    const double mass_, const double radius_,
		    double x_, double y_, double z_,
		    double u_, double v_, double w_,
                    const int bio_flag);
  void Remove_Particle(const int pid);
  void Remove_Particle(){};
  void Read_Particle_Field();
  void Read_Boundary_Conditions();
  void Write_Particle_FieldFile(FILE* out,int step,double time);
  void Read_Particle_FieldFile(FILE* out,int &step,double &time);
  double Total_Force_Function(const double x_, const double y_, const double z_);
  void GMultiplier(const int qt, const double sigma, double* r2, double* rez);
  inline void Update_Particle_Ids(){
    //Better not to use mvector.myid method for each particle:O(n^2) total time
    //Call this once after adding and removing particles!!
    for(int pid=0;pid<P_list.size();++pid) P_list[pid]->Set_id(pid);
  };
  void Add_Particle_Array(char ptype, int a, double m, double r,
			  double Lx, double Ly, double Lz,
			  double Dx, double Dy, double Dz,
			  int nx, int ny, int nz);
  int Empty_Region(double x, double y, double z, double r);
  void Set_uid_cnt(const int uid_cnt_){uid_cnt=uid_cnt_;};
  int Get_uid_cnt(){return uid_cnt;};
};//class FCM_Domain

//----******--- UTILITY FUNCTIONS ---******---//

inline double dist(const double x0, const double x1, const double x2,
		   const double y0, const double y1, const double y2){
  return sqrt((x0-y0)*(x0-y0)+(x1-y1)*(x1-y1)+(x2-y2)*(x2-y2));
};

inline double dist2(const double x0, const double x1, const double x2,
		   const double y0, const double y1, const double y2){
  return ((x0-y0)*(x0-y0)+(x1-y1)*(x1-y1)+(x2-y2)*(x2-y2));
};


static char *findSection (char *name, char *buf, FILE *fp){
  char *p;
  while (p = fgets (buf, BUFSIZ, fp))
    if (strstr (p, name))
      break;
  return p;
};

inline double randdbl(const double rmin, const double rmax){
#if 0 
  //return random double from rmin to rmax
  //stupid, NEED to replace!
  int rint = rand();
  double dr;
  if(rint<0) rint=-rint;
  dr=(rint%10000)/10000.;
  return rmin+(rmax-rmin)*dr;
#else
  return rmin+(rmax-rmin)*drand48();
#endif 
};

inline int randint(int maxint){
#if 0
  //return random integer from 0 to maxint-1
  int r = rand();
  if(r<0) r=-r;
  return r%maxint;
#else
  return (int)(maxint*drand48());
#endif
}


#endif
  
