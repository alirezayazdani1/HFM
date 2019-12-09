#include "mvector.H"
#include <string.h>
#include <stdio.h>

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
typedef Mvector<Vertx*> Vvector;


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

class PDR_Model{
 private:
 public:
  int uid_cnt;
  //need pointers here!
  Ivector SC_list;//Particle current bioState list
  Ivector SN_list;//Particle new bioState list
  Dvector ST_list;//When changed to current state
  PBvector PB_list;//list of bonds structures
  
  PDR_Model(){};
  ~PDR_Model(){};
  void Add_Particle(const int biostate_){
    SC_list += biostate_;//Particle current bioState list
    SN_list += biostate_;//Particle new bioState list
    ST_list += 0.0;//When changed to current state
    Particle_Bonds* pb = new Particle_Bonds();
    PB_list += pb;;//Part-Part-Bond list
    return;
  };
  int Get_Particle_State(const int pid){return SC_list[pid];};
  void Write_Particle_FieldFile(FILE* fout,int step,double time){
    //Dump particles
    int pid,bid;
    fprintf(fout,"PDR\n");
    fprintf(fout,"%ld  %ld %lf %ld\n",SC_list.size(),uid_cnt,time,step);
    for(pid=0;SC_list.size();++pid){
      fprintf(fout,"%ld %ld %lf %ld %ld\n",
	      SC_list[pid],SN_list[pid],ST_list[pid],
	      PB_list[pid]->PBP_list.size(),
	      PB_list[pid]->PBW_list.size());
      for(bid=0;bid<PB_list[pid]->PBP_list.size();++bid){
	fprintf(fout,"%ld\n",PB_list[pid]->PBP_list[bid]);
      }
      for(bid=0;bid<PB_list[pid]->PBW_list.size();++bid){
	fprintf(fout,"%lf %lf %lf\n",
		PB_list[pid]->PBW_list[bid]->x,
		PB_list[pid]->PBW_list[bid]->y,
		PB_list[pid]->PBW_list[bid]->z);
      }
    }
  };
  
  void Read_Particle_FieldFile(FILE* fin,int &step,double &time){
    char buf[BUFSIZ];
    int pid,bid,pSC,pSN,pPBP,pPBW,pb,pnum,puid_cnt;
    double pST,px,py,pz;
    fgets(buf, BUFSIZ, fin);
    if (strcmp(buf,"PDR") == 0){
      fprintf(stderr,"Wrong Bio Model!\n");
      exit(-1);
    }
    fgets(buf, BUFSIZ, fin);
    sscanf(buf,"%ld %ld %*lf %*ld\n",&pnum,&puid_cnt);
    fprintf(stderr,"%ld %ld \n",pnum,puid_cnt);
    
    for(pid=0;pid<pnum;++pid){
      fgets(buf, BUFSIZ, fin);
      sscanf(buf,"%ld %ld %lf %ld %ld ",
	     &pSC,&pSN,&pST,&pPBP,&pPBW);
      fprintf(stderr,"%ld %ld %lf %ld %ld ",
	      pSC,pSN,pST,pPBP,pPBW);
      Add_Particle(pSC);
      SC_list[pid]=pSC; SN_list[pid]=pSN;ST_list[pid]=pST;
      for(bid=0;bid<pPBP;++bid){
	fgets(buf, BUFSIZ, fin);
	sscanf(buf,"%ld ",&pb);
	fprintf(stderr,"%ld ",pb);
	PB_list[pid]->PBP_list.add_once(pb);
      }
      for(bid=0;bid<pPBW;++bid){
	fgets(buf, BUFSIZ, fin);
	sscanf(buf,"%lf %lf %lf ",&px,&py,&pz);
	fprintf(stderr,"%lf %lf %lf ",px,py,pz);
	Vertx* v = new Vertx(px,py,pz);
	(PB_list[pid]->PBW_list) += v;
      }
      fprintf(stderr,"\n");
    }
  };
};//class PDR_Model


