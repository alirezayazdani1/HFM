#ifndef H_ISOUTILS
#define H_ISOUTILS

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <list>
#include <string>
#include <algorithm>

#include <veclib.h>
#include <nektar.h>

void set_SQ_PNT_TOL(double val);

class Iso{
private:
  int m_condensed;

  int   m_ntris;           // number of triangles
  int   m_nvert;           // number of vertices
  int   **m_vid;           // Vertex ids of triangles starting from 1
  double *m_x;
  double *m_y;
  double *m_z;
  int     *m_vidregion;     //region to which vertex exists
  int     m_nregions;
  int     m_eid;           // elmt id from which iso is generated (if elmt based)

public:

  void    condense();
  void    write(std::ostream &out, double val, int mincells);
  void    readzone(std::ifstream &in, double &val);
  void    globalcondense(int n_iso, Iso **iso, Element_List *);
  void    separate_regions(void);
  void    smooth(int n_iter, double lambda, double mu);

  void set_ntris( int n){
    m_ntris = n;
  }
  int  get_ntris(){
    return m_ntris ;
  }

  void set_x(int loc, double val){
    m_x[loc] = val;
  }

  void set_y(int loc, double val){
    m_y[loc] = val;
  }

  void set_z(int loc, double val){
    m_z[loc] = val;
  }

  void set_eid(int eid){
    m_eid = eid;
  }

  int get_eid(){
    return m_eid;
  }

  double get_x(int loc){
    return m_x[loc];
  }

  double get_y(int loc){
    return m_y[loc];
  }

  double get_z(int loc){
    return m_z[loc];
  }

  double *get_x(){
    return m_x;
  }

  double *get_y(){
    return m_y;
  }

  double *get_z(){
    return m_x;
  }

  void xyz_realloc(int size){
    if(m_nvert > 0){
      m_x = (double *)std::realloc(m_x,size*sizeof(double));
      m_y = (double *)std::realloc(m_y,size*sizeof(double));
      m_z = (double *)std::realloc(m_z,size*sizeof(double));
    }
    else{
      m_x = (double *)std::malloc(size*sizeof(double));
      m_y = (double *)std::malloc(size*sizeof(double));
      m_z = (double *)std::malloc(size*sizeof(double));
    }
    m_nvert = size;
  }

  Iso()
  {
    m_condensed = 0;
    m_nregions  = 0;
    m_nvert     = 0; 
    m_eid       = 0;
    m_x = (double *) NULL;
    m_y = (double *) NULL;
    m_z = (double *) NULL;
    m_vid = (int **) NULL;
    m_vidregion = (int *) NULL;
  };

  ~Iso(){
    if(m_x) free(m_x);
    if(m_y) free(m_y);
    if(m_z) free(m_z);
    if(m_vid) free_imatrix(m_vid,0,0);
    if(m_vidregion) free(m_vidregion);
  }
};


class Vertex {
  friend class Iso;
 private:
  int m_id;
  double m_x,m_y,m_z;

 public:

  Vertex (){
    m_id = -1; 
    m_x = m_y = m_z = -99999;
  }

  ~Vertex(){};
  void writevert ();

  friend bool operator == (const Vertex& x, const Vertex& y);
  friend bool operator != (const Vertex& x, const Vertex& y);

};
#endif
