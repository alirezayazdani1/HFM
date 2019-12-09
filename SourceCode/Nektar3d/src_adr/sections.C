#include <stdio.h>
#include <string.h>
#include "nektar.h"
#include <stdlib.h>
#include <math.h>
#include <veclib.h>
#include <ctype.h>
#include <time.h>




//Extraction of sectional information
//First count the number of surface,'W', elements

int cnt_srf_elements(Domain *omega){ 
  
  FILE  *fp_out=fopen("n_srf","w");  
  Bndry *B;
  int srf_nel = 0;

  for(B=omega->Ubc;B;B=B->next){
    if(B->type == 'W') { 
      ++srf_nel;
    }
  }
  
  
  fprintf(fp_out,"surf edges on this partition = %d\n",srf_nel);
  fclose(fp_out); 
  
  return srf_nel;
}



Element_List *setup_surflist(Domain *omega, Bndry *Ubc, char type){
  Bndry *Bc;
  Element **selmt;
  int S_nbcs, i,j,k;  
  Coord X;
  iparam_set("FAMILIES", 0);

  Element_List *U = omega->U;

  X.x = dvector(0, 4);
  X.y = dvector(0, 4);
  
  for(S_nbcs = 0, Bc = Ubc ; Bc ; Bc=Bc->next)
    if(Bc->type == type){
      ++S_nbcs;
#ifdef TRI_SURF
      if(Bc->elmt->Nfverts(Bc->face) != 3){
	fprintf(stderr,"Quad face attached to surface: Check setup_trilist\n");
	exit(-1);
      }
#endif
    }

  if(!S_nbcs)
    return (Element_List *) NULL;    

  selmt = (Element **) calloc (S_nbcs, sizeof(Element *));

  i = 0;
  for(Bc = Ubc;Bc;Bc=Bc->next)
    if(Bc->type == type){  
      X.x[0] = -1.0;    X.y[0] = -1.0;
      X.x[1] =  1.0;    X.y[1] = -1.0;
      if(Bc->elmt->Nfverts(Bc->face) == 3){
	X.x[2] = -1.0;    X.y[2] = 1.0;
	
	int faceid = (Bc->elmt->identify() == Nek_Tet)?  0: 1;
	
	selmt[i++] = (Element*) new Tri(i, type,Bc->elmt->lmax, 
				      Bc->elmt->get_face_q1(faceid), 
				      Bc->elmt->get_face_q2(faceid), 0, &X);
      }
      else{
	X.x[2] =  1.0;	X.y[2] = 1.0;
	X.x[3] = -1.0;	X.y[3] = 1.0;
	selmt[i++] = (Element*) new Quad(i, type,Bc->elmt->lmax, 
				       Bc->elmt->get_face_q1(0), 
				       Bc->elmt->get_face_q2(0), 0, &X);
      }
    }
  
  for(i = 0; i < S_nbcs-1; ++i) 
    selmt[i]->next = selmt[i+1];
  selmt[i]->next = (Element*) NULL;
  
   Element_List *surf_list = (Element_List*) new Element_List(selmt, S_nbcs);
   surf_list->Cat_mem();
   
   i = 0;
   for(Bc = Ubc;Bc;Bc=Bc->next)
     if(Bc->type == type){ 
       selmt[i]->id = i;
       selmt[i]->set_geofac();
       selmt[i]->set_curved_elmt((Element_List*)NULL);
       for(j=0;j<Bc->elmt->Nfverts(Bc->face);++j){
	 selmt[i]->vert[j].gid = 
	   Bc->elmt->vert[Bc->elmt->vnum(Bc->face,j)].gid;
	 selmt[i]->edge[j].gid = 
	   Bc->elmt->edge[Bc->elmt->ednum(Bc->face,j)].gid;
	 selmt[i]->edge[j].con = 
	   Bc->elmt->edge[Bc->elmt->ednum(Bc->face,j)].con;
      }
       ++i;
     }
   
  
   free(X.x);  free(X.y);
   return surf_list;
}   




    








