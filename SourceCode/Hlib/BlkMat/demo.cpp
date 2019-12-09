#include "BlkMat.h"

/*
  g++ -I./ -I../../include -I./include -g -o demo demo.cpp -L./ -lblkmat -lblas -lg2c -llapack
*/

using namespace NekBlkMat;

main(){
  int i;
  BlkMat *A,*B,*C;
  double *mat;

  /* declare matrix
      | 1   2  -   - | 
  A = | 3   4  -   - |   using submatrix mat = | 1  2 |
      | -   -  1   2 |                         | 3  4 |
      | -   -  3   4 | 
  */

  mat = new double [4];
  mat[0] = 1; mat[1] = 2; mat[2] = 3; mat[3] = 4;
  
  cout << "A: " << endl;
  A = new BlkMat(2,2);
  A->GenBlk(0,0,2,2,mat);
  A->GenBlk(0,1,2,2,mat);
  A->GenBlk(1,1,2,2,mat);
  A->PrintBlks();  

  cout << endl << "B: " << endl;
  B = new BlkMat(2,2);
  B->GenBlk(0,0,2,2,mat);
  B->GenBlk(1,1,2,2,mat);
  B->PrintBlks();  

  cout << endl << "C=A+B: " << endl;
  C = new BlkMat(2,2);
  C->add(*A,*B);
  C->PrintBlks();


  cout << endl << "C=A-B: " << endl;
  C->sub(*A,*B);
  C->PrintBlks();

  cout << endl << "C=A*B: " << endl;
  C->MxM(*A,*B);
  C->PrintBlks();


  double *y = new(double)[4];
  double *v = new(double)[4];
  vmath::fill(4,1.0,v,1);
  vmath::zero(4,y,1);

  cout << endl << "y = A*v: "<< endl;
  C->Mxvpy(v,y);
  for(i = 0; i < 4; ++i)
    cout << y[i] << " ";
  cout << endl;

  vmath::zero(4,y,1);
  cout << endl << "y = A^T*v: "<< endl;
  //C->Mtxvpy(v,y);
  C->geMxv(ColMajor,1,v,1,y);
  for(i = 0; i < 4; ++i)
    cout << y[i] << " ";
  cout << endl;

  delete A;
  delete B;
  delete C;
  delete[] y;
  delete[] v;
  delete[] mat;

  return  0;
}
