#include "BlkMat.h"

/*
  g++ -I./ -I../../include -I./include -g -o SCdemo SCdemo.cpp -L./ -lblkmat -lblas -llapack -lg2c
*/

using namespace NekBlkMat;

  /* declare matrix

      | 1   2  -   - | 
  A = | 3   4  -   - |                      using submatrix mat = | 1  2 |
      | -   -  1   2 |                                            | 3  4 |
      | -   -  3   4 | 

      | 1   2   3   -   -   -   -   -   - | 
  B = | 4   5   6   -   -   -   -   -   - | using submatrix mat = | 1  2  3 |
      | 1   2   3   1   2   3   1   2   3 |                       | 4  5  6 | 
      | 4   5   6   4   5   6   4   5   6 | 

      | 1   2   3   -   -   -   -   -   - | 
      | 4   5   6   -   -   -   -   -   - | 
  C = | 7   0   0   -   -   -   -   -   - | using submatrix mat = | 1  2  3 |
      | -   -   -   1   2   3   -   -   - |                       | 4  5  6 | 
      | -   -   -   4   5   6   -   -   - |                       | 7  0  0 | 
      | -   -   -   7   0   0   -   -   - | 
      | -   -   -   -   -   -   1   2   3 | 
      | -   -   -   -   -   -   4   5   6 | 
      | -   -   -   -   -   -   7   0   0 | 

      | 1   2   1   2 | 
  D = | 3   4   3   4 |                     using submatrix mat = | 1  2 |
      | 5   6   5   6 |                                           | 3  4 | 
      | -   -   1   2 |                                           | 5  6 | 
      | -   -   3   4 | 
      | -   -   5   6 | 
      | -   -   1   2 |      
      | -   -   3   4 | 
      | -   -   5   6 | 

      Calculate DC = A - B*C*D

  */

main(){
  BlkMat *A,*B,*C,*D,*T,*SC;
  double *mat;

  mat = new double [9];
  mat[0] = 1; mat[1] = 2; mat[2] = 3; 
  mat[3] = 4; mat[4] = 5; mat[5] = 6;
  mat[6] = 7; mat[7] = 0; mat[8] = 0; 
  
  cout << "A: " << endl;
  A = new BlkMat(2,2);
  A->GenBlk(0,0,2,2,mat);
  A->GenBlk(1,1,2,2,mat);
  A->PrintBlks();  

  cout << endl << "B: " << endl;
  B = new BlkMat(2,3);
  B->GenBlk(0,0,2,3,mat);
  B->GenBlk(1,0,2,3,mat);
  B->GenBlk(1,1,2,3,mat);
  B->GenBlk(1,2,2,3,mat);
  B->PrintBlks();  

  cout << endl << "C: " << endl;
  C = new BlkMat(3,3);
  C->GenBlk(0,0,3,3,mat);
  C->GenBlk(1,1,3,3,mat);
  C->GenBlk(2,2,3,3,mat);
  C->PrintBlks();
  
  cout << endl << "C^{-1}: " << endl;
  C->invert_diag();
  C->PrintBlks();

  cout << endl << "D: " << endl;
  D = new BlkMat(3,2);
  D->GenBlk(0,0,3,2,mat);
  D->GenBlk(0,1,3,2,mat);
  D->GenBlk(1,1,3,2,mat);
  D->GenBlk(2,1,3,2,mat);
  D->PrintBlks();

  cout << endl << "SC=A-B*C*D: " << endl;
  SC = new BlkMat(2,2);
  T  = new BlkMat(3,2);
  // T->geMxM(RowMajor,RowMajor,1,*C,*D,0);  
  T->MxM(*C,*D);  
  SC->sub(*A,SC->MxM(*B,*T));
  SC->PrintBlks();

  cout << endl << "SC=A-D^T*C*D: " << endl;
  T->MxM(*C,*D); 
  SC->sub(*A,SC->MtxM(*D,*T));
  SC->PrintBlks();

  cout << endl << "SC=A-B*C*B^T: " << endl;
  T->MxMt(*C,*B); 
  //T->geMxM(RowMajor,ColMajor,1,*C,*B,0); 
  SC->sub(*A,SC->MxM(*B,*T));
  SC->PrintBlks();


  delete A;
  delete B;
  delete C;
  delete D;
  delete T;
  delete SC;
  delete[] mat;

  return 0;
}
