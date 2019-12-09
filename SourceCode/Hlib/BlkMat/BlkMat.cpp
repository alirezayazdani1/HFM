#include "BlkSubMat.h"
#include "BlkVec.h"
#include "BlkMat.h"

namespace NekBlkMat {

  void BlkVec::GenBlk(const int id, const int rows, const int cols, 
		      const double *mat){
    int en;
    BlkSubMat *B;
    
    if((en = _entries[id])+1){// entry already exists so add matrix 
      _Blk[en]->AddMat(rows,cols,mat);
    }
    else{ // add new entry
      B = new BlkSubMat(id,rows,cols,mat);
      
      _entries[id] = _Blk.size();
      
      _Blk.push_back(B);
    }
  }  

  void BlkVec::GenBlk(const int id, const int rows, const int cols, 
		      const int id1, const int id2, const double val){
    int en;
    BlkSubMat *B;
    
    if((en = _entries[id])+1){// entry already exists so add matrix 
	if((_Blk[en]->get_rows() != rows)||(_Blk[en]->get_cols() != cols))
	  NekError::error(warning,"BlkVec::GenBlk","add value to Blk which "
		  "does not have the same dimensions as calling argument");
      _Blk[en]->AddVal(id1,id2,val);
    }
    else{ // add new entry
      B = new BlkSubMat(id,rows,cols);
      B->AddVal(id1,id2,val);

      _entries[id] = _Blk.size();
      
      _Blk.push_back(B);
    }
  }  

  // (this) = a + b
  BlkVec& BlkVec::add(const BlkVec& a, const BlkVec& b){
    int i,id,en1,en,cnt=0;
    BlkSubMat *B;

    if(this == &b){ // need to to b-matrix first 

      // check over 'a' row entries first 
      for(i = 0; i < a._Blk.size(); ++i){
	if((en = _entries[id = a._Blk[i]->get_id()])+1){
	  *(_Blk[en]) += *(a._Blk[i]);
	}
	else{ // need to declare memory and copy  in a matrix
	  B = new BlkSubMat(*a._Blk[i]);
	  en = _entries[id] = _Blk.size();
	  _Blk.push_back(B);
	}
      }
    }
    else{

      // check over 'a' row entries first 
      for(i = 0; i < a._Blk.size(); ++i){
	if((en = _entries[id = a._Blk[i]->get_id()])+1)
	  *(_Blk[en]) = *(a._Blk[i]);
	else{ // need to declare memory and copy  in a matrix
	  B = new BlkSubMat(*a._Blk[i]);
	  en = _entries[id] = _Blk.size();
	  _Blk.push_back(B);
	}
	
	//check to see if b entry exists and if so add 
	if((en1 = b._entries[id])+1) 
	  *(_Blk[en]) += *(b._Blk[en1]);
      }
      
      cnt = a._Blk.size();
      
      // check over 'b' vector  entries for any independent entries not in 'a'
      for(i = 0; i < b._Blk.size(); ++i){
	if((a._entries[id = b._Blk[i]->get_id()]) == -1){
	  if((en = _entries[id])+1)
	    *(_Blk[en]) = *(b._Blk[i]);
	  else{
	    B = new BlkSubMat(*b._Blk[i]);
	    _entries[id] = _Blk.size();
	    _Blk.push_back(B); 
	  }
	  ++cnt;
	}
      }
    

     //If original row is larger than new contributions then delete extra terms
      if(_Blk.size() != cnt){ 
	for(i = 0; i < _Blk.size(); ++i){
	  id = _Blk[i]->get_id();
	  if((a._entries[id] == -1)&&(b._entries[id] == -1)){
	    // remove element
	    B = _Blk[i] ;
	    _Blk.erase(_Blk.begin()+i);
	    _entries[i] = -1;
	    delete B;
	  }
	}
	if(_Blk.size() != cnt)
	  NekError::error(fatal,"BlkVec::add","incorrect final size");
      } 
    }

    return *this;
  }

  // (this) = a - b
  BlkVec& BlkVec::sub(const BlkVec& a, const BlkVec& b){
    int i,id,en1,en,cnt;
    BlkSubMat *B;
    
    if(this == &b){ // need to do b-matrix first 

      // Negate b entries 
      for(i = 0; i < b._Blk.size(); ++i)
	_Blk[i]->neg();
	
      
      // add 'a'  entries 
      for(i = 0; i < a._Blk.size(); ++i){
	if((en = _entries[id = a._Blk[i]->get_id()])+1)
	  *(_Blk[en]) += *(a._Blk[i]);
	else{ // need to declare memory and copy  in a matrix
	  B = new BlkSubMat(*a._Blk[i]);
	  en = _entries[id] = _Blk.size();
	  _Blk.push_back(B);
	}
      }
      
    }
    else{
      // check over 'a' row entries first 
      for(i = 0; i < a._Blk.size(); ++i){
	if((en = _entries[id = a._Blk[i]->get_id()])+1)
	  *(_Blk[en]) = *(a._Blk[i]);
	else{ // need to declare memory and copy  in a matrix
	  B = new BlkSubMat(*a._Blk[i]);
	  en = _entries[id] = _Blk.size();
	  _Blk.push_back(B);
	}
	//check to see if b entry exists and if so add 
	if((en1 = b._entries[id])+1) 
	  *(_Blk[en]) -= *(b._Blk[en1]);
      }
      
      cnt = a._Blk.size();
      
      // check over 'b' row entries for any independent entries 
      for(i = 0; i < b._Blk.size(); ++i){
	if((a._entries[id = b._Blk[i]->get_id()]) == -1){
	  if((en = _entries[id])+1){
	    *(_Blk[en]) = *(b._Blk[i]);
	    _Blk[en]->neg();
	  }
	else{
	  B = new BlkSubMat(*b._Blk[i]);
	  B->neg();
	  _entries[id] = _Blk.size();
	  _Blk.push_back(B); 
	}
	  ++cnt;
	}
      }
      
      // If original row is larger than new contributions then delete extra terms
      if(_Blk.size() != cnt){ 
	for(i = 0; i < _Blk.size(); ++i){
	  id = _Blk[i]->get_id();
	  if((a._entries[id] == -1)&&(b._entries[id] == -1)){
	    // remove element
	    B = _Blk[i] ;
	    _Blk.erase(_Blk.begin()+i);
	    _entries[i] = -1;
	    delete B;
	  }
	}
	if(_Blk.size() != cnt)
	  NekError::error(fatal,"BlkVec::sub","incorrect final size");
      } 
    }

    return *this;
  }


  // (this) += a*v where 'a' is a submatrix 
  BlkVec& BlkVec::axpy(const double alpha, const BlkVec& A){
    int en,i,id;

    // check over 'A' row entries 
    for(i = 0; i < A._Blk.size(); ++i){
      if((en = _entries[id = A._Blk[i]->get_id()])+1)
	_Blk[en]->axpy(alpha,*(A._Blk[i]));
      else{ // need to declare memory and copy in a matrix
	BlkSubMat *B; 
	B = new BlkSubMat(*A._Blk[i]);
	B->scal(alpha);
	en = _entries[id] = _Blk.size();
	_Blk.push_back(B);
      }
    }
      
    return *this;
  }



  // y += a*v where a is a submatrix 
  double* BlkVec::geMxv(MatStorage form, const double& alpha, const double *v,
			const double& beta, double *y, const int *offset){
    int i,id;
    
    switch(form){
    case RowMajor:

      // Multiply a by v and add to  y 
      for(i = 0; i < _Blk.size(); ++i){
	id = _Blk[i]->get_id();
	_Blk[i]->geMxv(form,alpha,v+offset[id],beta,y);
      }
      break;
    case ColMajor:

      // Multiply a by v and add to  y 
      for(i = 0; i < _Blk.size(); ++i){
	id = _Blk[i]->get_id();
	_Blk[i]->geMxv(form,alpha,v,beta,y+offset[id]);
      }
      break;
    default:
      NekError::error(warning,"BlkMat::geMxv",
		      "matrix format not specified for matrix A");
      break;
    }
    
    return y;
  }


  // y += a*v where a is a submatrix 
  double* BlkVec::Mxvpy(const int *offset, const double *v, double *y){
    int i,id;
    
    // Multiply a by v and add to  y 
    for(i = 0; i < _Blk.size(); ++i){
      id = _Blk[i]->get_id();
      _Blk[i]->Mxvpy(v+offset[id],y);
    }
    
    return y;
  }

  // y += a^t*v ywhere a is a submatrix 
  double* BlkVec::Mtxvpy(const int *offset, const double *v, double *y){
    int i,id;
    
    // Multiply a by v and add to  y 
    for(i = 0; i < _Blk.size(); ++i){
      id = _Blk[i]->get_id();
      _Blk[i]->Mtxvpy(v,y+offset[id]);
    }
    
    return y;
  }


  // (this) += a*b where a is a submatrix 
  BlkVec& BlkVec::MxMpM(const BlkSubMat& a, const BlkVec& b){
    int i,id,en;
    
    // Multiply a by all b entries and put into (*this)
    for(i = 0; i < b._Blk.size(); ++i){
      if(b._Blk[i]->get_rows() != a.get_cols())
	NekError::error(warning,"BlkVec::MxMpM","cols and rows do not match");

      if((en = _entries[id = b._Blk[i]->get_id()])+1)
	_Blk[en]->MxMpM(a,*(b._Blk[i]));
      else{ // need to declare memory and copy  in a matrix
	BlkSubMat *B;
	B = new BlkSubMat(id,a.get_rows(),b._Blk[i]->get_cols());
	B->MxMpM(a,*(b._Blk[i]));
	en = _entries[id] = _Blk.size();
	_Blk.push_back(B);
      }
    }
    
    return *this;
  }


  // (this) += a^t*b where 'a' is a submatrix 
  BlkVec& BlkVec::MtxMpM(const BlkSubMat& a, const BlkVec& b){
    int i,id,en;
    
    // Multiply a by all b entries and put into (*this)
    for(i = 0; i < b._Blk.size(); ++i){
      if(b._Blk[i]->get_rows() != a.get_rows())
	NekError::error(warning,"BlkVec::axbpy","cols and rows do not match");

      if((en = _entries[id = b._Blk[i]->get_id()])+1)
	_Blk[en]->MtxMpM(a,*(b._Blk[i]));
      else{ // need to declare memory and copy  in a matrix
	BlkSubMat *B;
	B = new BlkSubMat(id,a.get_cols(),b._Blk[i]->get_cols());
	B->MtxMpM(a,*(b._Blk[i]));
	en = _entries[id] = _Blk.size();
	_Blk.push_back(B);
      }
    }

    return *this;
  }


  //----------------------------------------------------------------

  // (this) = a + b
  BlkMat& BlkMat::add(const BlkMat& a, const BlkMat& b){
    
    if((a._max_rowblk != _max_rowblk)||(b._max_rowblk != _max_rowblk)||
       (a._max_colblk != _max_colblk)||(b._max_colblk != _max_colblk))
      NekError::error(warning,"BlkMat::sub","Matrix rows/cols not compatible");
    
    for(int i=0; i < _max_rowblk; ++i)
      _Row[i].add(a._Row[i],b._Row[i]);
    
    return *this;
  }

  // (this) = a - b
  BlkMat& BlkMat::sub(const BlkMat& a, const BlkMat& b){
    
    if((a._max_rowblk != _max_rowblk)||(b._max_rowblk != _max_rowblk)||
       (a._max_colblk != _max_colblk)||(b._max_colblk != _max_colblk))
      NekError::error(warning,"BlkMat::sub","Matrix rows/cols not compatible");
    
    for(int i=0; i < _max_rowblk; ++i)
      _Row[i].sub(a._Row[i],b._Row[i]);
    
    return *this;
  }


  // (this) = (this) + alpha*A
  BlkMat& BlkMat::axpy(const double alpha, const BlkMat& A){
    
    if((A._max_rowblk != _max_rowblk)||(A._max_colblk != _max_colblk))
      NekError::error(warning,"BlkMat::axMpM","Matrix rows/cols not compatible");
    
    for(int i=0; i < _max_rowblk; ++i)
      _Row[i].axpy(alpha,A._Row[i]);
    
    return *this;
  }


  // y = alpha A*v + beta y
  double* BlkMat::geMxv(MatStorage form, const double& alpha, const double *v,
			const double& beta, double *y){
    int i,n;

    if(!_offset) setup_offset();

    switch(form){
    case RowMajor:
      // multiply every row by v 
      for(n=i=0; i < _max_rowblk; ++i){
	_Row[i].geMxv(form,alpha,v,beta,y+n,_offset);
	n += _Row[i]._Blk[0]->get_rows();
      }
      break;
    case ColMajor:
      // multiply every row by v 
      for(n=i=0; i < _max_rowblk; ++i){
	_Row[i].geMxv(form,alpha,v+n,beta,y,_offset);
	n += _Row[i]._Blk[0]->get_rows();
      }
      break;
    default:
      NekError::error(warning,"BlkMat::geMxv",
			"matrix format not specified for matrix A");
      break;
    }
    
    return y;
  }

  // y  = A * v + y
  double* BlkMat::Mxvpy(double* v, double* y){
    int n,i;
    
    if(!_offset) setup_offset();
    
    // multiply every row by v 
    for(n=i=0; i < _max_rowblk; ++i){
      _Row[i].Mxvpy(_offset,v,y+n);
      n += _Row[i]._Blk[0]->get_rows();
    }
    
    return y;
  }

  // y  = A^t * v + y
  double* BlkMat::Mtxvpy(double* v, double* y){
    int n,i;
    
    if(!_offset) setup_offset();
    
    // multiply every row by v 
    for(n=i=0; i < _max_rowblk; ++i){
      _Row[i].Mtxvpy(_offset,v+n,y);
      n += _Row[i]._Blk[0]->get_rows();
    }
    return y;
  }

  // (this) = alpha*A*B + beta (this)  
  BlkVec& BlkVec::geMxM(MatStorage formA, MatStorage formB,
			const double& alpha, const BlkSubMat& A, 
			const BlkVec& B, const double& beta){
    int i,id,en;
    
    // Multiply A by all b entries and put into (*this)
    for(i = 0; i < B._Blk.size(); ++i){
      if((en = _entries[id = B._Blk[i]->get_id()])+1)
	_Blk[en]->geMxM(formA,formB,alpha,A,*(B._Blk[i]),beta);
      else{ // need to declare memory and copy  in a matrix
	BlkSubMat *Btmp;
	Btmp = new BlkSubMat(id,A.get_rows(),B._Blk[i]->get_cols());
	Btmp->geMxM(formA,formB,alpha,A,*(B._Blk[i]),beta);
	en = _entries[id] = _Blk.size();
	_Blk.push_back(Btmp);
      }
    }
    return *this;
  }

  // (this) =  A *(this) where A is block diagonal 
  BlkMat& BlkMat::diagMxy(const BlkMat& A){
    int i,j,rows,cols;
    BlkSubMat *Btmp = (BlkSubMat *)NULL;

    rows = cols = -1;

    // multiply every entry in (this) by matrix diagonal of A
    for(i=0; i < _max_rowblk; ++i){
      // check rows are the same 
      if(A._Row[i]._Blk[0]->get_rows() != 
	 _Row[i]._Blk[0]->get_rows())
	NekError::error(fatal,"Blkmat::diagMxy","Rows are not the same");
      
      for(j=0;j < _Row[i]._Blk.size(); ++j){
	// set up temporary matrix 
	if((rows != _Row[i]._Blk[j]->get_rows())
	   ||(cols != _Row[i]._Blk[j]->get_cols())){
	  if(Btmp) delete Btmp;
	  rows = _Row[i]._Blk[j]->get_rows();
	  cols = _Row[i]._Blk[j]->get_cols();
	  Btmp = new BlkSubMat(j,rows,cols);
	}
	
	// do multiplication
	Btmp->geMxM(RowMajor,RowMajor,1.0,
		    A._Row[i]._Blk[A._Row[i]._entries[i]][0],
		    _Row[i]._Blk[j][0],0.0);

	// copy into original storage; 
	*(_Row[i]._Blk[j]) = *Btmp;
      }
    }
    if(Btmp) delete Btmp;
    
    return *this;
  }


  // (this) = alpha* A * B + beta (this)
  BlkMat& BlkMat::geMxM(MatStorage formA, MatStorage formB,const double& alpha,
			const BlkMat& A, const BlkMat& B, const double& beta){
    int id,ida,idb,i,j,k;
    BlkSubMat *Btmp;

    // if beta == 0 then initalise (*this) matrix by  deleting row blocks 
    if(!beta)
      Reset_Rows();
    
    switch(formB){
    case RowMajor:
      switch(formA){
      case RowMajor:
	// multiply every entry in row of B by A matrix and put in C[i]
	for(i=0; i < A._max_rowblk; ++i)
	  for(j=0;j < A._Row[i]._Blk.size(); ++j){
	    id = A._Row[i]._Blk[j]->get_id();
	    _Row[i].geMxM(formA,formB,alpha,A._Row[i]._Blk[j][0],
			  B._Row[id],beta);
	  }
	break;
      case ColMajor:
	// multiply every entry in row of B by A matrix and put in C[id]
	for(i=0; i < A._max_rowblk; ++i)
	  for(j=0;j < A._Row[i]._Blk.size(); ++j){
	    id = A._Row[i]._Blk[j]->get_id();
	    _Row[id].geMxM(formA,formB,alpha,A._Row[i]._Blk[j][0],
			   B._Row[i],beta);
	  }
	break;
      default:
	NekError::error(warning,"BlkMat::geMxM",
			"matrix format not specified for matrix A");
	break;
      }
      break;
    case ColMajor:
      switch(formA){
      case RowMajor:
	// Can not use BlkVec since Blk row ordering is not consistent
	for(i=0; i < A._max_rowblk; ++i)
	  for(j=0;j < B._max_rowblk; ++j)
	    
	    for(k = 0; k < A._max_colblk; ++k){
	      if(((ida=A._Row[i]._entries[k])+1)
		 &&((idb=B._Row[j]._entries[k])+1))
		if((id = _Row[i]._entries[j])+1) // add local contribution
		  _Row[i]._Blk[id]->geMxM(formA,formB,alpha,
		  A._Row[i]._Blk[ida][0],B._Row[j]._Blk[idb][0],beta);
		else{ // declare and put in matrix 
		  Btmp = new BlkSubMat(j,A._Row[i]._Blk[ida]->get_rows(),
				    B._Row[j]._Blk[idb]->get_rows());
		  Btmp->geMxM(formA,formB,alpha,A._Row[i]._Blk[ida][0],
			   B._Row[j]._Blk[idb][0],beta);

		  _Row[i]._entries[j] = _Row[i]._Blk.size();
		  _Row[i]._Blk.push_back(Btmp);
		}
	    }
	break;
      case ColMajor:
	// Can not use BlkVec since Blk row ordering is not consistent
	for(i=0; i < A._max_colblk; ++i)
	  for(j=0;j < B._max_rowblk; ++j)
	    
	    for(k = 0; k < A._max_rowblk; ++k){
	      if(((ida=A._Row[k]._entries[i])+1)
		 &&((idb=B._Row[j]._entries[k])+1))
		if((id = _Row[i]._entries[j])+1) // add local contribution
		  _Row[i]._Blk[id]->geMxM(formA,formB,alpha,
		  A._Row[k]._Blk[ida][0],B._Row[j]._Blk[idb][0],beta);
		else{ // declare and put in matrix 
		  Btmp = new BlkSubMat(j,A._Row[k]._Blk[ida]->get_cols(),
				    B._Row[j]._Blk[idb]->get_rows());
		  Btmp->geMxM(formA,formB,alpha,A._Row[k]._Blk[ida][0],
			   B._Row[j]._Blk[idb][0],beta);

		  _Row[i]._entries[j] = _Row[i]._Blk.size();
		  _Row[i]._Blk.push_back(Btmp);
		}
	    }
	break;
      default:
	NekError::error(warning,"BlkMat::geMxM",
			"matrix format not specified for matrix A");
	break;
      }
      break;
    default:
      NekError::error(warning,"BlkMat::geMxM",
		      "matrix format not specified for matrix B");
      break;
    }
      
    if(!beta)
      setup_offset();

    return *this;
  }


  // (this) = a * b
  BlkMat& BlkMat::MxM(const BlkMat& a, const  BlkMat& b){
    int id,i,j;

    if((a._max_colblk != b._max_rowblk)||(a._max_rowblk != _max_rowblk)||
       (b._max_colblk != _max_colblk))
      NekError::error(warning,"BlkMat::mult",
		      "Matrix rows/cols are not compatible");

    // initalise (*this) matrix - delete row blocks 
    Reset_Rows();
    
    // multiply every entry in row in each entry 
    for(i=0; i < a._max_rowblk; ++i)
      for(j=0;j < a._Row[i]._Blk.size(); ++j){
	id = a._Row[i]._Blk[j]->get_id();
	_Row[i].MxMpM(a._Row[i]._Blk[j][0],b._Row[id]);
      }

    setup_offset();
    return *this;
  }

  // (this) = a^T * b
  BlkMat& BlkMat::MtxM(const BlkMat& a, const BlkMat& b){
    int id,i,j;

    if((a._max_rowblk != b._max_rowblk)||(a._max_colblk != _max_rowblk)||
       (b._max_colblk != _max_colblk))
      NekError::error(warning,"BlkMat::mult",
		      "Matrix rows/cols are not compatible");

    // initalise (*this) matrix - delete row blocks 
    Reset_Rows();

    // multiply every entry in row in each entry 
    for(i=0; i < a._max_rowblk; ++i)
      for(j=0;j < a._Row[i]._Blk.size(); ++j){
	id = a._Row[i]._Blk[j]->get_id();
	_Row[id].MtxMpM(a._Row[i]._Blk[j][0],b._Row[i]);
      }

    setup_offset();
    return *this;
  }


  // (this) = a * b^t
  BlkMat& BlkMat::MxMt(const BlkMat& a, const BlkMat& b){
    int id,ida,idb,i,j,k;

    if((a._max_colblk != b._max_colblk)||(a._max_rowblk != _max_rowblk)||
       (b._max_rowblk != _max_colblk))
      NekError::error(warning,"BlkMat::mult",
		      "Matrix rows/cols are not compatible");
      

    // initalise (*this) matrix - delete row blocks 
    Reset_Rows();
    

    // Can not use BlkVec since Blk row ordering is not consistent
    for(i=0; i < a._max_rowblk; ++i)
      for(j=0;j < b._max_rowblk; ++j)

	for(k = 0; k < a._max_colblk; ++k){
	  if(((ida=a._Row[i]._entries[k])+1)&&((idb=b._Row[j]._entries[k])+1))
	    if((id = _Row[i]._entries[j])+1) // add local contribution
	      _Row[i]._Blk[id]->MxMtpM(a._Row[i]._Blk[ida][0],
				       b._Row[j]._Blk[idb][0]);
	    else{ // declare and put in matrix 
	      BlkSubMat *B;
	      B = new BlkSubMat(j,a._Row[i]._Blk[ida]->get_rows(),
				b._Row[j]._Blk[idb]->get_rows());
	      B->MxMtpM(a._Row[i]._Blk[ida][0],b._Row[j]._Blk[idb][0]);
	      _Row[i]._entries[j] = _Row[i]._Blk.size();
	      _Row[i]._Blk.push_back(B);
	    }
	}

    setup_offset();
    return *this;
  }

  // invert diagonal matrices 
  BlkMat& BlkMat::invert_diag(){
    int i,id;

    for(i = 0; i < _max_rowblk; ++i){
      if(_Row[i]._Blk.size() != 1)
	NekError::error(warning,"BlkMat::invert","matrix not diagonal");

      if((id = _Row[i]._entries[i])+1)
	_Row[i]._Blk[id]->invert();	
      else
	NekError::error(warning,"BlkMat::invert","No diagonal component");
    }
    return *this;
  }

  
  double BlkMat::AmAt(){
    int i,j,k,l,nr,nc;
    double  sum = 0;
    double *mat1,*mat2;

    for(i = 0; i < _max_rowblk; ++i)
      for(j = 0; j < _max_colblk; ++j)
	if(mat1 = get_mat(i,j,nr,nc))
	  if(mat2 = get_mat(j,i,nr,nc))
	    for(k = 0; k < nr; ++k)
	      for(l = 0; l < nc; ++l)
		sum += mat1[k*nc+l] - mat2[l*nc+k];
    return sum;
  }
}
  
