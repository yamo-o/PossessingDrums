

/*==============================================================================

    Copyright (c) 2014 Yamo

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
    
    ----------------------------------------------------------------------------
    
    See http://www.yamo-n.org/pd/possessingdrums.html & https://github.com/yamo-o/PossessingDrums 
    for documentation
    
==============================================================================*/


#include "CVector.h"

#ifndef CDENSEMATRIX_HH
#define CDENSEMATRIX_HH


template <class T,class S> class CDenseMatrixContainer{
public:

inline const S _eval(int i,int j) const {
return((static_cast<const T*>(this))->_eval(i,j));
}

};


template <class L,class Op,class R> class CDenseMatrixExpr:public CDenseMatrixContainer<CDenseMatrixExpr<L,Op,R>,typename R::value_type>{
const CDenseMatrixContainer<L,typename L::value_type> &_l;
const CDenseMatrixContainer<R,typename R::value_type> &_r;
public:

typedef typename R::value_type value_type;

CDenseMatrixExpr(const CDenseMatrixContainer<L,typename L::value_type>& l,const CDenseMatrixContainer<R,typename R::value_type>& r):_l(l),_r(r){}

inline const value_type _eval(int i,int j) const {
return(Op::_apply(_l._eval(i,j),_r._eval(i,j)));
}
  
};


template <class Op,class R> class CDenseMatrixExprR:public CDenseMatrixContainer<CDenseMatrixExprR<Op,R>,typename R::value_type>{
const double &_l;
const CDenseMatrixContainer<R,typename R::value_type> &_r;
public:

typedef typename R::value_type value_type;

CDenseMatrixExprR(const double& l,const CDenseMatrixContainer<R,typename R::value_type>& r):_l(l),_r(r){}

inline const value_type _eval(int i,int j) const {
return(Op::_apply(_l,_r._eval(i,j)));
}
};

template <class L,class Op> class CDenseMatrixExprL:public CDenseMatrixContainer<CDenseMatrixExprL<L,Op>,typename L::value_type>{
const CDenseMatrixContainer<L,typename L::value_type> &_l;
const double &_r;
public:

typedef typename L::value_type value_type;

CDenseMatrixExprL(const CDenseMatrixContainer<L,typename L::value_type>& l,const double& r):_l(l),_r(r){}

inline const value_type _eval(int i,int j) const {
return(Op::_apply(_l._eval(i,j),_r));
}
  
inline const int _get_size() const { return(_l._get_size()); }
};



// Matrix + Matrix
template <class L,class R> inline CDenseMatrixExpr<L,OpAdd<typename L::value_type,typename L::value_type,typename R::value_type>,R> operator+(const CDenseMatrixContainer<L,typename L::value_type>& lhs,const CDenseMatrixContainer<R,typename R::value_type>& rhs){
return(CDenseMatrixExpr<L,OpAdd<typename L::value_type,typename L::value_type,typename R::value_type>,R>(lhs,rhs));
}

// scalar + Matrix
template <class R> inline CDenseMatrixExprR<OpAdd<typename R::value_type,double,typename R::value_type>,R> operator+(const double& lhs,const CDenseMatrixContainer<R,typename R::value_type>& rhs){
return(CDenseMatrixExprR<OpAdd<typename R::value_type,double,typename R::value_type>,R>(lhs,rhs));
}

// vector + scalar
template <class L> inline CDenseMatrixExprL<L,OpAdd<typename L::value_type,typename L::value_type,double> > operator+(const CDenseMatrixContainer<L,typename L::value_type>& lhs,const double& rhs){
return(CDenseMatrixExprL<L,OpAdd<typename L::value_type,typename L::value_type,double> >(lhs,rhs));
}

// Matrix - Matrix
template <class L,class R> inline CDenseMatrixExpr<L,OpSub<typename L::value_type,typename L::value_type,typename R::value_type>,R> operator-(const CDenseMatrixContainer<L,typename L::value_type>& lhs,const CDenseMatrixContainer<R,typename R::value_type>& rhs){
return(CDenseMatrixExpr<L,OpSub<typename L::value_type,typename L::value_type,typename R::value_type>,R>(lhs,rhs));
}

// scalar - Matrix
template <class R> inline CDenseMatrixExprR<OpSub<typename R::value_type,double,typename R::value_type>,R> operator-(const double& lhs,const CDenseMatrixContainer<R,typename R::value_type>& rhs){
return(CDenseMatrixExprR<OpSub<typename R::value_type,double,typename R::value_type>,R>(lhs,rhs));
}

// vector - scalar
template <class L> inline CDenseMatrixExprL<L,OpSub<typename L::value_type,typename L::value_type,double> > operator-(const CDenseMatrixContainer<L,typename L::value_type>& lhs,const double& rhs){
return(CDenseMatrixExprL<L,OpSub<typename L::value_type,typename L::value_type,double> >(lhs,rhs));
}

// scalar * Matrix
template <class R> inline CDenseMatrixExprR<OpMul<typename R::value_type,double,typename R::value_type>,R> operator*(const double& lhs,const CDenseMatrixContainer<R,typename R::value_type>& rhs){
return(CDenseMatrixExprR<OpMul<typename R::value_type,double,typename R::value_type>,R>(lhs,rhs));
}

// vector * scalar
template <class L> inline CDenseMatrixExprL<L,OpMul<typename L::value_type,typename L::value_type,double> > operator*(const CDenseMatrixContainer<L,typename L::value_type>& lhs,const double& rhs){
return(CDenseMatrixExprL<L,OpMul<typename L::value_type,typename L::value_type,double> >(lhs,rhs));
}

// Matrix / Matrix
template <class L,class R> inline CDenseMatrixExpr<L,OpDiv<typename L::value_type,typename L::value_type,typename R::value_type>,R> operator/(const CDenseMatrixContainer<L,typename L::value_type>& lhs,const CDenseMatrixContainer<R,typename R::value_type>& rhs){
return(CDenseMatrixExpr<L,OpDiv<typename L::value_type,typename L::value_type,typename R::value_type>,R>(lhs,rhs));
}

// scalar / Matrix
template <class R> inline CDenseMatrixExprR<OpDiv<typename R::value_type,double,typename R::value_type>,R> operator/(const double& lhs,const CDenseMatrixContainer<R,typename R::value_type>& rhs){
return(CDenseMatrixExprR<OpDiv<typename R::value_type,double,typename R::value_type>,R>(lhs,rhs));
}

// Matrix / scalar
template <class L> inline CDenseMatrixExprL<L,OpAdd<typename L::value_type,typename L::value_type,double> > operator/(const CDenseMatrixContainer<L,typename L::value_type>& lhs,const double& rhs){
return(CDenseMatrixExprL<L,OpDiv<typename L::value_type,typename L::value_type,double> >(lhs,rhs));
}


template <class T> class CDenseMatrix:public CDenseMatrixContainer<CDenseMatrix<T>,T>{
public:

typedef T value_type;


int _row,_col;
T** _value;
CVector<T> _pool;

CDenseMatrix(){_value=NULL; }

CDenseMatrix(unsigned int row):_row(row),_col(row){
_pool._resize(this->_row*this->_col);
_pool=0;
_value=new T*[this->_row];
for(int i=0;i<this->_row;++i){
   _value[i]=&_pool[this->_col*i];
}
}

CDenseMatrix(unsigned int row,unsigned int col):_row(row),_col(col){
_pool._resize(this->_row*this->_col);
_pool=0;
_value=new T*[this->_row];
for(int i=0;i<this->_row;++i){
   _value[i]=&_pool[this->_col*i];
}
}

CDenseMatrix(CDenseMatrix<T>& mat):_row(mat._row),_col(mat._row){
_pool._resize(this->_row*this->_col);
_value=new T*[this->_row];
for(int i=0;i<this->_row;++i){
   _value[i]=&_pool[this->_col*i];
}
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]=mat[i][j];
   }
}
}


~CDenseMatrix(){
if(_value){
  delete [] _value;
}
}


virtual void _resize(unsigned int row){
if(_value){
  delete [] _value;
}

_pool._resize(row*row);
this->_row=row;
this->_col=row;
_value=new T*[row];
for(int i=0;i<row;++i){
   _value[i]=&_pool[row*i];
}

}

virtual void _resize(int row,int col){
if(_value){
  delete [] _value;
}

_pool._resize(row*col);
this->_row=row;
this->_col=col;
_value=new T*[row];
for(int i=0;i<row;++i){
   _value[i]=&_pool[col*i];
}

}


inline CVector<T>& _getVector(){ return(_pool); }
inline T* _getPointer(){ return(&_pool[0]); }

inline const T& _eval(int i,int j) const { return(_value[i][j]); }
inline const int _get_size() const { return(this->_row*this->_col); }

inline virtual T _get_value(int row,int col){
return(this->_value[row][col]);
}

inline T _maxDiag(){
T m=this->_value[0][0];
for(int i=1;i<this->_row;++i){
   if(this->_value[i][i]>m) m=this->_value[i][i];
}
return(m);
}

inline T _minDiag(){
T m=this->_value[0][0];
for(int i=1;i<this->_row;++i){
   if(this->_value[i][i]<m) m=this->_value[i][i];
}

return(m);
}

inline void _addDiag(T val){
for(int i=0;i<this->_row;++i){
   this->_value[i][i]+=val;
}
}


void _setBlock(int r,int c,CDenseMatrix<T>* mat){
for(int j=0;j<mat->_row;++j){
   for(int i=0;i<mat->_col;++i){
      _value[r+j][c+i]=(*mat)[j][i];
   }
}
}


inline CDenseMatrix<T>& operator=(const double& scalar){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]=scalar;
   }
}
return(*this);
}

template <class S> inline CDenseMatrix<T>& operator=(const CDenseMatrixContainer<S,typename S::value_type>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]=expr._eval(i,j);
   }
}
return(*this);
}

template <class L,class Op,class R> inline CDenseMatrix<T>& operator=(const CDenseMatrixExpr<L,Op,R>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]=expr._eval(i,j);
   }
}
return(*this);
}

template <class Op,class R> inline CDenseMatrix<T>& operator=(const CDenseMatrixExprR<Op,R>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]=expr._eval(i,j);
   }
}
return(*this);
}

template <class L,class Op> inline CDenseMatrix<T>& operator=(const CDenseMatrixExprL<L,Op>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]=expr._eval(i,j);
   }
}
return(*this);
}

inline CDenseMatrix<T>& operator=(const CDenseMatrix<T>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]=expr._eval(i,j);
   }
}
return(*this);
}


inline CDenseMatrix<T>& operator+=(const double& scalar){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]+=scalar;
   }
}
return(*this);
}

template <class S> inline CDenseMatrix<T>& operator+=(const CDenseMatrixContainer<S,typename S::value_type>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]+=expr._eval(i,j);
   }
}
return(*this);
}

template <class L,class Op,class R> inline CDenseMatrix<T>& operator+=(const CDenseMatrixExpr<L,Op,R>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]+=expr._eval(i,j);
   }
}
return(*this);
}

template <class Op,class R> inline CDenseMatrix<T>& operator+=(const CDenseMatrixExprR<Op,R>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]+=expr._eval(i,j);
   }
}
return(*this);
}

template <class L,class Op> inline CDenseMatrix<T>& operator+=(const CDenseMatrixExprL<L,Op>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]+=expr._eval(i,j);
   }
}
return(*this);
}

inline CDenseMatrix<T>& operator+=(const CDenseMatrix<T>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]+=expr._eval(i,j);
   }
}
return(*this);
}

inline CDenseMatrix<T>& operator-=(const double& scalar){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]-=scalar;
   }
}
return(*this);
}

template <class S> inline CDenseMatrix<T>& operator-=(const CDenseMatrixContainer<S,typename S::value_type>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]-=expr._eval(i,j);
   }
}
return(*this);
}

template <class L,class Op,class R> inline CDenseMatrix<T>& operator-=(const CDenseMatrixExpr<L,Op,R>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]-=expr._eval(i,j);
   }
}
return(*this);
}

template <class Op,class R> inline CDenseMatrix<T>& operator-=(const CDenseMatrixExprR<Op,R>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]-=expr._eval(i,j);
   }
}
return(*this);
}

template <class L,class Op> inline CDenseMatrix<T>& operator-=(const CDenseMatrixExprL<L,Op>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]-=expr._eval(i,j);
   }
}
return(*this);
}

inline CDenseMatrix<T>& operator-=(const CDenseMatrix<T>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]-=expr._eval(i,j);
   }
}
return(*this);
}

inline CDenseMatrix<T>& operator*=(const double scalar){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]*=scalar;
   }
}
return(*this);
}

template <class S> inline CDenseMatrix<T>& operator*=(const CDenseMatrixContainer<S,typename S::value_type>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]*=expr._eval(i,j);
   }
}
return(*this);
}

template <class L,class Op,class R> inline CDenseMatrix<T>& operator*=(const CDenseMatrixExpr<L,Op,R>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]*=expr._eval(i,j);
   }
}
return(*this);
}

template <class Op,class R> inline CDenseMatrix<T>& operator*=(const CDenseMatrixExprR<Op,R>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]*=expr._eval(i,j);
   }
}
return(*this);
}

template <class L,class Op> inline CDenseMatrix<T>& operator*=(const CDenseMatrixExprL<L,Op>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]*=expr._eval(i,j);
   }
}
return(*this);
}

inline CDenseMatrix<T>& operator*=(const CDenseMatrix<T>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]*=expr._eval(i,j);
   }
}
return(*this);
}

inline CDenseMatrix<T>& operator/=(const double& scalar){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]/=scalar;
   }
}
return(*this);
}

template <class S> inline CDenseMatrix<T>& operator/=(const CDenseMatrixContainer<S,typename S::value_type>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]/=expr._eval(i,j);
   }
}
}

template <class L,class Op,class R> inline CDenseMatrix<T>& operator/=(const CDenseMatrixExpr<L,Op,R>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]/=expr._eval(i,j);
   }
}
return(*this);
}

template <class Op,class R> inline CDenseMatrix<T>& operator/=(const CDenseMatrixExprR<Op,R>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]/=expr._eval(i,j);
   }
}
return(*this);
}

template <class L,class Op> inline CDenseMatrix<T>& operator/=(const CDenseMatrixExprL<L,Op>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]/=expr._eval(i,j);
   }
}
return(*this);
}

inline CDenseMatrix<T>& operator/=(const CDenseMatrix<T>& expr){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<this->_col;++j){
      _value[i][j]/=expr._eval(i,j);
   }
}
return(*this);
}

void _zero(){
for(int i=0;i<this->_row;i++){
   memset(_value[i],0,this->_col*sizeof(T));
}
}

virtual void _zeroRow(int row){
memset(_value[row],0,this->_col*sizeof(T));
}

void _identity(){
_zero();
for(int i=0;i<this->_row;i++){
   _value[i][i]=1;
}
}

bool _isSymmetric(){
for(int i=0;i<this->_row;++i){
   for(int j=0;j<i;++j){
      if(fabs(_value[i][j]-_value[j][i])>1.0e-3) return(false);
      
   }
}
return(true);
}

bool _isPositiveDiagonal(){
for(int i=0;i<this->_row;++i){
   if(_value[i][i]<0) return(false);
}
return(true);
}


inline T &operator()(int i,int j){
return(_value[i][j]);
}

inline const T &operator()(int i,int j) const {
return(_value[i][j]);
}

inline T* operator[](int i){
return(_value[i]);
}

inline const T* operator[](int i) const {
return(_value[i]);
}

inline void _set(int i,int j,T val){
_value[i][j]=val;
}

inline T *_get_vec(int row){
return(_value[row]);
}

inline void _copy_vec(int row,T* vec){
memcpy(_value[row],vec,this->_col*sizeof(T));
}


void _symmetrize(bool upper=true){

if(upper){
  for(int j=0;j<this->_row;++j){
     for(int i=j+1;i<this->_col;++i){
        _value[i][j]=_value[j][i];
     }
  }
}else{
  for(int j=1;j<this->_row;++j){
     for(int i=0;i<j;++i){
        _value[i][j]=_value[j][i];
     }
  }
}

}

void _transpose(CDenseMatrix<T>& m){
for(int j=0;j<this->_row;++j){
   for(int i=0;i<this->_col;++i){
      m[i][j]=_value[j][i];
   }
}
}

T _trace(){
T d=0;
for(int j=0;j<this->_row;++j) d+=_value[j][j];
return(d);
}

};



#endif
