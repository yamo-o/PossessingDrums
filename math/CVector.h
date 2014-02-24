

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


#include "CVector3d.h"

#ifndef CVECTOR_HH
#define CVECTOR_HH



template <class T,class S> class CVectorContainer{
public:

inline const S _eval(int i) const {
return((static_cast<const T*>(this))->_eval(i));
}

inline const int _get_size() const {
return((static_cast<const T*>(this))->_get_size());
}

};


template <class L,class Op,class R> class CVectorExpr:public CVectorContainer<CVectorExpr<L,Op,R>,typename R::value_type>{
const CVectorContainer<L,typename L::value_type> &_l;
const CVectorContainer<R,typename R::value_type> &_r;
public:

typedef typename R::value_type value_type;

CVectorExpr(const CVectorContainer<L,typename L::value_type>& l,const CVectorContainer<R,typename R::value_type>& r):_l(l),_r(r){}

inline const value_type _eval(int i) const {
return(Op::_apply(_l._eval(i),_r._eval(i)));
}
  
inline const int _get_size() const { return(_r._get_size()); }
};

template <class Op,class R> class CVectorExprR:public CVectorContainer<CVectorExprR<Op,R>,typename R::value_type>{
const double &_l;
const CVectorContainer<R,typename R::value_type> &_r;
public:

typedef typename R::value_type value_type;

CVectorExprR(const double& l,const CVectorContainer<R,typename R::value_type>& r):_l(l),_r(r){}

inline const value_type _eval(int i) const {
return(Op::_apply(_l,_r._eval(i)));
}
  
inline const int _get_size() const { return(_r._get_size()); }
};

template <class L,class Op> class CVectorExprL:public CVectorContainer<CVectorExprL<L,Op>,typename L::value_type>{
const CVectorContainer<L,typename L::value_type> &_l;
const double &_r;
public:

typedef typename L::value_type value_type;

CVectorExprL(const CVectorContainer<L,typename L::value_type>& l,const double& r):_l(l),_r(r){}

inline const value_type _eval(int i) const {
return(Op::_apply(_l._eval(i),_r));
}
  
inline const int _get_size() const { return(_l._get_size()); }
};


// vector + vector
template <class L,class R> inline CVectorExpr<L,OpAdd<typename L::value_type,typename L::value_type,typename R::value_type>,R> operator+(const CVectorContainer<L,typename L::value_type>& lhs,const CVectorContainer<R,typename R::value_type>& rhs){
return(CVectorExpr<L,OpAdd<typename L::value_type,typename L::value_type,typename R::value_type>,R>(lhs,rhs));
}

// scalar + vector
template <class R> inline CVectorExprR<OpAdd<typename R::value_type,double,typename R::value_type>,R> operator+(const double& lhs,const CVectorContainer<R,typename R::value_type>& rhs){
return(CVectorExprR<OpAdd<typename R::value_type,double,typename R::value_type>,R>(lhs,rhs));
}

// vector + scalar
template <class L> inline CVectorExprL<L,OpAdd<typename L::value_type,typename L::value_type,double> > operator+(const CVectorContainer<L,typename L::value_type>& lhs,const double& rhs){
return(CVectorExprL<L,OpAdd<typename L::value_type,typename L::value_type,double> >(lhs,rhs));
}

// vector - vector
template <class L,class R> inline CVectorExpr<L,OpSub<typename L::value_type,typename L::value_type,typename R::value_type>,R> operator-(const CVectorContainer<L,typename L::value_type>& lhs,const CVectorContainer<R,typename R::value_type>& rhs){
return(CVectorExpr<L,OpSub<typename L::value_type,typename L::value_type,typename R::value_type>,R>(lhs,rhs));
}

// scalar - vector
template <class R> inline CVectorExprR<OpSub<typename R::value_type,double,typename R::value_type>,R> operator-(const double& lhs,const CVectorContainer<R,typename R::value_type>& rhs){
return(CVectorExprR<OpSub<typename R::value_type,double,typename R::value_type>,R>(lhs,rhs));
}

// vector - scalar
template <class L> inline CVectorExprL<L,OpSub<typename L::value_type,typename L::value_type,double> > operator-(const CVectorContainer<L,typename L::value_type>& lhs,const double& rhs){
return(CVectorExprL<L,OpSub<typename L::value_type,typename L::value_type,double> >(lhs,rhs));
}

// vector * vector
template <class L,class R> inline CVectorExpr<L,OpMul<typename L::value_type,typename L::value_type,typename R::value_type>,R> operator*(const CVectorContainer<L,typename L::value_type>& lhs,const CVectorContainer<R,typename R::value_type>& rhs){
return(CVectorExpr<L,OpMul<typename L::value_type,typename L::value_type,typename R::value_type>,R>(lhs,rhs));
}

// scalar * vector
template <class R> inline CVectorExprR<OpMul<typename R::value_type,double,typename R::value_type>,R> operator*(const double& lhs,const CVectorContainer<R,typename R::value_type>& rhs){
return(CVectorExprR<OpMul<typename R::value_type,double,typename R::value_type>,R>(lhs,rhs));
}

// vector * scalar
template <class L> inline CVectorExprL<L,OpMul<typename L::value_type,typename L::value_type,double> > operator*(const CVectorContainer<L,typename L::value_type>& lhs,const double& rhs){
return(CVectorExprL<L,OpMul<typename L::value_type,typename L::value_type,double> >(lhs,rhs));
}

// vector / vector
template <class L,class R> inline CVectorExpr<L,OpDiv<typename L::value_type,typename L::value_type,typename R::value_type>,R> operator/(const CVectorContainer<L,typename L::value_type>& lhs,const CVectorContainer<R,typename R::value_type>& rhs){
return(CVectorExpr<L,OpDiv<typename L::value_type,typename L::value_type,typename R::value_type>,R>(lhs,rhs));
}

// scalar / vector
template <class R> inline CVectorExprR<OpDiv<typename R::value_type,double,typename R::value_type>,R> operator/(const double& lhs,const CVectorContainer<R,typename R::value_type>& rhs){
return(CVectorExprR<OpDiv<typename R::value_type,double,typename R::value_type>,R>(lhs,rhs));
}

// vector / scalar
template <class L> inline CVectorExprL<L,OpDiv<typename L::value_type,typename L::value_type,double> > operator/(const CVectorContainer<L,typename L::value_type>& lhs,const double& rhs){
return(CVectorExprL<L,OpDiv<typename L::value_type,typename L::value_type,double> >(lhs,rhs));
}



template <class T> class CVector:public CVectorContainer<CVector<T>,T>{
public:

typedef T value_type;

T* _value;
int _size;

CVector(){
_size=0;
_value=NULL;
}

explicit CVector(int size){
_size=size;
_value=new T[size];
}

template <class S> CVector(const CVectorContainer<S,typename S::value_type>& expr){
_size=expr._get_size();
_value=new T[_size];
for(int i=0;i<_size;++i){
   _value[i]=expr._eval(i);
}
}

~CVector(){
if(_value) delete [] _value;

}

inline T& operator[](int i){ return(_value[i]); }
inline const T& operator[](int i) const { return(_value[i]); }
inline T& operator()(int &i){ return(_value[i]); }
inline const T& operator()(int &i) const { return(_value[i]); }
inline const T& _eval(int i) const { return(_value[i]); }
inline const int _get_size() const { return(_size); }
inline T* _getPointer(){ return(_value); }


inline CVector<T>& operator=(const T& scalar){
for(int i=0;i<_size;++i){
   _value[i]=scalar;
}
return(*this);
}


template <class S> inline CVector<T>& operator=(const CVectorContainer<S,typename S::value_type>& expr){
if(!_value){
  _size=expr._get_size();
  _value=new T[_size];
}else if(_size!=expr._get_size()){
  delete [] _value;
  _size=expr._get_size();
  _value=new T[_size];
}

for(int i=0;i<_size;++i){
   _value[i]=expr._eval(i);
}
return(*this);
}

template <class L,class Op,class R> inline CVector<T>& operator=(const CVectorExpr<L,Op,R>& expr){
if(!_value){
  _size=expr._get_size();
  _value=new T[_size];
}else if(_size!=expr._get_size()){
  delete [] _value;
  _size=expr._get_size();
  _value=new T[_size];
}


for(int i=0;i<_size;++i){
   _value[i]=expr._eval(i);
}
return(*this);
}

template <class Op,class R> inline CVector<T>& operator=(const CVectorExprR<Op,R>& expr){
if(!_value){
  _size=expr._get_size();
  _value=new T[_size];
}else if(_size!=expr._get_size()){
  delete [] _value;
  _size=expr._get_size();
  _value=new T[_size];
}


for(int i=0;i<_size;++i){
   _value[i]=expr._eval(i);
}
return(*this);
}

template <class L,class Op> inline CVector<T>& operator=(const CVectorExprL<L,Op>& expr){
if(!_value){
  _size=expr._get_size();
  _value=new T[_size];
}else if(_size!=expr._get_size()){
  delete [] _value;
  _size=expr._get_size();
  _value=new T[_size];
}


for(int i=0;i<_size;++i){
   _value[i]=expr._eval(i);
}
return(*this);
}

inline CVector<T>& operator=(const CVector<T>& expr){

if(!_value){
  _size=expr._get_size();
  _value=new T[_size];
}else if(_size!=expr._get_size()){
  delete [] _value;
  _size=expr._get_size();
  _value=new T[_size];
}

for(int i=0;i<_size;++i){
   _value[i]=expr._eval(i);
}
return(*this);
}

inline CVector<T>& operator<<(int s){
for(int i=0;i<_size;++i){
   _value[i]<<s;
}
return(*this);
}


inline CVector<T>& operator>>(int s){
for(int i=0;i<_size;++i){
   _value[i]>>s;
}
return(*this);
}

template <class L,class Op,class R> inline CVector<T>& operator+=(const CVectorExpr<L,Op,R>& expr){
for(int i=0;i<_size;++i){
   _value[i]+=expr._eval(i);
}
return(*this);
}

template <class Op,class R> inline CVector<T>* operator+=(const CVectorExprR<Op,R>& expr){
for(int i=0;i<_size;++i){
   _value[i]+=expr._eval(i);
}
return(this);
}

template <class L,class Op> inline CVector<T>* operator+=(const CVectorExprL<L,Op>& expr){
for(int i=0;i<_size;++i){
   _value[i]+=expr._eval(i);
}
return(this);
}

inline CVector<T>* operator+=(const CVector<T>& v){ 
for(int i=0;i<_size;++i){
   _value[i]+=v[i];
}
return(this);
}

inline CVector<T>* operator+=(const T* v){ 
for(int i=0;i<_size;++i){
   _value[i]+=v[i];
}
return(this);
}

inline CVector<T>* operator+=(const T v){ 
for(int i=0;i<_size;++i){
   _value[i]+=v;
}
return(this);
}

template <class L,class Op,class R> inline CVector<T>* operator-=(const CVectorExpr<L,Op,R>& expr){ 
for(int i=0;i<_size;++i){
   _value[i]-=expr._eval(i);
}
return(this);
}

template <class Op,class R> inline CVector<T>* operator-=(const CVectorExprR<Op,R>& expr){
for(int i=0;i<_size;++i){
   _value[i]-=expr._eval(i);
}
return(this);
}

template <class L,class Op> inline CVector<T>* operator-=(const CVectorExprL<L,Op>& expr){
for(int i=0;i<_size;++i){
   _value[i]-=expr._eval(i);
}
return(this);
}

inline CVector<T>* operator-=(const CVector<T>& v){ 
for(int i=0;i<_size;++i){
   _value[i]-=v[i];
}
return(this);
}

inline CVector<T>* operator-=(const T* v){ 
for(int i=0;i<_size;++i){
   _value[i]-=v[i];
}
return(this);
}

inline CVector<T>* operator-=(const T v){ 
for(int i=0;i<_size;++i){
   _value[i]-=v;
}
return(this);
}

template <class L,class Op,class R> inline CVector<T>* operator*=(const CVectorExpr<L,Op,R>& expr){ 
for(int i=0;i<_size;++i){
   _value[i]*=expr._eval(i);
}
return(this);
}

template <class Op,class R> inline CVector<T>* operator*=(const CVectorExprR<Op,R>& expr){
for(int i=0;i<_size;++i){
   _value[i]*=expr._eval(i);
}
return(this);
}

template <class L,class Op> inline CVector<T>* operator*=(const CVectorExprL<L,Op>& expr){
for(int i=0;i<_size;++i){
   _value[i]*=expr._eval(i);
}
return(this);
}

inline CVector<T>* operator*=(const CVector<T>& v){ 
for(int i=0;i<_size;++i){
   _value[i]*=v[i];
}
return(this);
}

inline CVector<T>* operator*=(const T* v){ 
for(int i=0;i<_size;++i){
   _value[i]*=v[i];
}
return(this);
}

inline CVector<T>* operator*=(const T v){ 
for(int i=0;i<_size;++i){
   _value[i]*=v;
}
return(this);
}

template <class L,class Op,class R> inline CVector<T>* operator/=(const CVectorExpr<L,Op,R>& expr){ 
for(int i=0;i<_size;++i){
   _value[i]/=expr._eval(i);
}
return(this);
}

template <class Op,class R> inline CVector<T>* operator/=(const CVectorExprR<Op,R>& expr){
for(int i=0;i<_size;++i){
   _value[i]/=expr._eval(i);
}
return(this);
}

template <class L,class Op> inline CVector<T>* operator/=(const CVectorExprL<L,Op>& expr){
for(int i=0;i<_size;++i){
   _value[i]/=expr._eval(i);
}
return(this);
}

inline CVector<T>* operator/=(const CVector<T>& v){ 
for(int i=0;i<_size;++i){
   _value[i]/=v[i];
}
return(this);
}

inline CVector<T>* operator/=(const T* v){ 
for(int i=0;i<_size;++i){
   _value[i]/=v[i];
}
return(this);
}

inline CVector<T>* operator/=(const T v){ 
for(int i=0;i<_size;++i){
   _value[i]/=v;
}
return(this);
}

inline void _zero(){
for(int i=0;i<_size;++i){
   _value[i]=0;
}
}

void _resize(int size){
T* tmp=new T[size];
memset(tmp,0,sizeof(T));
if(_value){
  int s=std::min(size,_size);
  for(int i=0;i<s;++i) tmp[i]=_value[i];
  delete [] _value;
}
_size=size;
_value=tmp;
}

void _swap(CVector<T>& v){
T* tmp=_value;
_value=v._value;
v._value=tmp;
}

void _swap(T* v){
T* tmp=_value;
_value=v;
v=tmp;
}

void _remove(const unsigned int start,const unsigned int end){
int j=0;
for(int i=end+1;i<_size;++i){
   _value[start+j]=_value[i];
   ++j;
}
_size=start+j;
}

inline T _max(){
T m=_value[0];
for(int i=1;i<_size;++i){
   if(_value[i]>m) m=_value[i];
}
return(m);
}

inline T _min(){
T m=_value[0];
for(int i=1;i<_size;++i){
   if(_value[i]<m) m=_value[i];
}
return(m);
}


inline T _sum(){
T s=0;
for(int i=0;i<_size;++i) s+=_value[i];
return(s);
}

inline void _sum(T& s){
s=0;
for(int i=0;i<_size;++i) s+=_value[i];
}

inline T _norm(){
T d=0;
for(int i=0;i<_size;++i) d+=_value[i]*_value[i];
return(sqrt(d));
}

inline void _normalize(){
T d=0;
for(int i=0;i<_size;++i) d+=_value[i]*_value[i];
d=1.0/sqrt(d);
for(int i=0;i<_size;++i) _value[i]*=d;
}

inline T _dot(CVector<T>& v){
T d=0;
for(int i=0;i<_size;++i) d+=v[i]*_value[i];
return(d);
}

void _create_random_vector(){
for(int i=0;i<_size;++i) _value[i]=rand();
double s=0;
for(int i=0;i<_size;++i) s+=_value[i];
s=1.0/s;
for(int i=0;i<_size;++i) _value[i]*=s;
}

void _gram_schmidt(CVector<T>& v){
T d=_dot(v);
for(int i=0;i<_size;++i){
   v[i]-=d*_value[i];
}
}

};


template <class L,class R> inline double dot(CVectorContainer<L,typename L::value_type>& v1,CVectorContainer<R,typename R::value_type>& v2){
double d=0;
int s=v1._get_size();
for(int i=0;i<s;++i) d+=v1._eval(i)*v2._eval(i);
return(d);
}

template <class T> inline T dist(CVector<T>& v1,CVector<T>& v2){
T d=0;
int s=v1._get_size();
for(int i=0;i<s;++i) d+=(v1._eval(i)-v2._eval(i))*(v1._eval(i)-v2._eval(i));
return(sqrt(d));
}

template <class T> inline T distSquare(CVector<T>& v1,CVector<T>& v2){
T d=0;
int s=v1._get_size();
for(int i=0;i<s;++i) d+=(v1._eval(i)-v2._eval(i))*(v1._eval(i)-v2._eval(i));
return(d);
}


#endif
