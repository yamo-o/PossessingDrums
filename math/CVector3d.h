

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
    for more information
    
    contact: Kazuhiko Yamamoto [yamotulp(at)gmail.com] 
             Twitter: @yamo_o
    
==============================================================================*/



#include "CVector2d.h"

#ifndef CVEC3D_HH
#define CVEC3D_HH


template <class T,class S> class CVector3DContainer{
public:

inline const S _eval(int i) const {
return((static_cast<const T*>(this))->_eval(i));
}

};


template <class L,class Op,class R> class CVector3DExpr:public CVector3DContainer<CVector3DExpr<L,Op,R>,typename R::value_type>{
const CVector3DContainer<L,typename L::value_type> &_l;
const CVector3DContainer<R,typename R::value_type> &_r;
public:

typedef typename R::value_type value_type;

CVector3DExpr(const CVector3DContainer<L,typename L::value_type>& l,const CVector3DContainer<R,typename R::value_type>& r):_l(l),_r(r){}

inline const value_type _eval(int i) const {
return(Op::_apply(_l._eval(i),_r._eval(i)));
}
  
};

template <class Op,class R> class CVector3DExprR:public CVector3DContainer<CVector3DExprR<Op,R>,typename R::value_type>{
const double &_l;
const CVector3DContainer<R,typename R::value_type> &_r;
public:

typedef typename R::value_type value_type;

CVector3DExprR(const double& l,const CVector3DContainer<R,typename R::value_type>& r):_l(l),_r(r){}

inline const value_type _eval(int i) const {
return(Op::_apply(_l,_r._eval(i)));
}

};

template <class L,class Op> class CVector3DExprL:public CVector3DContainer<CVector3DExprL<L,Op>,typename L::value_type>{
const CVector3DContainer<L,typename L::value_type> &_l;
const double &_r;
public:

typedef typename L::value_type value_type;

CVector3DExprL(const CVector3DContainer<L,typename L::value_type>& l,const double& r):_l(l),_r(r){}

inline const value_type _eval(int i) const {
return(Op::_apply(_l._eval(i),_r));
}
  
};


// vector + vector
template <class L,class R> inline CVector3DExpr<L,OpAdd<typename L::value_type,typename L::value_type,typename R::value_type>,R> operator+(const CVector3DContainer<L,typename L::value_type>& lhs,const CVector3DContainer<R,typename R::value_type>& rhs){
return(CVector3DExpr<L,OpAdd<typename L::value_type,typename L::value_type,typename R::value_type>,R>(lhs,rhs));
}

// scholar + vector
template <class R> inline CVector3DExprR<OpAdd<typename R::value_type,double,typename R::value_type>,R> operator+(const double& lhs,const CVector3DContainer<R,typename R::value_type>& rhs){
return(CVector3DExprR<OpAdd<typename R::value_type,double,typename R::value_type>,R>(lhs,rhs));
}

// vector + scholar
template <class L> inline CVector3DExprL<L,OpAdd<typename L::value_type,typename L::value_type,double> > operator+(const CVector3DContainer<L,typename L::value_type>& lhs,const double& rhs){
return(CVector3DExprL<L,OpAdd<typename L::value_type,typename L::value_type,double> >(lhs,rhs));
}

// vector - vector
template <class L,class R> inline CVector3DExpr<L,OpSub<typename L::value_type,typename L::value_type,typename R::value_type>,R> operator-(const CVector3DContainer<L,typename L::value_type>& lhs,const CVector3DContainer<R,typename R::value_type>& rhs){
return(CVector3DExpr<L,OpSub<typename L::value_type,typename L::value_type,typename R::value_type>,R>(lhs,rhs));
}

// scholar - vector
template <class R> inline CVector3DExprR<OpSub<typename R::value_type,double,typename R::value_type>,R> operator-(const double& lhs,const CVector3DContainer<R,typename R::value_type>& rhs){
return(CVector3DExprR<OpSub<typename R::value_type,double,typename R::value_type>,R>(lhs,rhs));
}

// vector - scholar
template <class L> inline CVector3DExprL<L,OpSub<typename L::value_type,typename L::value_type,double> > operator-(const CVector3DContainer<L,typename L::value_type>& lhs,const double& rhs){
return(CVector3DExprL<L,OpSub<typename L::value_type,typename L::value_type,double> >(lhs,rhs));
}

// vector * vector
template <class L,class R> inline CVector3DExpr<L,OpMul<typename L::value_type,typename L::value_type,typename R::value_type>,R> operator*(const CVector3DContainer<L,typename L::value_type>& lhs,const CVector3DContainer<R,typename R::value_type>& rhs){
return(CVector3DExpr<L,OpMul<typename L::value_type,typename L::value_type,typename R::value_type>,R>(lhs,rhs));
}

// scholar * vector
template <class R> inline CVector3DExprR<OpMul<typename R::value_type,double,typename R::value_type>,R> operator*(const double& lhs,const CVector3DContainer<R,typename R::value_type>& rhs){
return(CVector3DExprR<OpMul<typename R::value_type,double,typename R::value_type>,R>(lhs,rhs));
}

// vector * scholar
template <class L> inline CVector3DExprL<L,OpMul<typename L::value_type,typename L::value_type,double> > operator*(const CVector3DContainer<L,typename L::value_type>& lhs,const double& rhs){
return(CVector3DExprL<L,OpMul<typename L::value_type,typename L::value_type,double> >(lhs,rhs));
}

// vector / vector
template <class L,class R> inline CVector3DExpr<L,OpDiv<typename L::value_type,typename L::value_type,typename R::value_type>,R> operator/(const CVector3DContainer<L,typename L::value_type>& lhs,const CVector3DContainer<R,typename R::value_type>& rhs){
return(CVector3DExpr<L,OpDiv<typename L::value_type,typename L::value_type,typename R::value_type>,R>(lhs,rhs));
}

// scholar / vector
template <class R> inline CVector3DExprR<OpDiv<typename R::value_type,double,typename R::value_type>,R> operator/(const double& lhs,const CVector3DContainer<R,typename R::value_type>& rhs){
return(CVector3DExprR<OpDiv<typename R::value_type,double,typename R::value_type>,R>(lhs,rhs));
}

// vector / scholar
template <class L> inline CVector3DExprL<L,OpDiv<typename L::value_type,typename L::value_type,double> > operator/(const CVector3DContainer<L,typename L::value_type>& lhs,const double& rhs){
return(CVector3DExprL<L,OpDiv<typename L::value_type,typename L::value_type,double> >(lhs,rhs));
}



template <class T> class CVector3d:public CVector3DContainer<CVector3d<T>,T>{
public:

typedef T value_type;
T _value[3];

CVector3d(){}

CVector3d(T x,T y,T z){
_value[0]=x;
_value[1]=y;
_value[2]=z;
}

CVector3d(CVector3d<T>& v){
_value[0]=v[0];
_value[1]=v[1];
_value[2]=v[2];
}

CVector3d(CVector2d<T>& v){
_value[0]=v[0];
_value[1]=v[1];
_value[2]=0;
}

explicit CVector3d(T v){
_value[0]=v;
_value[1]=v;
_value[2]=v;
}

explicit CVector3d(T* v){
_value[0]=v[0];
_value[1]=v[1];
_value[2]=v[2];
}

template <class S> CVector3d(const CVector3DContainer<S,T>& expr){
for(int i=0;i<3;++i){
   _value[i]=expr._eval(i);
}	
}

~CVector3d(){}

inline T& operator[](int i){ return(_value[i]); }
inline const T& operator[](int i) const { return(_value[i]); }
inline T& operator()(int i){ return(_value[i]); }
inline const T& operator()(int i) const { return(_value[i]); }
inline const T& _eval(int i) const { return(_value[i]); }

inline void _set(T x,T y,T z){
_value[0]=x; _value[1]=y; _value[2]=z;
}

inline void _set(const T* x){
_value[0]=x[0]; _value[1]=x[1]; _value[2]=x[2];
}

inline void _minSet(CVector3d<T>& v){
if(v[0]<_value[0]) _value[0]=v[0];
if(v[1]<_value[1]) _value[0]=v[1];
if(v[2]<_value[2]) _value[0]=v[2];
}

inline void _maxSet(CVector3d<T>& v){
if(v[0]>_value[0]) _value[0]=v[0];
if(v[1]>_value[1]) _value[0]=v[1];
if(v[2]>_value[2]) _value[0]=v[2];
}

inline CVector3d<T>& operator=(const double& scholar){
for(int i=0;i<3;++i){
   _value[i]=scholar;
}
return(*this);
}

inline CVector3d<T>& operator=(const CVector3d<T>& vec){
for(int i=0;i<3;++i){
   _value[i]=vec[i];
}
return(*this);
}

template <class L,class Op,class R> inline CVector3d<T>& operator=(const CVector3DExpr<L,Op,R>& expr){
for(int i=0;i<3;++i){
   _value[i]=expr._eval(i);
}
return(*this);
}

template <class Op,class R> inline CVector3d<T>& operator=(const CVector3DExprR<Op,R>& expr){
for(int i=0;i<3;++i){
   _value[i]=expr._eval(i);
}
return(*this);
}

template <class L,class Op> inline CVector3d<T>& operator=(const CVector3DExprL<L,Op>& expr){
for(int i=0;i<3;++i){
   _value[i]=expr._eval(i);
}
return(*this);
}

inline int operator==(const CVector3d<T>& vec) const {
return(_value[0]==vec[0] && _value[1]==vec[1] && _value[2]==vec[2]);
}

inline CVector3d<T>& operator<<(int s){
for(int i=0;i<3;++i){
   _value[i]<<s;
}
return(*this);
}


inline CVector3d<T>& operator>>(int s){
for(int i=0;i<3;++i){
   _value[i]>>s;
}
return(*this);
}

template <class L,class Op,class R> inline CVector3d<T>* operator+=(const CVector3DExpr<L,Op,R>& expr){ 
for(int i=0;i<3;++i){
   _value[i]+=expr._eval(i);
}
return(this);
}

template <class Op,class R> inline CVector3d<T>* operator+=(const CVector3DExprR<Op,R>& expr){
for(int i=0;i<3;++i){
   _value[i]+=expr._eval(i);
}
return(this);
}

template <class L,class Op> inline CVector3d<T>* operator+=(const CVector3DExprL<L,Op>& expr){
for(int i=0;i<3;++i){
   _value[i]+=expr._eval(i);
}
return(this);
}

inline CVector3d<T>* operator+=(const CVector3d<T>& v){ 
for(int i=0;i<3;++i){
   _value[i]+=v[i];
}
return(this);
}

inline CVector3d<T>* operator+=(const T* v){ 
for(int i=0;i<3;++i){
   _value[i]+=v[i];
}
return(this);
}

inline CVector3d<T>* operator+=(const T v){ 
for(int i=0;i<3;++i){
   _value[i]+=v;
}
return(this);
}

template <class L,class Op,class R> inline CVector3d<T>* operator-=(const CVector3DExpr<L,Op,R>& expr){ 
for(int i=0;i<3;++i){
   _value[i]-=expr._eval(i);
}
return(this);
}

template <class Op,class R> inline CVector3d<T>* operator-=(const CVector3DExprR<Op,R>& expr){
for(int i=0;i<3;++i){
   _value[i]-=expr._eval(i);
}
return(this);
}

template <class L,class Op> inline CVector3d<T>* operator-=(const CVector3DExprL<L,Op>& expr){
for(int i=0;i<3;++i){
   _value[i]-=expr._eval(i);
}
return(this);
}

inline CVector3d<T>* operator-=(const CVector3d<T>& v){ 
for(int i=0;i<3;++i){
   _value[i]-=v[i];
}
return(this);
}

inline CVector3d<T>* operator-=(const T* v){ 
for(int i=0;i<3;++i){
   _value[i]-=v[i];
}
return(this);
}

inline CVector3d<T>* operator-=(const T v){ 
for(int i=0;i<3;++i){
   _value[i]-=v;
}
return(this);
}

template <class L,class Op,class R> inline CVector3d<T>* operator*=(const CVector3DExpr<L,Op,R>& expr){ 
for(int i=0;i<3;++i){
   _value[i]*=expr._eval(i);
}
return(this);
}

template <class Op,class R> inline CVector3d<T>* operator*=(const CVector3DExprR<Op,R>& expr){
for(int i=0;i<3;++i){
   _value[i]*=expr._eval(i);
}
return(this);
}

template <class L,class Op> inline CVector3d<T>* operator*=(const CVector3DExprL<L,Op>& expr){
for(int i=0;i<3;++i){
   _value[i]*=expr._eval(i);
}
return(this);
}

inline CVector3d<T>* operator*=(const CVector3d<T>& v){ 
for(int i=0;i<3;++i){
   _value[i]*=v[i];
}
return(this);
}

inline CVector3d<T>* operator*=(const T* v){ 
for(int i=0;i<3;++i){
   _value[i]*=v[i];
}
return(this);
}

inline CVector3d<T>* operator*=(const T v){ 
for(int i=0;i<3;++i){
   _value[i]*=v;
}
return(this);
}

template <class L,class Op,class R> inline CVector3d<T>* operator/=(const CVector3DExpr<L,Op,R>& expr){ 
for(int i=0;i<3;++i){
   _value[i]/=expr._eval(i);
}
return(this);
}

template <class Op,class R> inline CVector3d<T>* operator/=(const CVector3DExprR<Op,R>& expr){
for(int i=0;i<3;++i){
   _value[i]/=expr._eval(i);
}
return(this);
}

template <class L,class Op> inline CVector3d<T>* operator/=(const CVector3DExprL<L,Op>& expr){
for(int i=0;i<3;++i){
   _value[i]/=expr._eval(i);
}
return(this);
}

inline CVector3d<T>* operator/=(const CVector3d<T>& v){ 
for(int i=0;i<3;++i){
   _value[i]/=v[i];
}
return(this);
}

inline CVector3d<T>* operator/=(const T* v){ 
for(int i=0;i<3;++i){
   _value[i]/=v[i];
}
return(this);
}

inline CVector3d<T>* operator/=(const T v){ 
for(int i=0;i<3;++i){
   _value[i]/=v;
}
return(this);
}

inline void _zero(){
for(int i=0;i<3;++i){
   _value[i]=0;
}
}

inline T _norm(){
T d=0;
for(int i=0;i<3;++i) d+=_value[i]*_value[i];
return(sqrt(d));
}

inline T _squareNorm(){
T d=0;
for(int i=0;i<3;++i) d+=_value[i]*_value[i];
return(d);
}

inline void _normalize(){
T d=0;
for(int i=0;i<3;++i) d+=_value[i]*_value[i];
d=1.0/sqrt(d);
for(int i=0;i<3;++i) _value[i]*=d;
}

inline T _dot(CVector3d<T>& v){
T d=0;
for(int i=0;i<3;++i) d+=v[i]*_value[i];
return(d);
}

void _create_random_vector(){
for(int i=0;i<3;++i) _value[i]=rand();
double s=0;
for(int i=0;i<3;++i) s+=_value[i];
s=1.0/s;
for(int i=0;i<3;++i) _value[i]*=s;
}

void _gram_schmidt(CVector3d<T>& v){
T d=_dot(v);
for(int i=0;i<3;++i){
   v[i]-=d*_value[i];
}
}

};



template <class L,class R> inline typename L::value_type dist(CVector3DContainer<L,typename L::value_type>& v1,CVector3DContainer<R,typename R::value_type>& v2){
return(sqrt(square(v1._eval(0)-v2._eval(0))+square(v1._eval(1)-v2._eval(1))+square(v1._eval(2)-v2._eval(2))));
}

template <class L,class R> inline typename L::value_type distSquare(CVector3DContainer<L,typename L::value_type>& v1,CVector3DContainer<R,typename R::value_type>& v2){
return(square(v1._eval(0)-v2._eval(0))+square(v1._eval(1)-v2._eval(1))+square(v1._eval(2)-v2._eval(2)));
}

template <class L,class R> inline typename L::value_type dot(CVector3DContainer<L,typename L::value_type>& v1,CVector3DContainer<R,typename R::value_type>& v2){
typename L::value_type d=0;
for(int i=0;i<3;++i) d+=v1._eval(i)*v2._eval(i);
return(d);
}

template <class L,class R> inline void cross(CVector3DContainer<L,typename L::value_type>& v1,CVector3DContainer<R,typename R::value_type>& v2,CVector3d<typename L::value_type>& res){
CVector3d<typename L::value_type> vec1(v1._eval(0),v1._eval(1),v1._eval(2));
CVector3d<typename L::value_type> vec2(v2._eval(0),v2._eval(1),v2._eval(2));
res[0]=vec1[1]*vec2[2]-vec1[2]*vec2[1];
res[1]=vec1[2]*vec2[0]-vec1[0]*vec2[2];
res[2]=vec1[0]*vec2[1]-vec1[1]*vec2[0];
}

template <class L,class R> inline void cross_add(CVector3DContainer<L,typename L::value_type>& v1,CVector3DContainer<R,typename R::value_type>& v2,CVector3d<typename L::value_type>& res){
CVector3d<typename L::value_type> vec1(v1._eval(0),v1._eval(1),v1._eval(2));
CVector3d<typename L::value_type> vec2(v2._eval(0),v2._eval(1),v2._eval(2));
res[0]+=vec1[1]*vec2[2]-vec1[2]*vec2[1];
res[1]+=vec1[2]*vec2[0]-vec1[0]*vec2[2];
res[2]+=vec1[0]*vec2[1]-vec1[1]*vec2[0];
}




#endif

