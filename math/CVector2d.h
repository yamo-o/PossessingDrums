

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


#include <ostream>
#include "CTypes.h"

#ifndef CVEC2D_HH
#define CVEC2D_HH


template <class T,class S> class CVector2DContainer{
public:

inline const S _eval(int i) const {
return((static_cast<const T*>(this))->_eval(i));
}

};


template <class L,class Op,class R> class CVector2DExpr:public CVector2DContainer<CVector2DExpr<L,Op,R>,typename R::value_type>{
const CVector2DContainer<L,typename L::value_type> &_l;
const CVector2DContainer<R,typename R::value_type> &_r;
public:

typedef typename R::value_type value_type;

CVector2DExpr(const CVector2DContainer<L,typename L::value_type>& l,const CVector2DContainer<R,typename R::value_type>& r):_l(l),_r(r){}

inline const value_type _eval(int i) const {
return(Op::_apply(_l._eval(i),_r._eval(i)));
}
  
};

template <class Op,class R> class CVector2DExprR:public CVector2DContainer<CVector2DExprR<Op,R>,typename R::value_type>{
const double &_l;
const CVector2DContainer<R,typename R::value_type> &_r;
public:

typedef typename R::value_type value_type;

CVector2DExprR(const double& l,const CVector2DContainer<R,typename R::value_type>& r):_l(l),_r(r){}

inline const value_type _eval(int i) const {
return(Op::_apply(_l,_r._eval(i)));
}

};

template <class L,class Op> class CVector2DExprL:public CVector2DContainer<CVector2DExprL<L,Op>,typename L::value_type>{
const CVector2DContainer<L,typename L::value_type> &_l;
const double &_r;
public:

typedef typename L::value_type value_type;

CVector2DExprL(const CVector2DContainer<L,typename L::value_type>& l,const double& r):_l(l),_r(r){}

inline const value_type _eval(int i) const {
return(Op::_apply(_l._eval(i),_r));
}
  
};


// vector + vector
template <class L,class R> inline CVector2DExpr<L,OpAdd<typename L::value_type,typename L::value_type,typename R::value_type>,R> operator+(const CVector2DContainer<L,typename L::value_type>& lhs,const CVector2DContainer<R,typename R::value_type>& rhs){
return(CVector2DExpr<L,OpAdd<typename L::value_type,typename L::value_type,typename R::value_type>,R>(lhs,rhs));
}

// scholar + vector
template <class R> inline CVector2DExprR<OpAdd<typename R::value_type,double,typename R::value_type>,R> operator+(const double& lhs,const CVector2DContainer<R,typename R::value_type>& rhs){
return(CVector2DExprR<OpAdd<typename R::value_type,double,typename R::value_type>,R>(lhs,rhs));
}

// vector + scholar
template <class L> inline CVector2DExprL<L,OpAdd<typename L::value_type,typename L::value_type,double> > operator+(const CVector2DContainer<L,typename L::value_type>& lhs,const double& rhs){
return(CVector2DExprL<L,OpAdd<typename L::value_type,typename L::value_type,double> >(lhs,rhs));
}

// vector - vector
template <class L,class R> inline CVector2DExpr<L,OpSub<typename L::value_type,typename L::value_type,typename R::value_type>,R> operator-(const CVector2DContainer<L,typename L::value_type>& lhs,const CVector2DContainer<R,typename R::value_type>& rhs){
return(CVector2DExpr<L,OpSub<typename L::value_type,typename L::value_type,typename R::value_type>,R>(lhs,rhs));
}

// scholar - vector
template <class R> inline CVector2DExprR<OpSub<typename R::value_type,double,typename R::value_type>,R> operator-(const double& lhs,const CVector2DContainer<R,typename R::value_type>& rhs){
return(CVector2DExprR<OpSub<typename R::value_type,double,typename R::value_type>,R>(lhs,rhs));
}

// vector - scholar
template <class L> inline CVector2DExprL<L,OpSub<typename L::value_type,typename L::value_type,double> > operator-(const CVector2DContainer<L,typename L::value_type>& lhs,const double& rhs){
return(CVector2DExprL<L,OpSub<typename L::value_type,typename L::value_type,double> >(lhs,rhs));
}

// vector * vector
template <class L,class R> inline CVector2DExpr<L,OpMul<typename L::value_type,typename L::value_type,typename R::value_type>,R> operator*(const CVector2DContainer<L,typename L::value_type>& lhs,const CVector2DContainer<R,typename R::value_type>& rhs){
return(CVector2DExpr<L,OpMul<typename L::value_type,typename L::value_type,typename R::value_type>,R>(lhs,rhs));
}

// scholar * vector
template <class R> inline CVector2DExprR<OpMul<typename R::value_type,double,typename R::value_type>,R> operator*(const double& lhs,const CVector2DContainer<R,typename R::value_type>& rhs){
return(CVector2DExprR<OpMul<typename R::value_type,double,typename R::value_type>,R>(lhs,rhs));
}

// vector * scholar
template <class L> inline CVector2DExprL<L,OpMul<typename L::value_type,typename L::value_type,double> > operator*(const CVector2DContainer<L,typename L::value_type>& lhs,const double& rhs){
return(CVector2DExprL<L,OpMul<typename L::value_type,typename L::value_type,double> >(lhs,rhs));
}

// vector / vector
template <class L,class R> inline CVector2DExpr<L,OpDiv<typename L::value_type,typename L::value_type,typename R::value_type>,R> operator/(const CVector2DContainer<L,typename L::value_type>& lhs,const CVector2DContainer<R,typename R::value_type>& rhs){
return(CVector2DExpr<L,OpDiv<typename L::value_type,typename L::value_type,typename R::value_type>,R>(lhs,rhs));
}

// scholar / vector
template <class R> inline CVector2DExprR<OpDiv<typename R::value_type,double,typename R::value_type>,R> operator/(const double& lhs,const CVector2DContainer<R,typename R::value_type>& rhs){
return(CVector2DExprR<OpDiv<typename R::value_type,double,typename R::value_type>,R>(lhs,rhs));
}

// vector / scholar
template <class L> inline CVector2DExprL<L,OpDiv<typename L::value_type,typename L::value_type,double> > operator/(const CVector2DContainer<L,typename L::value_type>& lhs,const double& rhs){
return(CVector2DExprL<L,OpDiv<typename L::value_type,typename L::value_type,double> >(lhs,rhs));
}



template <class T> class CVector2d:public CVector2DContainer<CVector2d<T>,T>{
public:

typedef T value_type;
T _value[2];

CVector2d(){}

CVector2d(T x,T y){
_value[0]=x;
_value[1]=y;
}

CVector2d(CVector2d<T>& v){
_value[0]=v[0];
_value[1]=v[1];
}

explicit CVector2d(T v){
_value[0]=v;
_value[1]=v;
}


template <class S> CVector2d(const CVector2DContainer<S,T>& expr){
for(int i=0;i<2;++i){
   _value[i]=expr._eval(i);
}	
}

~CVector2d(){}

inline T& operator[](int i){ return(_value[i]); }
inline const T& operator[](int i) const { return(_value[i]); }
inline T& operator()(int i){ return(_value[i]); }
inline const T& operator()(int i) const { return(_value[i]); }
inline const T& _eval(int i) const { return(_value[i]); }


inline bool operator==(const CVector2d<T>& vec){
return((_value[0]==vec[0]) && (_value[1]==vec[1]));
}


inline CVector2d<T>& operator=(const T scholar){
for(int i=0;i<2;++i){
   _value[i]=scholar;
}
return(*this);
}


inline CVector2d<T>& operator=(const CVector2d<T>& vec){
for(int i=0;i<2;++i){
   _value[i]=vec[i];
}
return(*this);
}

template <class L,class Op,class R> inline CVector2d<T>& operator=(const CVector2DExpr<L,Op,R>& expr){
for(int i=0;i<2;++i){
   _value[i]=expr._eval(i);
}
return(*this);
}

template <class Op,class R> inline CVector2d<T>& operator=(const CVector2DExprR<Op,R>& expr){
for(int i=0;i<2;++i){
   _value[i]=expr._eval(i);
}
return(*this);
}

template <class L,class Op> inline CVector2d<T>& operator=(const CVector2DExprL<L,Op>& expr){
for(int i=0;i<2;++i){
   _value[i]=expr._eval(i);
}
return(*this);
}

inline CVector2d<T>& operator<<(int s){
for(int i=0;i<2;++i){
   _value[i]<<s;
}
return(*this);
}


inline CVector2d<T>& operator>>(int s){
for(int i=0;i<2;++i){
   _value[i]>>s;
}
return(*this);
}

template <class L,class Op,class R> inline CVector2d<T>* operator+=(const CVector2DExpr<L,Op,R>& expr){ 
for(int i=0;i<2;++i){
   _value[i]+=expr._eval(i);
}
return(this);
}

template <class Op,class R> inline CVector2d<T>* operator+=(const CVector2DExprR<Op,R>& expr){
for(int i=0;i<2;++i){
   _value[i]+=expr._eval(i);
}
return(this);
}

template <class L,class Op> inline CVector2d<T>* operator+=(const CVector2DExprL<L,Op>& expr){
for(int i=0;i<2;++i){
   _value[i]+=expr._eval(i);
}
return(this);
}

inline CVector2d<T>* operator+=(const CVector2d<T>& v){ 
for(int i=0;i<2;++i){
   _value[i]+=v[i];
}
return(this);
}

inline CVector2d<T>* operator+=(const T* v){ 
for(int i=0;i<2;++i){
   _value[i]+=v[i];
}
return(this);
}

inline CVector2d<T>* operator+=(const T v){ 
for(int i=0;i<2;++i){
   _value[i]+=v;
}
return(this);
}

template <class L,class Op,class R> inline CVector2d<T>* operator-=(const CVector2DExpr<L,Op,R>& expr){ 
for(int i=0;i<2;++i){
   _value[i]-=expr._eval(i);
}
return(this);
}

template <class Op,class R> inline CVector2d<T>* operator-=(const CVector2DExprR<Op,R>& expr){
for(int i=0;i<2;++i){
   _value[i]-=expr._eval(i);
}
return(this);
}

template <class L,class Op> inline CVector2d<T>* operator-=(const CVector2DExprL<L,Op>& expr){
for(int i=0;i<2;++i){
   _value[i]-=expr._eval(i);
}
return(this);
}

inline CVector2d<T>* operator-=(const CVector2d<T>& v){ 
for(int i=0;i<2;++i){
   _value[i]-=v[i];
}
return(this);
}

inline CVector2d<T>* operator-=(const T* v){ 
for(int i=0;i<2;++i){
   _value[i]-=v[i];
}
return(this);
}

inline CVector2d<T>* operator-=(const T v){ 
for(int i=0;i<2;++i){
   _value[i]-=v;
}
return(this);
}

template <class L,class Op,class R> inline CVector2d<T>* operator*=(const CVector2DExpr<L,Op,R>& expr){ 
for(int i=0;i<2;++i){
   _value[i]*=expr._eval(i);
}
return(this);
}

template <class Op,class R> inline CVector2d<T>* operator*=(const CVector2DExprR<Op,R>& expr){
for(int i=0;i<2;++i){
   _value[i]*=expr._eval(i);
}
return(this);
}

template <class L,class Op> inline CVector2d<T>* operator*=(const CVector2DExprL<L,Op>& expr){
for(int i=0;i<2;++i){
   _value[i]*=expr._eval(i);
}
return(this);
}

inline CVector2d<T>* operator*=(const CVector2d<T>& v){ 
for(int i=0;i<2;++i){
   _value[i]*=v[i];
}
return(this);
}

inline CVector2d<T>* operator*=(const T* v){ 
for(int i=0;i<2;++i){
   _value[i]*=v[i];
}
return(this);
}

inline CVector2d<T>* operator*=(const T v){ 
for(int i=0;i<2;++i){
   _value[i]*=v;
}
return(this);
}

template <class L,class Op,class R> inline CVector2d<T>* operator/=(const CVector2DExpr<L,Op,R>& expr){ 
for(int i=0;i<2;++i){
   _value[i]/=expr._eval(i);
}
return(this);
}

template <class Op,class R> inline CVector2d<T>* operator/=(const CVector2DExprR<Op,R>& expr){
for(int i=0;i<2;++i){
   _value[i]/=expr._eval(i);
}
return(this);
}

template <class L,class Op> inline CVector2d<T>* operator/=(const CVector2DExprL<L,Op>& expr){
for(int i=0;i<2;++i){
   _value[i]/=expr._eval(i);
}
return(this);
}

inline CVector2d<T>* operator/=(const CVector2d<T>& v){ 
for(int i=0;i<2;++i){
   _value[i]/=v[i];
}
return(this);
}

inline CVector2d<T>* operator/=(const T* v){ 
for(int i=0;i<2;++i){
   _value[i]/=v[i];
}
return(this);
}

inline CVector2d<T>* operator/=(const T v){ 
for(int i=0;i<2;++i){
   _value[i]/=v;
}
return(this);
}


inline bool operator!=(const CVector2d<T>& v){
return((this->_value[0]!=v[0]) || (this->_value[1]!=v[1]));
}


inline void _zero(){
for(int i=0;i<2;++i){
   _value[i]=0;
}
}

inline T _norm(){
T d=0;
for(int i=0;i<2;++i) d+=_value[i]*_value[i];
return(sqrt(d));
}

inline T _normSquare(){
T d=0;
for(int i=0;i<2;++i) d+=_value[i]*_value[i];
return(d);
}

inline void _normalize(){
T d=0;
for(int i=0;i<2;++i) d+=_value[i]*_value[i];
d=1.0/sqrt(d);
for(int i=0;i<2;++i) _value[i]*=d;
}

inline T _dot(CVector2d<T>& v){
T d=0;
for(int i=0;i<2;++i) d+=v[i]*_value[i];
return(d);
}

void _create_random_vector(){
for(int i=0;i<2;++i) _value[i]=rand();
double s=0;
for(int i=0;i<2;++i) s+=_value[i];
s=1.0/s;
for(int i=0;i<2;++i) _value[i]*=s;
}

void _gram_schmidt(CVector2d<T>& v){
T d=_dot(v);
for(int i=0;i<2;++i){
   v[i]-=d*_value[i];
}
}

inline void _conjugate(CVector2d<T>& c){
c[0]=_value[0];
c[1]=-_value[1];
}

};

template <class L,class R> inline double cross(CVector2DContainer<L,typename L::value_type>& v1,CVector2DContainer<R,typename R::value_type>& v2){
CVector2d<typename L::value_type> l,r;
l[0]=v1._eval(0); l[1]=v1._eval(1);
r[0]=v2._eval(0); r[1]=v2._eval(1);
return(l[0]*r[1]-l[1]*r[0]);
}

template <class L,class R> inline typename L::value_type dist(CVector2DContainer<L,typename L::value_type>& v1,CVector2DContainer<R,typename R::value_type>& v2){
return(sqrt(square(v1._eval(0)-v2._eval(0),v1._eval(1)-v2._eval(1))));
}

template <class L,class R> inline typename L::value_type dot(CVector2DContainer<L,typename L::value_type>& v1,CVector2DContainer<R,typename R::value_type>& v2){
typename L::value_type d=0;
for(int i=0;i<2;++i) d+=v1._eval(i)*v2._eval(i);
return(d);
}




#endif

