

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


#ifndef TYPE_HH
#define TYPE_HH
#pragma once


#include <map>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <list>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <time.h>
#include <pthread.h>
#include <queue>
#include <memory.h>
#include <float.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#define M_TWOPI 2.0*3.14159265358979323846264338327950288
#define EPSILON 1.0e-6
#define MIN(a,b) (a<b?a:b)
#define MAX(a,b) (a>b?a:b)
#define SIGN(a,b) ((b)>0.0?fabs(a):-fabs(a))


template <class R,class T,class S> struct OpAdd{
static inline R _apply(T i,S j){
return(i+j);
}
};

template <class R,class T,class S> struct OpSub{
static inline R _apply(T i,S j){
return(i-j);
}
};

template <class R,class T,class S> struct OpMul{
static inline R _apply(T i,S j){
return(i*j);
}
};

template <class R,class T,class S> struct OpDiv{
static inline R _apply(T i,S j){
return(i/j);
}
};

template <class T> inline T square(T s){ return(s*s); }
template <class T> inline T square(T s1,T s2){ return(s1*s1+s2*s2); }
inline float random(float minf,float maxf){
return(minf+(maxf-minf)*(rand()/(float)INT_MAX));
}


#endif
