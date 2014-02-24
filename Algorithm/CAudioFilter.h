

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

#ifndef CAUDIOFILTER_HH
#define CAUDIOFILTER_HH

#include "FFT.h"
#include "CRingBufferMatrix.h"


template <class T> class COnePoleFilter{
public:

T _a,_b;
T _gain;
T _tmp[2];
T _lastFrame;

COnePoleFilter(){
_tmp[0]=_tmp[1]=0.0;
_gain=1.0;
_setPole(0.9);
}

void _setPole(T pole){
if(pole>0.0) _b=1.0-pole;
else _b=1.0+pole;
_a=-pole;
}

inline T _render(T input){
_tmp[0]=_gain*input;
_lastFrame=_b*_tmp[0]-_a*_tmp[1];
_tmp[1]=_lastFrame;
return(_lastFrame);
}

virtual void _render(int size,T* buf){
for(int i=0;i<size;i++){
   _tmp[0]=_gain*buf[i];
   buf[i]=_b*_tmp[0]-_a*_tmp[1];
   _tmp[1]=buf[i];
}
}

void _clear(){
_tmp[0]=0;
_tmp[1]=0;
_lastFrame=0; 
}

};


template <class T> class CLowPassFilter{
public:

T _b,_y;
T _gain;
T _newGain;

CLowPassFilter(){
_b=_y=0.0;
_gain=1.0;
_newGain=1.0;
}

void _setCutoff(T freq,int sampleRate){
T f=freq*2.0*M_PI/(T)sampleRate;
if(f>M_PI) f=M_PI;
T t=2.0-cos(f);
_b=1.0-(1.0-(t-sqrt(t*t-1.0)))*1.2071;
}


T _render(T input){
T g=_gain*(1.0-_b);
T gIncr=(_newGain-_gain)*(1.0-_b);
T buf=g*input+_b*_y;
_y=buf;
g+=gIncr;
_gain=_newGain;

return(buf);
}

void _render(int size,T* buf){
T g=_gain*(1.0-_b);
T gIncr=(_newGain-_gain)/(T)size*(1.0-_b);

for(int i=0;i<size;i++){
   buf[i]=g*buf[i]+_b*_y;
   _y=buf[i];
   g+=gIncr;
}

_gain=_newGain;
}

};


template <class T> class CAudioModeFilter{

CVector<T> _a,_b;
CVector<T> _inputs,_outputs;
int _num_modes;
T _gain;

public:

CAudioModeFilter(){
_gain=0.6;
}

~CAudioModeFilter(){}

void _constructInverseFilter(int num_modes,T* freq,T* gain,int sampleRate){
_num_modes=num_modes;
_a._resize(3*num_modes);
_b._resize(3*num_modes);
_a=0;
_b=0;
_inputs._resize(3*num_modes);
_inputs=0;
_outputs._resize(3*num_modes);
_outputs=0;

for(int i=0;i<num_modes;i++){
   _setNotch(i,freq[i],gain[i],&_b[0],sampleRate);
   _a[3*i]=1.0;
   _b[3*i]=1.0;
}

}

void _setResonance(T frequency,T radius,bool normalize,T* a,T* b,int sampleRate){
a[1]=radius*radius;
a[0]=-2.0*radius*cos(2.0*M_PI*frequency/(T)sampleRate);
if(normalize){
  b[0]=0.5-0.5*a[2];
  b[1]=0.0;
  b[2]=-b[0];
}
}

void _setNotch(int mode,T frequency,T radius,T* b,int sampleRate){
b[3*mode+2]=radius*radius;
b[3*mode+1]=-2.0*radius*cos(2.0*M_PI*frequency/(T)sampleRate);
}

virtual void _render(int size,T* buf){
for(int m=0;m<_num_modes;m++){
   for(int i=0;i<size;i++){
      _inputs[0]=_gain*buf[i];
      buf[i]=_b[3*m]*_inputs[0]+_b[3*m+1]*_inputs[1]+_b[3*m+2]*_inputs[2];
      buf[i]-=_a[3*m+2]*_outputs[1]+_a[3*m+1]*_outputs[1];
      _inputs[2]=_inputs[1];
      _inputs[1]=_inputs[0];
      _outputs[2]=_outputs[1];
      _outputs[1]=buf[i];
   }
}
}

};


#endif