

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


#include "CModalRespnator.h"

#ifndef CSCMOD_HH
#define CSCMOD_HH

#pragma once 



template <class T> class CScratchModulation:public CModulation<T>{
public:

T _gammma;
T _freqFract;
T _time;
T _intVel;
CQueue<T> _vel;
T _prevVel;
T _scale;
T _threthold;
int _sampleRate;
T _attackTime;

CRelatedBumpNoise<T> _noise;

void _init(T attackTime,T freqFract,unsigned int sampleRate,T threthold){
_sampleRate=sampleRate;
_attackTime=attackTime;
_freqFract=freqFract;
_vel._resize((int)(attackTime*sampleRate*2));
_noise._init(1,sampleRate,20,0.1);
_threthold=sampleRate*attackTime*threthold;
_scale=sampleRate*attackTime*0.25;
_reset();
}

void _reset(){
_time=0;
_intVel=0;
_prevVel=0;
_gammma=0;
}

//vel must be scaled from 0.0 to 1.0
void _setVelocity(T vel){
_vel._enqueue(vel);
}


virtual inline T _render(int param,T* input,void* ref){

CModalResonator<T>* reson=(CModalResonator<T>*)ref;

T v;
if(_vel._empty()) v=0;
else v=_vel._dequeue();

_intVel+=0.5*(_prevVel+v);
if(_intVel>=_attackTime*_sampleRate) _intVel=_attackTime*_sampleRate;
_prevVel=v;

if(_intVel>_threthold){
  _gammma=(_intVel-_threthold)/(_attackTime*_sampleRate-_threthold);
  if(_gammma>1.0) _gammma=1.0;
}else{
  if(_gammma<=0) return(0);
  _gammma=0;
}

const T freq=(1+_freqFract*(_noise._render(_gammma)+0.05)*_gammma)*reson->_freq[param];
const T invWavenumber=M_TWOPI/freq;
const T amp=invWavenumber*sin(M_TWOPI*_intVel/_scale);

const T ds=2.0*reson->_contactDamping*reson->_auxDampScale/(T)reson->_sampleRate;
const T fs=reson->_auxFreqScale/(T)reson->_sampleRate;
T t1=exp(-reson->_damp[param]*ds);
T t2=M_TWOPI*freq*fs;
T creal=cos(t2)*t1;
T cimg=sin(t2)*t1;
T t3=sqrt(1.0-cimg*cimg);
reson->_cplus[param]=creal+t3;
reson->_cminus[param]=creal-t3;
reson->_aa[param]=reson->_amp[param]*cimg*reson->_auxAmpScale; 

_time+=1.0/(T)reson->_sampleRate;
_intVel*=exp(-reson->_damp[param]*1.0/(T)reson->_sampleRate);

*input*=amp;

return(0);
}

};



#endif