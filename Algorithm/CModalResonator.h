

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


#include "CModulation.h"

#ifndef CMODALRESONATOR_HH
#define CMODALRESONATOR_HH

#pragma once 



template <class T>class CModalResonator{

public:

int _sampleRate;

int _num_modes;
CVector<T> _freq;
CVector<T> _amp;
CVector<T> _damp;

CVector<T> _cplus;
CVector<T> _cminus;
CVector<T> _aa;
CVector<T> _u;	
CVector<T> _v;	
CVector<T> _in;

COnePoleFilter<T> _onepole;

T _contactDamping;
T _auxDampScale;
T _auxFreqScale;
T _auxAmpScale;
T _gain;

CVector<int> _transition;
CVector<T> _dcp;
CVector<T> _dcm;
CVector<T> _da;

bool _playing;


public:

int _nActiveModes;
CModulation<T>** _mod;


CModalResonator():_playing(false),_mod(NULL){}

~CModalResonator(){
if(_mod) delete [] _mod;
}


void _init(int bufsize,int sampleRate,T contactDamping,T auxDampScale=1,T auxFreqScale=1.0,T auxAmpScale=0.001,T gain=0.5){
_contactDamping=contactDamping;
_auxDampScale=auxDampScale;
_auxFreqScale=auxFreqScale;
_auxAmpScale=auxAmpScale;
_sampleRate=sampleRate;
_in._resize(bufsize);
_in=0;
_gain=gain;
}

void _setModes(int num_modes,T* amp,T*freq,T* damp){
_nActiveModes=num_modes;
_cplus._resize(_nActiveModes);
_cminus._resize(_nActiveModes);
_aa._resize(_nActiveModes);
_u._resize(_nActiveModes);
_v._resize(_nActiveModes);
_amp._resize(num_modes);
_freq._resize(num_modes);
_damp._resize(num_modes);
for(int i=0;i<num_modes;++i){
   _amp[i]=amp[i];
   _freq[i]=freq[i];
   _damp[i]=damp[i];
}
_transition._resize(_nActiveModes);
_transition=0;
_dcp._resize(_nActiveModes);
_dcm._resize(_nActiveModes);
_da._resize(_nActiveModes);
if(_mod) delete [] _mod;
_mod=new CModulation<T>*[_nActiveModes];
for(int i=0;i<_nActiveModes;++i) _mod[i]=NULL;
_zero();
}


void _setPitch(T freq,int transition=256){
if(_playing){
  
  const T ds=2.0*_contactDamping*_auxDampScale/(T)_sampleRate;
  const T fs=_auxFreqScale/(T)_sampleRate;
  
  for(int i=0;i<_freq._size;++i){
     _transition[i]=0;
  }

  T prev=_freq[0];
  T t1=exp(-_damp[0]*ds);
  T t2=M_TWOPI*freq*fs;
  T creal=cos(t2)*t1;
  T cimg=sin(t2)*t1;
  T t3=sqrt(1.0-cimg*cimg);
  _dcp[0]=(creal+t3-_cplus[0])/(T)transition;
  _dcm[0]=(creal-t3-_cminus[0])/(T)transition;
  _da[0]=(_amp[0]*cimg*_auxAmpScale-_aa[0])/(T)transition;
  
  for(int i=1;i<_freq._size;++i){
     T target=_freq[i]/prev*freq;
     t1=exp(-_damp[i]*ds);
     t2=M_TWOPI*target*fs;
     creal=cos(t2)*t1;
     cimg=sin(t2)*t1;
     t3=sqrt(1.0-cimg*cimg);
     _dcp[i]=(creal+t3-_cplus[i])/(T)transition;
     _dcm[i]=(creal-t3-_cminus[i])/(T)transition;
     _da[i]=(_amp[i]*cimg*_auxAmpScale-_aa[i])/(T)transition;
  }
  
  for(int i=0;i<_freq._size;++i){
     _transition[i]=transition;
  }

}else{
  T prev=_freq[0];
  _freq[0]=freq;
  for(int i=1;i<_freq._size;++i){
     _freq[i]=_freq[i]/prev*freq;
  }
  _calcCoefficients();
}
}


void _setModulation(int mode,CModulation<T>* m){
_mod[mode]=m;
}


void _zero(){
_calcCoefficients();
for(int i=0;i<_nActiveModes;i++){
   _u[i]=0; 
   _v[i]=0;
}
}


void _calcCoefficients(){
T t1;
T t2;
T t3;
T creal;
T cimg;
const T ds=2.0*_contactDamping*_auxDampScale/(T)_sampleRate;
const T fs=_auxFreqScale/(T)_sampleRate;

T fTemp;
for(int m=0;m<_nActiveModes;m++){
   fTemp=_freq[m]*fs;
   if(fTemp>0.5) _aa[m]=0.0;
   else{   
     t1=exp(-_damp[m]*ds);
     t2=M_TWOPI*fTemp;  
     creal=cos(t2)*t1;
     cimg=sin(t2)*t1;
     t3=sqrt(1.0-cimg*cimg);
     _cplus[m]=creal+t3;
     _cminus[m]=creal-t3;
     _aa[m]=_amp[m]*cimg*_auxAmpScale; 
   }
}
}


void _strike(T strength){

}

void _input(int size,T* buf){
for(int i=0;i<size;i++) _in[i]=buf[i];
}


void _render(int size,T* buf){

T uPrev;
T u;
T v;
T cp,cm,amp;

_onepole._render(size,_in._getPointer());

for(int mode=0;mode<_nActiveModes;++mode){
   
   uPrev=_u[mode];
   v=_v[mode]; 
   cp=_cplus[mode];
   cm=_cminus[mode];
   amp=_aa[mode];
   
   if(_transition[mode]>0){
   
     int i;
     for(i=0;i<size;i++){
        if(_mod[mode]){
          _mod[mode]->_render(mode,&_in[i],this);
          cp=_cplus[mode];
          cm=_cminus[mode];
          amp=_aa[mode];
        }
        u=cm*uPrev-v+amp*_in[i]+1e-20;
        v=cp*v+uPrev;
        uPrev=u;
        buf[i]+=_gain*v;
        cp+=_dcp[mode];
        cm+=_dcm[mode];
        amp+=_da[mode];
        --_transition[mode];
        if(_transition[mode]<=0) break;
     }
     
     for(;i<size;i++){
        if(_mod[mode]){
          _mod[mode]->_render(mode,&_in[i],this);
          cp=_cplus[mode];
          cm=_cminus[mode];
          amp=_aa[mode];
        }
        u=cm*uPrev-v+amp*_in[i]+1e-20;
        v=cp*v+uPrev;
        uPrev=u;
        buf[i]+=_gain*v;
     }
     
     _aa[mode]=amp;
     _cplus[mode]=cp;
     _cminus[mode]=cm;
     
   }else{
   
     for(int i=0;i<size;i++){
        if(_mod[mode]){
          _mod[mode]->_render(mode,&_in[i],this);
          cp=_cplus[mode];
          cm=_cminus[mode];
          amp=_aa[mode];
        }
        u=cm*uPrev-v+amp*_in[i]+1e-20;
        v=cp*v+uPrev;
        uPrev=u;
        buf[i]+=_gain*v;
     }
     
   }
   
   _u[mode]=u;
   _v[mode]=v;
}

_in=0;

}


};


#endif