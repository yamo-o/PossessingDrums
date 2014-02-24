

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


#include "CAudioFilter.h"

#ifndef CBUMPNOISE_HH
#define CBUMPNOISE_HH

#pragma once 



template <class T> class CRelatedBumpNoise{
public:

T* _in;
T _gain;
T _inputGain;
int _sampleRate;
CLowPassFilter<T> _filter;

CRelatedBumpNoise(){
_in=NULL;
}

~CRelatedBumpNoise(){
if(_in) free(_in);
}

void _init(unsigned int size,unsigned int sampleRate,T cutoff,T gain){
_in=(T*)malloc(size*sizeof(T));
_sampleRate=sampleRate;
_filter._setCutoff(cutoff,sampleRate);
_gain=gain;
_inputGain=1.0;
}

void _input(int size,T* buf){
for(int i=0;i<size;i++){
   _in[i]=_inputGain*buf[i];
}
}

T _render(T input){
T y=_gain*_inputGain*input*random(-1.0,1.0);
return(_filter._render(y));
}


void _render(int size,T* buf){
T y;
for(int i=0;i<size;i++){
   y=random(-1.0,1.0);
   buf[i]+=_gain*_in[i]*y;
}
_filter._render(size,buf);
}

};



#endif