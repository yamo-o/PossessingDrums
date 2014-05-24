

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


#ifndef POSSESSINGDRUMS_HH
#define POSSESSINGDRUMS_HH

#include "CRealtimeAudioSeparator.h"



template <class T> class PossessingDrums{

private:


unsigned int _num_timbres;
unsigned int _windowSize;
unsigned int _sampleRate;


CFFT_Analyzer<T> _fft;
CFFTOverlapedInputBuffer<T> _inputBuf;
CFFTOverlapedOutputBuffer<T>* _outputBuf;

CVector<T> _power;
CVector< CVector2d<T> > _spectrum;
CVector< CVector2d<T> > _timbres_spectrum;
CVector<T> _tmp;


CRealtimeAudioSeparator<T> _separator;
CAudioModeFilter<T>* _inverseFilters;
CModalResonator<T>* _resonators;
CVector<bool> _active;


public:



PossessingDrums():_outputBuf(NULL),_inverseFilters(NULL),_resonators(NULL){}

~PossessingDrums(){
if(_outputBuf) delete [] _outputBuf;
if(_inverseFilters) delete [] _inverseFilters;
if(_resonators) delete [] _resonators;
}


void _init(unsigned int num_timbres,unsigned int windowSize,unsigned int maxBufferSize,unsigned int sampleRate){

_num_timbres=num_timbres;
_windowSize=windowSize;
_sampleRate=sampleRate;

_fft._init(windowSize);

_inputBuf._init(windowSize,windowSize>>1);
_outputBuf=new CFFTOverlapedOutputBuffer<T>[num_timbres];
for(int i=0;i<num_timbres;++i){
   _outputBuf[i]._init(windowSize,maxBufferSize,windowSize>>1);
}

_fft._init(windowSize);
_power._resize(windowSize);
_spectrum._resize(windowSize);
_timbres_spectrum._resize(windowSize);
_tmp._resize(windowSize);


_separator._init(10,_windowSize,_num_timbres,8,0.01);
_inverseFilters=new CAudioModeFilter<T>[num_timbres];
_resonators=new CModalResonator<T>[num_timbres];

_active._resize(num_timbres);
_active=false;

}


void _setPair(const int timbre,const char* drivefile,const char* timbrefile){

CWaveData<T> wav1;
wav1._open(drivefile);
wav1._cache_from_file();
wav1._stereo_to_mono();

// calcurate the time average power spectrum
unsigned int num_sp=calc_average_power_spectrum(_fft,wav1._data,_power._getPointer(),0,wav1._length,_windowSize/2);
_separator._registerTeacher(timbre,_power._getPointer(),num_sp); // register audio source


CAudioModalizer<T> analysis1; // analyzer
analysis1._multiLevelEstimation(wav1); // estimate the modal parameters
_inverseFilters[timbre]._constructInverseFilter(analysis1._num_modes/4,&(analysis1._freq[0]),&(analysis1._ampcoupling[0]),_sampleRate);


_separator._registerTeacher(timbre,_power._getPointer(),num_sp);

CWaveData<T> wav2;
wav2._open(timbrefile);
wav2._cache_from_file();
wav2._stereo_to_mono();

CAudioModalizer<T> analysis2; // analyzer
analysis2._multiLevelEstimation(wav2); // analyze the timbre
_resonators[timbre]._init(_windowSize,_sampleRate,1,1,1,1); // initialize the resonator
_resonators[timbre]._setModes(analysis2._num_modes,&analysis2._ampcoupling[0],&analysis2._freq[0],&analysis2._damping[0]);

_active[timbre]=true;

}



void _ready(){

_separator._calcWeight(); // pre-processing

}



void _record(unsigned int size,T* buf){

int remains=size; // recording buffer size
while(remains>0){
     int writeSize=_inputBuf._record(remains,buf); //buffering
     buf+=writeSize;
     remains-=writeSize;
     if(_inputBuf._filled()){
       _fft._process(&_inputBuf[0],&_spectrum[0][0]); // convert to frequency domain.
       for(int j=0;j<(_windowSize>>1);++j){
          _power[j]=_spectrum[j+1]._normSquare(); // power spectrum
       }
       _separator._process(_power._getPointer()); // sound source separation in spectral domain
       for(int t=0;t<_num_timbres;++t){
          _separator._get_current_result(t,_spectrum._getPointer(),_timbres_spectrum._getPointer()); // get each timbre
          _fft._inverse(&_timbres_spectrum[0][0],_tmp._getPointer()); // convert to time domain signal
          _outputBuf[t]._record(_tmp._getPointer()); // buffering because I use overlapped fft.
       }
     }
}

}



void _render(unsigned int size,T* buf){

for(int t=0;t<_num_timbres;++t){ // for each timbre
   if(_active[t]){
     if(_outputBuf[t]._get_size()>=size){
       _tmp=0;
       _outputBuf[t]._render(size,_tmp._getPointer()); // take separated signals from buffer
       _inverseFilters[t]._render(size,_tmp._getPointer()); // extract the driving signal
       _resonators[t]._input(size,_tmp._getPointer()); // input to resonator
       _resonators[t]._render(size,buf); // modal sound synthesis
     }
   }
}

}


};




#endif
