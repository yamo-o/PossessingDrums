

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

#ifndef CAUDIOMODALIZER_HH
#define CAUDIOMODALIZER_HH

#include "CModalResonator.h"

template <class T> class CModalResidualTransfer{

public:

unsigned int _windowSize;
unsigned int _numLevel;
unsigned int _currentPos;
unsigned int _currentLevel;
unsigned int _prevPos;
unsigned int _prevLevel;
T _accum;

T _gain;
CDenseMatrix<T> _signal;


CModalResidualTransfer(){}

void _init(unsigned int windowSize,unsigned int numLevels){
_windowSize=windowSize;
_numLevel=numLevels;
_signal._resize(_numLevel,_windowSize);
_signal=3;
_reset();
_gain=0.01;
}

void _setGain(T g){
_gain=g;
}

void _setSignal(T refAmp,T* signal){
const T stride=1.0/(T)_numLevel;
unsigned int lv=(refAmp/(T)_windowSize)/stride;
if(fabs(_signal[lv][0])<2.0){
  for(int i=0;i<_windowSize;++i) _signal[lv][i]=signal[i];
}else{
  for(int i=0;i<_windowSize;++i) _signal[lv][i]=0.5*(_signal[lv][i]+signal[i]);
}
}

void _ready(){

unsigned int prev=0;

for(int i=0;i<_windowSize;++i) _signal[0][i]=0;

for(int i=1;i<_signal._row;++i){

   unsigned int j=0;
   while(fabs(_signal[i][0])>=2.0){
        ++j;
        ++i;
        if(i>=_signal._row) break;
   }

   if(j>0){
     if(i<_signal._row){
       for(int k=0;k<j;++k){
          for(int l=0;l<_windowSize;++l){
             _signal[prev+k+1][l]=(1-(k+1)/(T)(j+1))*_signal[prev][l]+(k+1)/(T)(j+1)*_signal[i][l];
          }  
       }
       prev=i;
     }else{
       for(int k=0;k<j;++k){
          for(int l=0;l<_windowSize;++l){
             _signal[prev+k+1][l]=_signal[prev][l];
          }
        }
     }
   }
   
}

}


void _reset(){
_currentLevel=0;
_currentPos=0;
_prevPos=_windowSize/2;
_prevPos=0;
_accum=0;
_ready();
}


inline T _render(T ref){
const T stride=1.0/(T)_numLevel;
T res=_signal[_currentLevel][_currentPos]+_signal[_prevLevel][_prevPos];

_accum+=ref;
++_currentPos;
if(_currentPos>=_windowSize){
  unsigned int idx=_accum/stride;
  _currentLevel=idx;
  _currentPos=0;
  _accum=0;
}

++_prevPos;
if(_prevPos>=_windowSize){
  _prevLevel=_currentLevel;
  _prevPos=0;
}

return(_gain*ref*res);
}


};



template <class T> class CAudioModalizer{

public:

unsigned int _num_modes;
unsigned int _windowSize;
unsigned int _sampleRate;
unsigned int _step;

CVector<int> _ampNum;
CVector<T> _ampSum;
CVector<int> _nPeaks;
CVector<T> _amp;
CVector<T> _ampMoment;
std::vector<T> _ampcoupling;
std::vector<T> _damping;
std::vector<T> _freq;
int _min_peaks;
T _threthold;
int _current_frame;
CModalResidualTransfer<T>* _residual;



CAudioModalizer():_residual(NULL){}
~CAudioModalizer(){
if(_residual) delete [] _residual;
}

void _init(unsigned int windowSize,unsigned int sampleRate,T fractional_window_gap=1,T threthold=0,int min_peaks=0){

_windowSize=windowSize;
_sampleRate=sampleRate;

_step=(int)(_windowSize*fractional_window_gap);
_amp._resize(_windowSize/2);
_ampSum._resize(_windowSize/2);
_ampMoment._resize(_windowSize/2);
_ampNum._resize(_windowSize/2);
_nPeaks._resize(_windowSize/2);


_threthold=0.1;
_min_peaks=4;
_num_modes=0;

}


void _clear(){
_ampcoupling.clear();
_damping.clear();
_freq.clear();
_num_modes=0;
if(_residual) delete [] _residual;
_residual=NULL;
}


void _save(const char* filename){

FILE* fp;
fp=fopen(filename,"wb");

fwrite(&_num_modes,sizeof(unsigned int),1,fp);
for(int i=0;i<_num_modes;++i){
   fwrite(&_ampcoupling[i],sizeof(T),1,fp);
   fwrite(&_damping[i],sizeof(T),1,fp);
   fwrite(&_freq[i],sizeof(T),1,fp);
}

fclose(fp);

}

void _load(const char* filename){

FILE* fp;
fp=fopen(filename,"rb");

T tmp;
fread(&_num_modes,sizeof(unsigned int),1,fp);
for(int i=0;i<_num_modes;++i){
   fread(&tmp,sizeof(T),1,fp);
   _ampcoupling._push(tmp);
   fread(&tmp,sizeof(T),1,fp);
   _damping._push(tmp);
   fread(&tmp,sizeof(T),1,fp);
   _freq._push(tmp);
}

fclose(fp);

}


int _multiLevelEstimation(CWaveData<T>& wav,int maxIter=0,bool calctransfer=false){

const int windowSize[]={128,256,512,1024,2048,4096,8192,16384};

_sampleRate=wav._sampleRate;

int best=0;
int num_modes=0;

for(int i=0;i<8;++i){
   CFFT_Analyzer<T> fft;
   fft._init(windowSize[i]);
   _init(windowSize[i],_sampleRate);
   CDenseMatrix<T> sp;
   fft._spectrogram(wav._data,wav._length,windowSize[i]/2,sp);
   _current_frame=0;
   _amp=0;
   _ampSum=0;
   _ampNum=0;
   _nPeaks=0;
   _ampMoment=0;
   _ampcoupling.clear();
   _freq.clear();
   _damping.clear();
   for(int j=0;j<sp._row;++j){
      _addFrame(windowSize[i],sp[j]);
   }
   _process();
   if(_num_modes>num_modes){
     best=i;
     num_modes=_num_modes;
   }
}

CFFT_Analyzer<T> fft;
fft._init(windowSize[best]);
_init(windowSize[best],_sampleRate);
CDenseMatrix<T> sp;
fft._spectrogram(wav._data,wav._length,windowSize[best]/2,sp);
_current_frame=0;
_amp=0;
_ampSum=0;
_ampNum=0;
_nPeaks=0;
_ampMoment=0;
_ampcoupling.clear();
_freq.clear();
_damping.clear();
for(int j=0;j<sp._row;++j){
   _addFrame(windowSize[best],sp[j]);
}
_process();

_residual=new CModalResidualTransfer<T>[_num_modes];

T step=_sampleRate/(T)_windowSize*0.1;
T stride=_sampleRate/(T)_windowSize;
T* estimated=new T[wav._length];

//greedy approach
for(int i=0;i<_num_modes;++i){
   
   T plus=-1;
   T minus=-1;
   int iter=0;
   T limit[2];
   if(i>0) limit[0]=0.5*(_freq[i]+_freq[i-1]);
   else limit[0]=0;
   if(i<_num_modes-1) limit[1]=0.5*(_freq[i]+_freq[i+1]);
   else limit[1]=0.5*(_freq[i]+_sampleRate*0.5);
   while(iter<maxIter){
        _createAudio(estimated,std::max((int)wav._length/4,4),i);
        CDenseMatrix<T> esp;
        fft._spectrogram(estimated,std::max((int)wav._length/4,4),windowSize[best]/2,esp);
        T start=(i>0)?(0.5*(_freq[i-1]+_freq[i])):0.5*_freq[i];
        T end=(i<_num_modes-1)?(0.5*(_freq[i]+_freq[i+1])):0.5*(_freq[i]+_sampleRate*0.5);
        T current=_matchCost(sp,esp,0,esp._row,start/stride,end/stride);
       
        if(plus<0){
          _freq[i]+=step;
          _createAudio(estimated,std::max((int)wav._length/4,4),i);
          fft._spectrogram(estimated,std::max((int)wav._length/4,4),windowSize[best]/2,esp);
          plus=_matchCost(sp,esp,0,esp._row,start/stride,end/stride);
          _freq[i]-=step;
        }
        if(minus<0){
          _freq[i]-=step;
          _createAudio(estimated,std::max((int)wav._length/4,4),i);
          fft._spectrogram(estimated,std::max((int)wav._length/4,4),windowSize[best]/2,esp);
          minus=_matchCost(sp,esp,0,esp._row,start/stride,end/stride);
          _freq[i]+=step;
        }

        if(current<plus && current<minus) break;
        if(current<sp._col*1.0e-3) break;
        if(_freq[i]<=limit[0]) break;
        else if(_freq[i]>=limit[1]) break;
       
        if(plus<minus){
          _freq[i]+=step;
          plus=-1;
          minus=current;
        }else{
          _freq[i]-=step;
          plus=current;
          minus=-1;
        }
       
        if(i>0){
          
        }

        ++iter;
   }

   
   plus=minus=-1;
   iter=0;
   while(iter<maxIter){
        T ampstep=_ampcoupling[i]*0.01;
        _createAudio(estimated,std::max((int)wav._length/10,4),i);
        CDenseMatrix<T> esp;
        fft._spectrogram(estimated,std::max((int)wav._length/10,4),windowSize[best]/2,esp);
        T start=(i>0)?(0.5*(_freq[i-1]+_freq[i])):_freq[i];
        T end=(i<_num_modes-1)?(0.5*(_freq[i]+_freq[i+1])):_sampleRate*0.5;
        T current=_matchCost(sp,esp,0,esp._row,start/stride,end/stride);
       
        if(plus<0){
          _ampcoupling[i]+=ampstep;
          _createAudio(estimated,std::max((int)wav._length/10,4),i);
          fft._spectrogram(estimated,std::max((int)wav._length/10,4),windowSize[best]/2,esp);
          plus=_matchCost(sp,esp,0,esp._row,start/stride,end/stride);
          _ampcoupling[i]-=ampstep;
        }
       
        if(minus<0){
          _ampcoupling[i]-=ampstep;
          _createAudio(estimated,std::max((int)wav._length/10,4),i);
          fft._spectrogram(estimated,std::max((int)wav._length/10,4),windowSize[best]/2,esp);
          minus=_matchCost(sp,esp,0,esp._row,start/stride,end/stride);
          _ampcoupling[i]+=ampstep;
        }

        if(current<plus && current<minus) break;
        if(current<sp._col*1.0e-3) break;
        
        if(plus<current){
          _ampcoupling[i]+=ampstep;
          plus=-1;
          minus=current;
        }else{
          _ampcoupling[i]-=ampstep;
          plus=current;
          minus=-1;
        }

        ++iter;
   }

   
   plus=minus=-1;
   iter=0;
   while(iter<maxIter){
        T step=_damping[i]*0.01;
        _createAudio(estimated,wav._length,i);
        CDenseMatrix<T> esp;
        fft._spectrogram(estimated,wav._length,windowSize[best]/2,esp);
        T start=(i>0)?(0.5*(_freq[i-1]+_freq[i])):_freq[i];
        T end=(i<_num_modes-1)?(0.5*(_freq[i]+_freq[i+1])):_sampleRate*0.5;
        T current=_matchCost(sp,esp,0,esp._row,start/stride,end/stride);
       
        if(plus<0){
          _damping[i]+=0.01;
          _createAudio(estimated,wav._length,i);
          fft._spectrogram(estimated,wav._length,windowSize[best]/2,esp);
          plus=_matchCost(sp,esp,0,esp._row,start/stride,end/stride);
          _damping[i]-=0.01;
        }
       
        if(minus<0){
          _damping[i]-=0.01;
          _createAudio(estimated,wav._length,i);
          fft._spectrogram(estimated,wav._length,windowSize[best]/2,esp);
          minus=_matchCost(sp,esp,0,esp._row,start/stride,end/stride);
          _damping[i]+=0.01;
        }

        if(current<plus && current<minus) break;
        if(current<sp._col*1.0e-3) break;
        
        if(plus<current){
          _damping[i]+=step;
          plus=-1;
          minus=current;
        }else{
          _damping[i]-=step;
          plus=current;
          minus=-1;
        }
        if(_damping[i]<1.0){
          _damping[i]=0.1; 
          break;
        }

        ++iter;
   }

}


if(calctransfer){
  CDenseMatrix< CVector2d<T> > original;
  fft._spectrogram(wav._data,wav._length,windowSize[best]/2,original);
  CDenseMatrix< CVector2d<T> > final;
  _createAudio(estimated,wav._length);
  fft._spectrogram(estimated,wav._length,windowSize[best]/2,final);

  CVector< CVector2d<T> > residual(windowSize[best]);
  CVector<T> sig(windowSize[best]);

  const CVector2d<T> zero(0);

  for(int m=0;m<_num_modes;++m){
     T start=(m>0)?(0.5*(_freq[m-1]+_freq[m])):1;
     T end=(m<_num_modes-1)?(0.5*(_freq[m]+_freq[m+1])):_sampleRate*0.5;   
     residual=zero;
     T total=0;
     T power_res=0;
     for(int i=0;i<original._col;++i){
        T power_mode=0;
        for(int j=start/stride;j<end/stride;++j){
           T a=original[j][i]._norm()-final[j][i]._norm();
           if(a>0) residual[j]=a*original[j][i]/original[j][i]._norm();
           else residual[j]=0;
           residual[windowSize[best]-j]=residual[j];
           power_mode+=final[j][i]._normSquare();
           power_res+=a*a;
        }
        total+=power_mode;
        fft._inverse(&residual[0][0],sig._getPointer());
        _residual[m]._setSignal(power_mode,sig._getPointer());
     }
     _residual[m]._setGain(sqrt(power_res/total));
  }
}

        
delete [] estimated;


return(best);
}



T _matchCost(CDenseMatrix<T>& m1,CDenseMatrix<T>& m2,int rowstart,int rowend,int colstart,int colend){

T residual=0;
for(int j=rowstart;j<rowend;++j){
   for(int i=colstart;i<colend;++i){
      residual+=fabs(m2[j][i]-m1[j][i]);
   }
}

return(residual);
}

void _createAudio(T* data,int frames){

T time_resolution=1.0/(T)_sampleRate;

memset(data,0,frames*sizeof(T));

for(int i=0;i<_num_modes;++i){
   T t=0;
   for(int j=0;j<frames;++j){
      data[j]+=_ampcoupling[i]*exp(-_damping[i]*t)*sin(2*M_PI*_freq[i]*t);
      t+=time_resolution;
   }
}

}

void _createAudio(T* data,int frames,int mode){

T time_resolution=1.0/(T)_sampleRate;

T t=0;
for(int j=0;j<frames;++j){
   data[j]=_ampcoupling[mode]*exp(-_damping[mode]*t)*sin(2*M_PI*_freq[mode]*t);
   t+=time_resolution;
}

}

void _addFrame(int windowSize,T* source){

for(int n=0;n<(windowSize/2);n++){
   _amp[n]=0.5*log(source[n]);
 
   if(_amp[n]>_threthold){
     _ampSum[n]+=_amp[n];
     _ampMoment[n]+=_current_frame*_amp[n];
     _ampNum[n]++;
   }
}

for(int n=1;n<((windowSize/2)-1);n++){
   if(_amp[n]>_threthold && _amp[n-1]<_amp[n] && _amp[n]>_amp[n+1]){
     _nPeaks[n]++;
   }
}

_current_frame++;

}


void _process(){

T gradient;
T offset;

int nModesFound=0;
T overlapFactor=(T)_windowSize/(T)_step;
T maxAmpcoupling=0.0;
int maxPeaks;

//McAaulay-Quatieri Analysis
//First, I roughly estimate the modes.
int maxidx=0;
do{
  maxPeaks=0;
  int maxPeakIndex=0;
  
  for(int i=1;i<((_windowSize/2)-1);i++){
      if(_nPeaks[i]>maxPeaks){ maxPeaks=_nPeaks[i]; maxPeakIndex=i; }
  }
  
  if(maxPeaks<_min_peaks) break;
  
  _nPeaks[maxPeakIndex]=0;
  int idx=maxPeakIndex;

  T n=_ampNum[idx]+2;
        
  T a=(n-1.0)*0.5;
  T b=12.0/(n*(n-1.0)*(n-1.0));
  T c=n*(n-1.0)*0.5;
  T d=1.0/n;

  gradient=2.0*(_ampMoment[idx]-a*_ampSum[idx])*b;
  offset=(_ampSum[idx]-c*gradient)*d;
  
  if(gradient<0.0){
    T freq=idx/(T)_windowSize*_sampleRate;
    _freq.push_back(freq);
    T damp=_sampleRate*(-gradient)/(T)_step;
    _damping.push_back(damp);
    T p=overlapFactor*(-gradient);
    T amp=exp(offset)*p/(1.0-exp(-p));
    _ampcoupling.push_back(amp);

    if(amp>maxAmpcoupling){
      maxAmpcoupling=amp;
      maxidx=idx;
    }
   
    nModesFound++;
  }
  
      
}while(1);


for(int i=0;i<nModesFound;i++){
  _ampcoupling[i]/=maxAmpcoupling;
}

T record;
int recordIndex;
T t;

for(int i=0;i<nModesFound-1;i++){
   record=_freq[i];
   recordIndex=i;
   for(int j=i+1;j<nModesFound;j++){
     if(_freq[j]<record){
      record=_freq[j];
      recordIndex=j;
     }
   }
   if(recordIndex>i){
    t=_ampcoupling[i];
    _ampcoupling[i]=_ampcoupling[recordIndex];
    _ampcoupling[recordIndex]=t;

    t=_damping[i];
    _damping[i]=_damping[recordIndex];
    _damping[recordIndex]=t;

    t=_freq[i];
    _freq[i]=_freq[recordIndex];
    _freq[recordIndex]=t;
   }
}

_num_modes=nModesFound;

}


};


#endif