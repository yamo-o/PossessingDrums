
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


#ifndef FFT_HH
#define FFT_HH

#include "CWaveData.h"
#include "CRingBufferMatrix.h"


template <class T> void create_sinbell_window(int windowSize,T* w){
const T f=M_PI/(T)(windowSize-1.0);
for(int n=0;n<windowSize;++n) w[n]=sin(f*n);
}

template <class T> void create_hanning_window(int windowSize,T* w){
const T f=2.0*M_PI/(T)(windowSize-1.0);
for(int n=0;n<windowSize;n++) w[n]=0.5*(1.0-cos(f*n));
}

template <class T> void create_hamming_window(int windowSize,T* w){
const T f=2.0*M_PI/(T)(windowSize-1.0);
for(int n=0;n<windowSize;n++) w[n]=0.54-0.46*cos(f*n);
}

template <class T> void create_blackman_window(int windowSize,T* w){
const T f=2.0*M_PI/(T)(windowSize-1.0);
for(int n=0;n<windowSize;n++) w[n]=0.42-0.5*cos(f*n)+0.08*cos(2*f*n);
}


typedef enum FFT_WINDOW_TYPE{
FFT_WINDOW_TYPE_HANNING,
FFT_WINDOW_TYPE_HAMMING,
FFT_WINDOW_TYPE_SINBELL,
FFT_WINDOW_TYPE_BLACKMAN,
} FFT_WINDOW_TYPE;



template <class T> class CFFTOverlapedInputBuffer{
public:

CVector<T> _input;
unsigned int _nFFT;
unsigned int _currentPos;
unsigned int _step;


void _init(unsigned int nFFT,unsigned int overlap){
_nFFT=nFFT;
_step=_nFFT-overlap;
_input._resize(_nFFT);
_input=0;
_currentPos=_nFFT-_step;
}

void _reset(){
_input=0;
_currentPos=_nFFT-_step;
}


inline T& operator[](int i){ return(_input[i]); }
inline const T& operator[](int i) const { return(_input[i]); }

inline bool _filled(){ return(_currentPos==(_nFFT-_step)); }
inline unsigned int _getLatency(){ return(_nFFT-_step); }
inline unsigned int _getSpectrumBufferNum(const unsigned int bufferSize){ int mod=(bufferSize%_step)?1:0; return(bufferSize/_step+mod); }


//return written size
int _record(unsigned int size,const T* buf){

const unsigned int step=_step;

if(_currentPos==(_nFFT-step)){
  for(int i=0;i<step;++i){
     _input[i]=_input[i+step];
  }
}

T* in=&_input[_currentPos];
if(_currentPos+size>=_nFFT){
  const int remain=_nFFT-_currentPos;
  for(int i=0;i<remain;++i){
     in[i]=buf[i];
  }
  _currentPos=_nFFT-_step;
  return(remain);
}

for(int i=0;i<size;++i){
   in[i]=buf[i];
}
_currentPos+=size;

return(size);
}


};



template <class T> class CFFTOverlapedOutputBuffer{
public:

CVector<T> _output;
unsigned int _nFFT;
unsigned int _currentPos;
unsigned int _step;


inline T& operator[](int i){ return(_output[i]); }
inline const T& operator[](int i) const { return(_output[i]); }


void _init(unsigned int nFFT,unsigned int maxBufferSize,unsigned int overlap){
_nFFT=nFFT;
_step=_nFFT-overlap;
if(maxBufferSize>_nFFT) _output._resize(_nFFT*(maxBufferSize/_nFFT+2)-overlap);
else _output._resize(2*_nFFT-overlap);
_output=0;
_currentPos=0;
}

inline unsigned int _get_size(){ return(_currentPos); }


// the size of input must be equal to _nFFT
void _record(T* input){
const unsigned int half=(_nFFT>>1);
T* out=&_output[_currentPos];

for(int i=0;i<half;++i){
   out[i]+=input[i];
}
for(int i=half;i<_nFFT;++i){
   out[i]=input[i];
}
_currentPos+=_step;
}


//return written size
int _render(unsigned int size,T* buf){
const unsigned int half=(_nFFT>>1);
const T* out=&_output[_currentPos];

if(_currentPos>=size){
  for(int i=0;i<size;++i){
     buf[i]+=_output[i];
  }
  _currentPos-=size;
  for(int i=0;i<_currentPos+half;++i){
     _output[i]=_output[i+size];
  }
  return(size);
}

for(int i=0;i<_currentPos;++i){
   buf[i]+=_output[i];
}
for(int i=0;i<half;++i){
   _output[i]=out[i];
}
size-=_currentPos;
unsigned int write=_currentPos;
_currentPos=0;

return(write);
}


};



template <class T> void cdft(int n,int isgn,T *a,int *ip,T *w);

template <class T> class CFFT_Analyzer{

public:

int _nFFT;
int* _ip;
T* _sinfunc;
T* _window;
T* _buffer;
unsigned int _currentPos;
CVector< CVector<T> > _realtimeBuffer;


CFFT_Analyzer(){
_ip=NULL;
_sinfunc=NULL;
_window=NULL;
_buffer=NULL;
}

~CFFT_Analyzer(){
_cleanup();
}

void _cleanup(){
if(_ip) delete [] _ip;
_ip=NULL;
if(_sinfunc) delete [] _sinfunc;
_sinfunc=NULL;
if(_window) delete [] _window;
_window=NULL;
if(_buffer) delete [] _buffer;
_buffer=NULL;  
}

void _init(int size){
_cleanup();
_ip=new int[2+(int)sqrt(2*size)+1];
_ip[0]=0;
_sinfunc=new T[size];
_window=new T[size];
_nFFT=size;
for(int i=0;i<size;++i) _window[i]=1;
}

inline int _getWindowSize(){ return(_nFFT); }

void _process(T* data,T* result){
for(int i=0;i<_nFFT;++i){
   result[2*i]=data[i]*_window[i];
   result[2*i+1]=0;
}
cdft<T>(2*_nFFT,-1,result,_ip,_sinfunc);
}

void _process(T* data){
if(!_buffer) _buffer=new T[2*_nFFT];
for(int i=0;i<_nFFT;++i){
   _buffer[2*i]=data[i]*_window[i];
   _buffer[2*i+1]=0;
}
cdft<T>(2*_nFFT,-1,_buffer,_ip,_sinfunc);
}

inline T* _getResult(){ return(_buffer); }

inline void _bufferingComplex(T* data){
if(!_buffer) _buffer=new T[2*_nFFT];
memcpy(_buffer,data,2*_nFFT*sizeof(T));
}

void _inverse(T* data,T* res){
if(!_buffer) _buffer=new T[2*_nFFT];
for(int i=0;i<_nFFT;++i){
   _buffer[2*i]=data[2*i];
   _buffer[2*i+1]=data[2*i+1];
}
cdft<T>(2*_nFFT,1,_buffer,_ip,_sinfunc);
const T invN=1.0/(T)_nFFT;
for(int i=0;i<_nFFT;++i) res[i]=_buffer[2*i]*_window[i]*invN;
}

void _inverse(){
cdft(2*_nFFT,1,_buffer,_ip,_sinfunc);
T invN=1.0/(T)_nFFT;
for(int i=0;i<_nFFT;++i) _buffer[2*i]*=_window[i]*invN;
}

void _setWindow(FFT_WINDOW_TYPE type){
switch(type){
      case FFT_WINDOW_TYPE_HANNING: create_hanning_window(_nFFT,_window); break;
      case FFT_WINDOW_TYPE_HAMMING: create_hamming_window(_nFFT,_window); break;
      case FFT_WINDOW_TYPE_SINBELL: create_sinbell_window(_nFFT,_window); break;
      case FFT_WINDOW_TYPE_BLACKMAN: create_blackman_window(_nFFT,_window); break;
}
}


void _overlap_add(int overlap,int n_frames,T* complex,T* real){
const T invN=1.0/(T)_nFFT;

memset(real,0,(_nFFT+(n_frames-1)*(_nFFT-overlap))*sizeof(T));

for(int i=0;i<_nFFT;++i) real[i]=complex[2*i]*_window[i]*invN;

for(int i=1;i<n_frames;++i){
   for(int w=0;w<_nFFT;++w){
      real[(_nFFT-overlap)*i+w]+=complex[2*(_nFFT*i+w)]*_window[w]*invN;
   }
}

}

int _create_overlapped_frame(int overlap,int length,T* real,T* complex){
int n_frames=ceil((length+overlap)/(T)(_nFFT-overlap));

T val;
for(int j=0;j<n_frames;j++){
   for(int i=0;i<_nFFT;i++){
      if((((_nFFT-overlap)*j+i)<overlap) || ((_nFFT-overlap)*j+i)>=(length+overlap)) val=0.0;
      else val=real[(_nFFT-overlap)*j+i-overlap];
      complex[2*(_nFFT*j+i)]=_window[i]*val;
      complex[2*(_nFFT*j+i)+1]=0.0;
   }
}

return(n_frames);
}

int _calc_num_frames(int overlap,int length){
return(ceil((length+overlap)/(T)(_nFFT-overlap)));
}

int _overlapped_fft(int overlap,int length,T* real,T* complex){
int n_frames=_create_overlapped_frame(overlap,length,real,complex);
for(int i=0;i<n_frames;i++){
   cdft<T>(2*_nFFT,-1,&complex[2*_nFFT*i],_ip,_sinfunc);
}
return(n_frames);
}

void _overlapped_ifft(int overlap,int num_frames,T* complex,T* real){
const T invN=1.0/(T)_nFFT;
memset(real,0,(_nFFT+(num_frames-1)*(_nFFT-overlap))*sizeof(T));
for(int i=0;i<num_frames;i++){
   _bufferingComplex(&complex[2*_nFFT*i]);
   cdft(2*_nFFT,1,_buffer,_ip,_sinfunc);
   for(int w=0;w<_nFFT;++w) real[(_nFFT-overlap)*i]+=_buffer[2*w]*_window[w]*invN;
}
}



void _spectrogram(T* data,int length,int overlap,CDenseMatrix<T>& spectrogram){
int num_frames=_calc_num_frames(overlap,length);
spectrogram._resize(num_frames,_nFFT/2+1);
if(!_buffer){
  _buffer=new T[2*_nFFT];
}

for(int j=0;j<num_frames;++j){
   for(int i=0;i<_nFFT;++i){
      T val=((((_nFFT-overlap)*j+i)<overlap) || ((_nFFT-overlap)*j+i)>=(length+overlap))?0.0:data[(_nFFT-overlap)*j+i-overlap];
      _buffer[2*i]=_window[i]*val;
      _buffer[2*i+1]=0.0;
   }
   cdft<T>(2*_nFFT,-1,_buffer,_ip,_sinfunc);
   for(int i=0;i<=(_nFFT>>1);++i){
      spectrogram[j][i]=_buffer[2*i]*_buffer[2*i]+_buffer[2*i+1]*_buffer[2*i+1];
   }
}
}


void _spectrogram(T* data,int length,int overlap,CDenseMatrix< CVector2d<T> >& spectrogram){
int num_frames=_calc_num_frames(overlap,length);
spectrogram._resize(num_frames,_nFFT/2+1);
if(!_buffer){
  _buffer=new T[2*_nFFT];
}

for(int j=0;j<num_frames;++j){
   for(int i=0;i<_nFFT;++i){
      T val=((((_nFFT-overlap)*j+i)<overlap) || ((_nFFT-overlap)*j+i)>=(length+overlap))?0.0:data[(_nFFT-overlap)*j+i-overlap];
      _buffer[2*i]=_window[i]*val;
      _buffer[2*i+1]=0.0;
   }
   cdft<T>(2*_nFFT,-1,_buffer,_ip,_sinfunc);
   for(int i=0;i<=(_nFFT>>1);++i){
      spectrogram[j][i][0]=_buffer[2*i];
      spectrogram[j][i][1]=_buffer[2*i+1];
   }
}
}


void _inverse(CDenseMatrix< CVector2d<T> >& spectrogram,int overlap,T* data){
_overlapped_ifft(overlap,spectrogram._col,&spectrogram[0][0],data);
}


};


template <class T> unsigned int calc_average_power_spectrum(CFFT_Analyzer<T>& fft,T* data,T* res,int start_pos,int end_pos,int step,T threthold=0.1){
const int windowSize=fft._getWindowSize();
const int HalfWindowSize=windowSize/2;
T* spectrum=new T[2*windowSize];
T* power_spectrum=new T[HalfWindowSize];
int length;

memset(res,0,HalfWindowSize*sizeof(T));

int num_sp=0;
length=end_pos-start_pos;
const int num_frames=length/step-windowSize/step-1;

for(int i=0;i<num_frames;++i){
   fft._process(&data[start_pos+step*i],spectrum);

   T sum=0;
   T power=0.0;
   for(int j=0;j<HalfWindowSize;++j){
      power_spectrum[j]=spectrum[2*(j+1)]*spectrum[2*(j+1)]+spectrum[2*(j+1)+1]*spectrum[2*(j+1)+1];
      if(power_spectrum[j]>power) power=power_spectrum[j];
      sum+=power_spectrum[j];
   }
   
   if(sum>threthold){
     for(int j=0;j<HalfWindowSize;++j){
        res[j]+=sqrt(power_spectrum[j]/power);
     }
     ++num_sp;
   }
   
}


T maxf=0;
if(num_sp>0){
  for(int i=0;i<HalfWindowSize;++i){
     res[i]/=(T)num_sp;
     if(res[i]>maxf) maxf=res[i]; 
  }
  for(int i=0;i<HalfWindowSize;++i){
     res[i]/=maxf;
  }
}

delete [] spectrum;
delete [] power_spectrum;

return(num_sp);
}




template <class T> void sonogram(CDenseMatrix<T>& spectrogram,int sampleRate,int outerear,bool spread,CDenseMatrix<T>& sone){

int windowSize=spectrogram._col*2;
int frames=spectrogram._row;

const int bark_upper[]={100,200,300,400,510,630,770,920,1080,1270,1480,1720,2000,2320,2700,3150,3700,4400,5300,6400,7700,9500,12000,15500};
int num_barks=0;
while(bark_upper[num_barks]<sampleRate*0.5){
     ++num_barks;
     if(num_barks>=24) break;
}

T step=sampleRate/(T)windowSize;
CVector<T> W_Adb(windowSize/2);
T f=0;
switch(outerear){
      case 1: // terhardt 1979 (calculating virtual pitch, hearing research #1, pp 155-182)
        for(int i=0;i<windowSize/2;++i){
           W_Adb[i]=square(pow(10,((-3.64*pow(f*0.001,-0.8)
                    +6.5*exp(square(-0.6*(f*0.001-3.3)))
                    -0.001*pow(f*0.001,4)/20.0))));
           f+=step;
        }
        break;
      case 2: // less emph around 4Hz, more emphasis on low freqs
        for(int i=0;i<windowSize/2;++i){        
           W_Adb[i]=square(pow(10,((0.6*-3.64*pow(f*0.001,-0.8)
                    +0.5*exp(-0.6*square(f*0.001-3.3))
                    -0.001*pow(f*0.001,4)/20.0))));
           f+=step;
        }
        break;
      default: // all weighted equally
        W_Adb=1;
        break;
}


CDenseMatrix<T> dlinear(windowSize/2,frames); // data from fft (linear freq scale)
sone._resize(num_barks,frames); // data after bark scale
for(int j=0;j<frames;++j){
   for(int i=0;i<windowSize/2;++i){
      dlinear[j][i]=W_Adb[i]*spectrogram[j][i];
   }
}

sone=0;
for(int j=0;j<frames;++j){
   f=0;
   int w=0;
   for(int i=0;i<num_barks;++i){
      for(;f<bark_upper[i] && w<windowSize/2;++w,f+=step){
         sone[j][i]+=dlinear[j][w];
      }
   }
}

// spreading function: schroeder et al., 1979, JASA, optimizing digital speech coders by exploiting masking properties of the human ear
if(spread){
  for(int j=0;j<sone._row;++j){
     for(int i=0;i<sone._col;++i){
        T s=pow(10,((15.81+7.5*((j-i)+0.474)-17.5*sqrt(1+square((j-i)+0.474)))*0.1));
        sone[j][i]*=s;
     }
  }
}

//convert to dB
for(int j=0;j<sone._row;++j){
   for(int i=0;i<sone._col;++i){
      if(sone[j][i]<1) sone[j][i]=1;
      sone[j][i]=10*log10(sone[j][i]);
      // bladon and lindblom, 1981, JASA, modelling the judment of vowel quality differences
      if(sone[j][i]>=40){
        sone[j][i]=pow(2,0.1*(sone[j][i]-40));
      }else{
        sone[j][i]=pow(sone[j][i]/40.0,2.642);
      }
   }
}

}



//This project uses modified version of Ooura's FFT code(http://www.kurims.kyoto-u.ac.jp/~ooura/index.html). 


inline void makeipt(int nw,int *ip){
int j,l,m,m2,p,q;
    
ip[2]=0;
ip[3]=16;
m=2;
for(l=nw;l>32;l>>=2){
   m2=m<<1;
   q=m2<<3;
   for(j=m;j<m2;j++){
       p=ip[j]<<2;
       ip[m+j]=p;
       ip[m2+j]=p+q;
   }
   m=m2;
}
}


template <class T> void makewt(int nw,int *ip,T *w){

int j,nwh,nw0,nw1;
T delta,wn4r,wk1r,wk1i,wk3r,wk3i;
    
ip[0]=nw;
ip[1]=1;
if(nw>2){
  nwh=nw>>1;
  delta=atan(1.0)/nwh;
  wn4r=cos(delta*nwh);
  w[0]=1;
  w[1]=wn4r;
  if(nwh==4){
    w[2]=cos(delta*2);
    w[3]=sin(delta*2);
  }else if(nwh>4){
    makeipt(nw,ip);
    w[2]=0.5/cos(delta*2);
    w[3]=0.5/cos(delta*6);
    for(j=4;j<nwh;j+=4){
       w[j]=cos(delta*j);
       w[j+1]=sin(delta*j);
       w[j+2]=cos(3*delta*j);
       w[j+3]=-sin(3*delta*j);
    }
  }
  nw0=0;
  while(nwh>2){
       nw1=nw0+nwh;
       nwh>>=1;
       w[nw1]=1;
       w[nw1+1]=wn4r;
       if(nwh==4){
         wk1r=w[nw0+4];
         wk1i=w[nw0+5];
         w[nw1+2]=wk1r;
         w[nw1+3]=wk1i;
       }else if(nwh>4){
         wk1r=w[nw0+4];
         wk3r=w[nw0+6];
         w[nw1+2]=0.5/wk1r;
         w[nw1+3]=0.5/wk3r;
         for(j=4;j<nwh;j+=4){
            wk1r=w[nw0+2*j];
            wk1i=w[nw0+2*j+1];
            wk3r=w[nw0+2*j+2];
            wk3i=w[nw0+2*j+3];
            w[nw1+j]=wk1r;
            w[nw1+j+1]=wk1i;
            w[nw1+j+2]=wk3r;
            w[nw1+j+3]=wk3i;
         }
       }
       nw0=nw1;
  }
}

}

template <class T> void cftmdl1(int n,T *a,T *w){
int j,j0,j1,j2,j3,k,m,mh;
T wn4r,wk1r,wk1i,wk3r,wk3i;
T x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i;
    
mh=n>>3;
m=2*mh;
j1=m;
j2=j1+m;
j3=j2+m;
x0r=a[0]+a[j2];
x0i=a[1]+a[j2+1];
x1r=a[0]-a[j2];
x1i=a[1]-a[j2+1];
x2r=a[j1]+a[j3];
x2i=a[j1+1]+a[j3+1];
x3r=a[j1]-a[j3];
x3i=a[j1+1]-a[j3+1];
a[0]=x0r+x2r;
a[1]=x0i+x2i;
a[j1]=x0r-x2r;
a[j1+1]=x0i-x2i;
a[j2]=x1r-x3i;
a[j2+1]=x1i+x3r;
a[j3]=x1r+x3i;
a[j3+1]=x1i-x3r;
wn4r=w[1];
k=0;
for(j=2;j<mh;j+=2){
   k+=4;
   wk1r=w[k];
   wk1i=w[k+1];
   wk3r=w[k+2];
   wk3i=w[k+3];
   j1=j+m;
   j2=j1+m;
   j3=j2+m;
   x0r=a[j]+a[j2];
   x0i=a[j+1]+a[j2+1];
   x1r=a[j]-a[j2];
   x1i=a[j+1]-a[j2+1];
   x2r=a[j1]+a[j3];
   x2i=a[j1+1]+a[j3+1];
   x3r=a[j1]-a[j3];
   x3i=a[j1+1]-a[j3+1];
   a[j]=x0r+x2r;
   a[j+1]=x0i+x2i;
   a[j1]=x0r-x2r;
   a[j1+1]=x0i-x2i;
   x0r=x1r-x3i;
   x0i=x1i+x3r;
   a[j2]=wk1r*x0r-wk1i*x0i;
   a[j2+1]=wk1r*x0i+wk1i*x0r;
   x0r=x1r+x3i;
   x0i=x1i-x3r;
   a[j3]=wk3r*x0r+wk3i*x0i;
   a[j3+1]=wk3r*x0i-wk3i*x0r;
   j0=m-j;
   j1=j0+m;
   j2=j1+m;
   j3=j2+m;
   x0r=a[j0]+a[j2];
   x0i=a[j0+1]+a[j2+1];
   x1r=a[j0]-a[j2];
   x1i=a[j0+1]-a[j2+1];
   x2r=a[j1]+a[j3];
   x2i=a[j1+1]+a[j3+1];
   x3r=a[j1]-a[j3];
   x3i=a[j1+1]-a[j3+1];
   a[j0]=x0r+x2r;
   a[j0+1]=x0i+x2i;
   a[j1]=x0r-x2r;
   a[j1+1]=x0i-x2i;
   x0r=x1r-x3i;
   x0i=x1i+x3r;
   a[j2]=wk1i*x0r-wk1r*x0i;
   a[j2+1]=wk1i*x0i+wk1r*x0r;
   x0r=x1r+x3i;
   x0i=x1i-x3r;
   a[j3]=wk3i*x0r+wk3r*x0i;
   a[j3+1]=wk3i*x0i-wk3r*x0r;
}
j0=mh;
j1=j0+m;
j2=j1+m;
j3=j2+m;
x0r=a[j0]+a[j2];
x0i=a[j0+1]+a[j2+1];
x1r=a[j0]-a[j2];
x1i=a[j0+1]-a[j2+1];
x2r=a[j1]+a[j3];
x2i=a[j1+1]+a[j3+1];
x3r=a[j1]-a[j3];
x3i=a[j1+1]-a[j3+1];
a[j0]=x0r+x2r;
a[j0+1]=x0i+x2i;
a[j1]=x0r-x2r;
a[j1+1]=x0i-x2i;
x0r=x1r-x3i;
x0i=x1i+x3r;
a[j2]=wn4r*(x0r-x0i);
a[j2+1]=wn4r*(x0i+x0r);
x0r=x1r+x3i;
x0i=x1i-x3r;
a[j3]=-wn4r*(x0r+x0i);
a[j3+1]=-wn4r*(x0i-x0r);
}


template <class T> void cftmdl2(int n,T *a,T *w){
int j,j0,j1,j2,j3,k,kr,m,mh;
T wn4r,wk1r,wk1i,wk3r,wk3i,wd1r,wd1i,wd3r,wd3i;
T x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,y0r,y0i,y2r,y2i;
    
mh=n>>3;
m=2*mh;
wn4r=w[1];
j1=m;
j2=j1+m;
j3=j2+m;
x0r=a[0]-a[j2+1];
x0i=a[1]+a[j2];
x1r=a[0]+a[j2+1];
x1i=a[1]-a[j2];
x2r=a[j1]-a[j3+1];
x2i=a[j1+1]+a[j3];
x3r=a[j1]+a[j3+1];
x3i=a[j1+1]-a[j3];
y0r=wn4r*(x2r-x2i);
y0i=wn4r*(x2i+x2r);
a[0]=x0r+y0r;
a[1]=x0i+y0i;
a[j1]=x0r-y0r;
a[j1+1]=x0i-y0i;
y0r=wn4r*(x3r-x3i);
y0i=wn4r*(x3i+x3r);
a[j2]=x1r-y0i;
a[j2+1]=x1i+y0r;
a[j3]=x1r+y0i;
a[j3+1]=x1i-y0r;
k=0;
kr=2*m;
for(j=2;j<mh;j+=2){
   k+=4;
   wk1r=w[k];
   wk1i=w[k+1];
   wk3r=w[k+2];
   wk3i=w[k+3];
   kr-=4;
   wd1i=w[kr];
   wd1r=w[kr+1];
   wd3i=w[kr+2];
   wd3r=w[kr+3];
   j1=j+m;
   j2=j1+m;
   j3=j2+m;
   x0r=a[j]-a[j2+1];
   x0i=a[j+1]+a[j2];
   x1r=a[j]+a[j2+1];
   x1i=a[j+1]-a[j2];
   x2r=a[j1]-a[j3+1];
   x2i=a[j1+1]+a[j3];
   x3r=a[j1]+a[j3+1];
   x3i=a[j1+1]-a[j3];
   y0r=wk1r*x0r-wk1i*x0i;
   y0i=wk1r*x0i+wk1i*x0r;
   y2r=wd1r*x2r-wd1i*x2i;
   y2i=wd1r*x2i+wd1i*x2r;
   a[j]=y0r+y2r;
   a[j+1]=y0i+y2i;
   a[j1]=y0r-y2r;
   a[j1+1]=y0i-y2i;
   y0r=wk3r*x1r+wk3i*x1i;
   y0i=wk3r*x1i-wk3i*x1r;
   y2r=wd3r*x3r+wd3i*x3i;
   y2i=wd3r*x3i-wd3i*x3r;
   a[j2]=y0r+y2r;
   a[j2+1]=y0i+y2i;
   a[j3]=y0r-y2r;
   a[j3+1]=y0i-y2i;
   j0=m-j;
   j1=j0+m;
   j2=j1+m;
   j3=j2+m;
   x0r=a[j0]-a[j2+1];
   x0i=a[j0+1]+a[j2];
   x1r=a[j0]+a[j2+1];
   x1i=a[j0+1]-a[j2];
   x2r=a[j1]-a[j3+1];
   x2i=a[j1+1]+a[j3];
   x3r=a[j1]+a[j3+1];
   x3i=a[j1+1]-a[j3];
   y0r=wd1i*x0r-wd1r*x0i;
   y0i=wd1i*x0i+wd1r*x0r;
   y2r=wk1i*x2r-wk1r*x2i;
   y2i=wk1i*x2i+wk1r*x2r;
   a[j0]=y0r+y2r;
   a[j0+1]=y0i+y2i;
   a[j1]=y0r-y2r;
   a[j1+1]=y0i-y2i;
   y0r=wd3i*x1r+wd3r*x1i;
   y0i=wd3i*x1i-wd3r*x1r;
   y2r=wk3i*x3r+wk3r*x3i;
   y2i=wk3i*x3i-wk3r*x3r;
   a[j2]=y0r+y2r;
   a[j2+1]=y0i+y2i;
   a[j3]=y0r-y2r;
   a[j3+1]=y0i-y2i;
}
wk1r=w[m];
wk1i=w[m+1];
j0=mh;
j1=j0+m;
j2=j1+m;
j3=j2+m;
x0r=a[j0]-a[j2+1];
x0i=a[j0+1]+a[j2];
x1r=a[j0]+a[j2+1];
x1i=a[j0+1]-a[j2];
x2r=a[j1]-a[j3+1];
x2i=a[j1+1]+a[j3];
x3r=a[j1]+a[j3+1];
x3i=a[j1+1]-a[j3];
y0r=wk1r*x0r-wk1i*x0i;
y0i=wk1r*x0i+wk1i*x0r;
y2r=wk1i*x2r-wk1r*x2i;
y2i=wk1i*x2i+wk1r*x2r;
a[j0]=y0r+y2r;
a[j0+1]=y0i+y2i;
a[j1]=y0r-y2r;
a[j1+1]=y0i-y2i;
y0r=wk1i*x1r-wk1r*x1i;
y0i=wk1i*x1i+wk1r*x1r;
y2r=wk1r*x3r-wk1i*x3i;
y2i=wk1r*x3i+wk1i*x3r;
a[j2]=y0r-y2r;
a[j2+1]=y0i-y2i;
a[j3]=y0r+y2r;
a[j3+1]=y0i+y2i;
}

template <class T> void cftf161(T *a,T *w){
T wn4r,wk1r,wk1i,
  x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,
  y0r,y0i,y1r,y1i,y2r,y2i,y3r,y3i,
  y4r,y4i,y5r,y5i,y6r,y6i,y7r,y7i,
  y8r,y8i,y9r,y9i,y10r,y10i,y11r,y11i,
  y12r,y12i,y13r,y13i,y14r,y14i,y15r,y15i;
    
wn4r=w[1];
wk1r=w[2];
wk1i=w[3];
x0r=a[0]+a[16];
x0i=a[1]+a[17];
x1r=a[0]-a[16];
x1i=a[1]-a[17];
x2r=a[8]+a[24];
x2i=a[9]+a[25];
x3r=a[8]-a[24];
x3i=a[9]-a[25];
y0r=x0r+x2r;
y0i=x0i+x2i;
y4r=x0r-x2r;
y4i=x0i-x2i;
y8r=x1r-x3i;
y8i=x1i+x3r;
y12r=x1r+x3i;
y12i=x1i-x3r;
x0r=a[2]+a[18];
x0i=a[3]+a[19];
x1r=a[2]-a[18];
x1i=a[3]-a[19];
x2r=a[10]+a[26];
x2i=a[11]+a[27];
x3r=a[10]-a[26];
x3i=a[11]-a[27];
y1r=x0r+x2r;
y1i=x0i+x2i;
y5r=x0r-x2r;
y5i=x0i-x2i;
x0r=x1r-x3i;
x0i=x1i+x3r;
y9r=wk1r*x0r-wk1i*x0i;
y9i=wk1r*x0i+wk1i*x0r;
x0r=x1r+x3i;
x0i=x1i-x3r;
y13r=wk1i*x0r-wk1r*x0i;
y13i=wk1i*x0i+wk1r*x0r;
x0r=a[4]+a[20];
x0i=a[5]+a[21];
x1r=a[4]-a[20];
x1i=a[5]-a[21];
x2r=a[12]+a[28];
x2i=a[13]+a[29];
x3r=a[12]-a[28];
x3i=a[13]-a[29];
y2r=x0r+x2r;
y2i=x0i+x2i;
y6r=x0r-x2r;
y6i=x0i-x2i;
x0r=x1r-x3i;
x0i=x1i+x3r;
y10r=wn4r*(x0r-x0i);
y10i=wn4r*(x0i+x0r);
x0r=x1r+x3i;
x0i=x1i-x3r;
y14r=wn4r*(x0r+x0i);
y14i=wn4r*(x0i-x0r);
x0r=a[6]+a[22];
x0i=a[7]+a[23];
x1r=a[6]-a[22];
x1i=a[7]-a[23];
x2r=a[14]+a[30];
x2i=a[15]+a[31];
x3r=a[14]-a[30];
x3i=a[15]-a[31];
y3r=x0r+x2r;
y3i=x0i+x2i;
y7r=x0r-x2r;
y7i=x0i-x2i;
x0r=x1r-x3i;
x0i=x1i+x3r;
y11r=wk1i*x0r-wk1r*x0i;
y11i=wk1i*x0i+wk1r*x0r;
x0r=x1r+x3i;
x0i=x1i-x3r;
y15r=wk1r*x0r-wk1i*x0i;
y15i=wk1r*x0i+wk1i*x0r;
x0r=y12r-y14r;
x0i=y12i-y14i;
x1r=y12r+y14r;
x1i=y12i+y14i;
x2r=y13r-y15r;
x2i=y13i-y15i;
x3r=y13r+y15r;
x3i=y13i+y15i;
a[24]=x0r+x2r;
a[25]=x0i+x2i;
a[26]=x0r-x2r;
a[27]=x0i-x2i;
a[28]=x1r-x3i;
a[29]=x1i+x3r;
a[30]=x1r+x3i;
a[31]=x1i-x3r;
x0r=y8r+y10r;
x0i=y8i+y10i;
x1r=y8r-y10r;
x1i=y8i-y10i;
x2r=y9r+y11r;
x2i=y9i+y11i;
x3r=y9r-y11r;
x3i=y9i-y11i;
a[16]=x0r+x2r;
a[17]=x0i+x2i;
a[18]=x0r-x2r;
a[19]=x0i-x2i;
a[20]=x1r-x3i;
a[21]=x1i+x3r;
a[22]=x1r+x3i;
a[23]=x1i-x3r;
x0r=y5r-y7i;
x0i=y5i+y7r;
x2r=wn4r*(x0r-x0i);
x2i=wn4r*(x0i+x0r);
x0r=y5r+y7i;
x0i=y5i-y7r;
x3r=wn4r*(x0r-x0i);
x3i=wn4r*(x0i+x0r);
x0r=y4r-y6i;
x0i=y4i+y6r;
x1r=y4r+y6i;
x1i=y4i-y6r;
a[8]=x0r+x2r;
a[9]=x0i+x2i;
a[10]=x0r-x2r;
a[11]=x0i-x2i;
a[12]=x1r-x3i;
a[13]=x1i+x3r;
a[14]=x1r+x3i;
a[15]=x1i-x3r;
x0r=y0r+y2r;
x0i=y0i+y2i;
x1r=y0r-y2r;
x1i=y0i-y2i;
x2r=y1r+y3r;
x2i=y1i+y3i;
x3r=y1r-y3r;
x3i=y1i-y3i;
a[0]=x0r+x2r;
a[1]=x0i+x2i;
a[2]=x0r-x2r;
a[3]=x0i-x2i;
a[4]=x1r-x3i;
a[5]=x1i+x3r;
a[6]=x1r+x3i;
a[7]=x1i-x3r;
}


template <class T> void cftf162(T *a,T *w){
T wn4r,wk1r,wk1i,wk2r,wk2i,wk3r,wk3i,
  x0r,x0i,x1r,x1i,x2r,x2i,
  y0r,y0i,y1r,y1i,y2r,y2i,y3r,y3i,
  y4r,y4i,y5r,y5i,y6r,y6i,y7r,y7i,
  y8r,y8i,y9r,y9i,y10r,y10i,y11r,y11i,
  y12r,y12i,y13r,y13i,y14r,y14i,y15r,y15i;
    
wn4r=w[1];
wk1r=w[4];
wk1i=w[5];
wk3r=w[6];
wk3i=-w[7];
wk2r=w[8];
wk2i=w[9];
x1r=a[0]-a[17];
x1i=a[1]+a[16];
x0r=a[8]-a[25];
x0i=a[9]+a[24];
x2r=wn4r*(x0r-x0i);
x2i=wn4r*(x0i+x0r);
y0r=x1r+x2r;
y0i=x1i+x2i;
y4r=x1r-x2r;
y4i=x1i-x2i;
x1r=a[0]+a[17];
x1i=a[1]-a[16];
x0r=a[8]+a[25];
x0i=a[9]-a[24];
x2r=wn4r*(x0r-x0i);
x2i=wn4r*(x0i+x0r);
y8r=x1r-x2i;
y8i=x1i+x2r;
y12r=x1r+x2i;
y12i=x1i-x2r;
x0r=a[2]-a[19];
x0i=a[3]+a[18];
x1r=wk1r*x0r-wk1i*x0i;
x1i=wk1r*x0i+wk1i*x0r;
x0r=a[10]-a[27];
x0i=a[11]+a[26];
x2r=wk3i*x0r-wk3r*x0i;
x2i=wk3i*x0i+wk3r*x0r;
y1r=x1r+x2r;
y1i=x1i+x2i;
y5r=x1r-x2r;
y5i=x1i-x2i;
x0r=a[2]+a[19];
x0i=a[3]-a[18];
x1r=wk3r*x0r-wk3i*x0i;
x1i=wk3r*x0i+wk3i*x0r;
x0r=a[10]+a[27];
x0i=a[11]-a[26];
x2r=wk1r*x0r+wk1i*x0i;
x2i=wk1r*x0i-wk1i*x0r;
y9r=x1r-x2r;
y9i=x1i-x2i;
y13r=x1r+x2r;
y13i=x1i+x2i;
x0r=a[4]-a[21];
x0i=a[5]+a[20];
x1r=wk2r*x0r-wk2i*x0i;
x1i=wk2r*x0i+wk2i*x0r;
x0r=a[12]-a[29];
x0i=a[13]+a[28];
x2r=wk2i*x0r-wk2r*x0i;
x2i=wk2i*x0i+wk2r*x0r;
y2r=x1r+x2r;
y2i=x1i+x2i;
y6r=x1r-x2r;
y6i=x1i-x2i;
x0r=a[4]+a[21];
x0i=a[5]-a[20];
x1r=wk2i*x0r-wk2r*x0i;
x1i=wk2i*x0i+wk2r*x0r;
x0r=a[12]+a[29];
x0i=a[13]-a[28];
x2r=wk2r*x0r-wk2i*x0i;
x2i=wk2r*x0i+wk2i*x0r;
y10r=x1r-x2r;
y10i=x1i-x2i;
y14r=x1r+x2r;
y14i=x1i+x2i;
x0r=a[6]-a[23];
x0i=a[7]+a[22];
x1r=wk3r*x0r-wk3i*x0i;
x1i=wk3r*x0i+wk3i*x0r;
x0r=a[14]-a[31];
x0i=a[15]+a[30];
x2r=wk1i*x0r-wk1r*x0i;
x2i=wk1i*x0i+wk1r*x0r;
y3r=x1r+x2r;
y3i=x1i+x2i;
y7r=x1r-x2r;
y7i=x1i-x2i;
x0r=a[6]+a[23];
x0i=a[7]-a[22];
x1r=wk1i*x0r+wk1r*x0i;
x1i=wk1i*x0i-wk1r*x0r;
x0r=a[14]+a[31];
x0i=a[15]-a[30];
x2r=wk3i*x0r-wk3r*x0i;
x2i=wk3i*x0i+wk3r*x0r;
y11r=x1r+x2r;
y11i=x1i+x2i;
y15r=x1r-x2r;
y15i=x1i-x2i;
x1r=y0r+y2r;
x1i=y0i+y2i;
x2r=y1r+y3r;
x2i=y1i+y3i;
a[0]=x1r+x2r;
a[1]=x1i+x2i;
a[2]=x1r-x2r;
a[3]=x1i-x2i;
x1r=y0r-y2r;
x1i=y0i-y2i;
x2r=y1r-y3r;
x2i=y1i-y3i;
a[4]=x1r-x2i;
a[5]=x1i+x2r;
a[6]=x1r+x2i;
a[7]=x1i-x2r;
x1r=y4r-y6i;
x1i=y4i+y6r;
x0r=y5r-y7i;
x0i=y5i+y7r;
x2r=wn4r*(x0r-x0i);
x2i=wn4r*(x0i+x0r);
a[8]=x1r+x2r;
a[9]=x1i+x2i;
a[10]=x1r-x2r;
a[11]=x1i-x2i;
x1r=y4r+y6i;
x1i=y4i-y6r;
x0r=y5r+y7i;
x0i=y5i-y7r;
x2r=wn4r*(x0r-x0i);
x2i=wn4r*(x0i+x0r);
a[12]=x1r-x2i;
a[13]=x1i+x2r;
a[14]=x1r+x2i;
a[15]=x1i-x2r;
x1r=y8r+y10r;
x1i=y8i+y10i;
x2r=y9r-y11r;
x2i=y9i-y11i;
a[16]=x1r+x2r;
a[17]=x1i+x2i;
a[18]=x1r-x2r;
a[19]=x1i-x2i;
x1r=y8r-y10r;
x1i=y8i-y10i;
x2r=y9r+y11r;
x2i=y9i+y11i;
a[20]=x1r-x2i;
a[21]=x1i+x2r;
a[22]=x1r+x2i;
a[23]=x1i-x2r;
x1r=y12r-y14i;
x1i=y12i+y14r;
x0r=y13r+y15i;
x0i=y13i-y15r;
x2r=wn4r*(x0r-x0i);
x2i=wn4r*(x0i+x0r);
a[24]=x1r+x2r;
a[25]=x1i+x2i;
a[26]=x1r-x2r;
a[27]=x1i-x2i;
x1r=y12r+y14i;
x1i=y12i-y14r;
x0r=y13r-y15i;
x0i=y13i+y15r;
x2r=wn4r*(x0r-x0i);
x2i=wn4r*(x0i+x0r);
a[28]=x1r-x2i;
a[29]=x1i+x2r;
a[30]=x1r+x2i;
a[31]=x1i-x2r;
}


template <class T> void cftf081(T *a,T *w){
T wn4r,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,
  y0r,y0i,y1r,y1i,y2r,y2i,y3r,y3i,
  y4r,y4i,y5r,y5i,y6r,y6i,y7r,y7i;
    
wn4r=w[1];
x0r=a[0]+a[8];
x0i=a[1]+a[9];
x1r=a[0]-a[8];
x1i=a[1]-a[9];
x2r=a[4]+a[12];
x2i=a[5]+a[13];
x3r=a[4]-a[12];
x3i=a[5]-a[13];
y0r=x0r+x2r;
y0i=x0i+x2i;
y2r=x0r-x2r;
y2i=x0i-x2i;
y1r=x1r-x3i;
y1i=x1i+x3r;
y3r=x1r+x3i;
y3i=x1i-x3r;
x0r=a[2]+a[10];
x0i=a[3]+a[11];
x1r=a[2]-a[10];
x1i=a[3]-a[11];
x2r=a[6]+a[14];
x2i=a[7]+a[15];
x3r=a[6]-a[14];
x3i=a[7]-a[15];
y4r=x0r+x2r;
y4i=x0i+x2i;
y6r=x0r-x2r;
y6i=x0i-x2i;
x0r=x1r-x3i;
x0i=x1i+x3r;
x2r=x1r+x3i;
x2i=x1i-x3r;
y5r=wn4r*(x0r-x0i);
y5i=wn4r*(x0r+x0i);
y7r=wn4r*(x2r-x2i);
y7i=wn4r*(x2r+x2i);
a[8]=y1r+y5r;
a[9]=y1i+y5i;
a[10]=y1r-y5r;
a[11]=y1i-y5i;
a[12]=y3r-y7i;
a[13]=y3i+y7r;
a[14]=y3r+y7i;
a[15]=y3i-y7r;
a[0]=y0r+y4r;
a[1]=y0i+y4i;
a[2]=y0r-y4r;
a[3]=y0i-y4i;
a[4]=y2r-y6i;
a[5]=y2i+y6r;
a[6]=y2r+y6i;
a[7]=y2i-y6r;
}


template <class T> void cftf082(T *a,T *w){
T wn4r,wk1r,wk1i,x0r,x0i,x1r,x1i,
  y0r,y0i,y1r,y1i,y2r,y2i,y3r,y3i,
  y4r,y4i,y5r,y5i,y6r,y6i,y7r,y7i;
    
wn4r=w[1];
wk1r=w[2];
wk1i=w[3];
y0r=a[0]-a[9];
y0i=a[1]+a[8];
y1r=a[0]+a[9];
y1i=a[1]-a[8];
x0r=a[4]-a[13];
x0i=a[5]+a[12];
y2r=wn4r*(x0r-x0i);
y2i=wn4r*(x0i+x0r);
x0r=a[4]+a[13];
x0i=a[5]-a[12];
y3r=wn4r*(x0r-x0i);
y3i=wn4r*(x0i+x0r);
x0r=a[2]-a[11];
x0i=a[3]+a[10];
y4r=wk1r*x0r-wk1i*x0i;
y4i=wk1r*x0i+wk1i*x0r;
x0r=a[2]+a[11];
x0i=a[3]-a[10];
y5r=wk1i*x0r-wk1r*x0i;
y5i=wk1i*x0i+wk1r*x0r;
x0r=a[6]-a[15];
x0i=a[7]+a[14];
y6r=wk1i*x0r-wk1r*x0i;
y6i=wk1i*x0i+wk1r*x0r;
x0r=a[6]+a[15];
x0i=a[7]-a[14];
y7r=wk1r*x0r-wk1i*x0i;
y7i=wk1r*x0i+wk1i*x0r;
x0r=y0r+y2r;
x0i=y0i+y2i;
x1r=y4r+y6r;
x1i=y4i+y6i;
a[0]=x0r+x1r;
a[1]=x0i+x1i;
a[2]=x0r-x1r;
a[3]=x0i-x1i;
x0r=y0r-y2r;
x0i=y0i-y2i;
x1r=y4r-y6r;
x1i=y4i-y6i;
a[4]=x0r-x1i;
a[5]=x0i+x1r;
a[6]=x0r+x1i;
a[7]=x0i-x1r;
x0r=y1r-y3i;
x0i=y1i+y3r;
x1r=y5r-y7r;
x1i=y5i-y7i;
a[8]=x0r+x1r;
a[9]=x0i+x1i;
a[10]=x0r-x1r;
a[11]=x0i-x1i;
x0r=y1r+y3i;
x0i=y1i-y3r;
x1r=y5r+y7r;
x1i=y5i+y7i;
a[12]=x0r-x1i;
a[13]=x0i+x1r;
a[14]=x0r+x1i;
a[15]=x0i-x1r;
}


template <class T> void cftf040(T *a){
T x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i;
    
x0r=a[0]+a[4];
x0i=a[1]+a[5];
x1r=a[0]-a[4];
x1i=a[1]-a[5];
x2r=a[2]+a[6];
x2i=a[3]+a[7];
x3r=a[2]-a[6];
x3i=a[3]-a[7];
a[0]=x0r+x2r;
a[1]=x0i+x2i;
a[2]=x1r-x3i;
a[3]=x1i+x3r;
a[4]=x0r-x2r;
a[5]=x0i-x2i;
a[6]=x1r+x3i;
a[7]=x1i-x3r;
}


template <class T> void cftb040(T *a){
T x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i;
    
x0r=a[0]+a[4];
x0i=a[1]+a[5];
x1r=a[0]-a[4];
x1i=a[1]-a[5];
x2r=a[2]+a[6];
x2i=a[3]+a[7];
x3r=a[2]-a[6];
x3i=a[3]-a[7];
a[0]=x0r+x2r;
a[1]=x0i+x2i;
a[2]=x1r+x3i;
a[3]=x1i-x3r;
a[4]=x0r-x2r;
a[5]=x0i-x2i;
a[6]=x1r-x3i;
a[7]=x1i+x3r;
}


template <class T> void cftx020(T *a){
T x0r,x0i;
    
x0r=a[0]-a[2];
x0i=a[1]-a[3];
a[0]+=a[2];
a[1]+=a[3];
a[2]=x0r;
a[3]=x0i;
}


template <class T> void rftfsub(int n,T *a,int nc,T *c){
int j,k,kk,ks,m; 
T wkr,wki,xr,xi,yr,yi;
    
m=n>>1;
ks=2*nc/m;
kk=0;
for(j=2;j<m;j+=2){
   k=n-j;
   kk+=ks;
   wkr=0.5-c[nc-kk];
   wki=c[kk];
   xr=a[j]-a[k];
   xi=a[j+1]+a[k+1];
   yr=wkr*xr-wki*xi;
   yi=wkr*xi+wki*xr;
   a[j]-=yr;
   a[j+1]-=yi;
   a[k]+=yr;
   a[k+1]-=yi;
}
}


template <class T> void rftbsub(int n,T *a,int nc,T *c){
int j,k,kk,ks,m;
T wkr,wki,xr,xi,yr,yi;
    
m=n>>1;
ks=2*nc/m;
kk=0;
for(j=2;j<m;j+=2){
   k=n-j;
   kk+=ks;
   wkr=0.5-c[nc-kk];
   wki=c[kk];
   xr=a[j]-a[k];
   xi=a[j+1]+a[k+1];
   yr=wkr*xr+wki*xi;
   yi=wkr*xi-wki*xr;
   a[j]-=yr;
   a[j+1]-=yi;
   a[k]+=yr;
   a[k+1]-=yi;
}
}


template <class T> void dctsub(int n,T *a,int nc,T *c){
int j,k,kk,ks,m;
T wkr,wki,xr;
    
m=n>>1;
ks=nc/n;
kk=0;
for(j=1;j<m;j++){
   k=n-j;
   kk+=ks;
   wkr=c[kk]-c[nc-kk];
   wki=c[kk]+c[nc-kk];
   xr=wki*a[j]-wkr*a[k];
   a[j]=wkr*a[j]+wki*a[k];
   a[k]=xr;
}
a[m]*=c[0];
}


template <class T> void dstsub(int n,T *a,int nc,T *c){
int j,k,kk,ks,m;
T wkr,wki,xr;
    
m=n>>1;
ks=nc/n;
kk=0;
for(j=1;j<m;j++){
   k=n-j;
   kk+=ks;
   wkr=c[kk]-c[nc-kk];
   wki=c[kk]+c[nc-kk];
   xr=wki*a[k]-wkr*a[j];
   a[k]=wkr*a[k]+wki*a[j];
   a[j]=xr;
}
a[m]*=c[0];
}


template <class T> void cftf1st(int n,T *a,T *w){
int j,j0,j1,j2,j3,k,m,mh;
T wn4r,csc1,csc3,wk1r,wk1i,wk3r,wk3i,wd1r,wd1i,wd3r,wd3i;
T x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,y0r,y0i,y1r,y1i,y2r,y2i,y3r,y3i;
    
mh=n>>3;
m=2*mh;
j1=m;
j2=j1+m;
j3=j2+m;
x0r=a[0]+a[j2];
x0i=a[1]+a[j2+1];
x1r=a[0]-a[j2];
x1i=a[1]-a[j2+1];
x2r=a[j1]+a[j3];
x2i=a[j1+1]+a[j3+1];
x3r=a[j1]-a[j3];
x3i=a[j1+1]-a[j3+1];
a[0]=x0r+x2r;
a[1]=x0i+x2i;
a[j1]=x0r-x2r;
a[j1+1]=x0i-x2i;
a[j2]=x1r-x3i;
a[j2+1]=x1i+x3r;
a[j3]=x1r+x3i;
a[j3+1]=x1i-x3r;
wn4r=w[1];
csc1=w[2];
csc3=w[3];
wd1r=1;
wd1i=0;
wd3r=1;
wd3i=0;
k=0;
for(j=2;j<mh-2;j+=4){
   k+=4;
   wk1r=csc1*(wd1r+w[k]);
   wk1i=csc1*(wd1i+w[k+1]);
   wk3r=csc3*(wd3r+w[k+2]);
   wk3i=csc3*(wd3i+w[k+3]);
   wd1r=w[k];
   wd1i=w[k+1];
   wd3r=w[k+2];
   wd3i=w[k+3];
   j1=j+m;
   j2=j1+m;
   j3=j2+m;
   x0r=a[j]+a[j2];
   x0i=a[j+1]+a[j2+1];
   x1r=a[j]-a[j2];
   x1i=a[j+1]-a[j2+1];
   y0r=a[j+2]+a[j2+2];
   y0i=a[j+3]+a[j2+3];
   y1r=a[j+2]-a[j2+2];
   y1i=a[j+3]-a[j2+3];
   x2r=a[j1]+a[j3];
   x2i=a[j1+1]+a[j3+1];
   x3r=a[j1]-a[j3];
   x3i=a[j1+1]-a[j3+1];
   y2r=a[j1+2]+a[j3+2];
   y2i=a[j1+3]+a[j3+3];
   y3r=a[j1+2]-a[j3+2];
   y3i=a[j1+3]-a[j3+3];
   a[j]=x0r+x2r;
   a[j+1]=x0i+x2i;
   a[j+2]=y0r+y2r;
   a[j+3]=y0i+y2i;
   a[j1]=x0r-x2r;
   a[j1+1]=x0i-x2i;
   a[j1+2]=y0r-y2r;
   a[j1+3]=y0i-y2i;
   x0r=x1r-x3i;
   x0i=x1i+x3r;
   a[j2]=wk1r*x0r-wk1i*x0i;
   a[j2+1]=wk1r*x0i+wk1i*x0r;
   x0r=y1r-y3i;
   x0i=y1i+y3r;
   a[j2+2]=wd1r*x0r-wd1i*x0i;
   a[j2+3]=wd1r*x0i+wd1i*x0r;
   x0r=x1r+x3i;
   x0i=x1i-x3r;
   a[j3]=wk3r*x0r+wk3i*x0i;
   a[j3+1]=wk3r*x0i-wk3i*x0r;
   x0r=y1r+y3i;
   x0i=y1i-y3r;
   a[j3+2]=wd3r*x0r+wd3i*x0i;
   a[j3+3]=wd3r*x0i-wd3i*x0r;
   j0=m-j;
   j1=j0+m;
   j2=j1+m;
   j3=j2+m;
   x0r=a[j0]+a[j2];
   x0i=a[j0+1]+a[j2+1];
   x1r=a[j0]-a[j2];
   x1i=a[j0+1]-a[j2+1];
   y0r=a[j0-2]+a[j2-2];
   y0i=a[j0-1]+a[j2-1];
   y1r=a[j0-2]-a[j2-2];
   y1i=a[j0-1]-a[j2-1];
   x2r=a[j1]+a[j3];
   x2i=a[j1+1]+a[j3+1];
   x3r=a[j1]-a[j3];
   x3i=a[j1+1]-a[j3+1];
   y2r=a[j1-2]+a[j3-2];
   y2i=a[j1-1]+a[j3-1];
   y3r=a[j1-2]-a[j3-2];
   y3i=a[j1-1]-a[j3-1];
   a[j0]=x0r+x2r;
   a[j0+1]=x0i+x2i;
   a[j0-2]=y0r+y2r;
   a[j0-1]=y0i+y2i;
   a[j1]=x0r-x2r;
   a[j1+1]=x0i-x2i;
   a[j1-2]=y0r-y2r;
   a[j1-1]=y0i-y2i;
   x0r=x1r-x3i;
   x0i=x1i+x3r;
   a[j2]=wk1i*x0r-wk1r*x0i;
   a[j2+1]=wk1i*x0i+wk1r*x0r;
   x0r=y1r-y3i;
   x0i=y1i+y3r;
   a[j2-2]=wd1i*x0r-wd1r*x0i;
   a[j2-1]=wd1i*x0i+wd1r*x0r;
   x0r=x1r+x3i;
   x0i=x1i-x3r;
   a[j3]=wk3i*x0r+wk3r*x0i;
   a[j3+1]=wk3i*x0i-wk3r*x0r;
   x0r=y1r+y3i;
   x0i=y1i-y3r;
   a[j3-2]=wd3i*x0r+wd3r*x0i;
   a[j3-1]=wd3i*x0i-wd3r*x0r;
}

wk1r=csc1*(wd1r+wn4r);
wk1i=csc1*(wd1i+wn4r);
wk3r=csc3*(wd3r-wn4r);
wk3i=csc3*(wd3i-wn4r);
j0=mh;
j1=j0+m;
j2=j1+m;
j3=j2+m;
x0r=a[j0-2]+a[j2-2];
x0i=a[j0-1]+a[j2-1];
x1r=a[j0-2]-a[j2-2];
x1i=a[j0-1]-a[j2-1];
x2r=a[j1-2]+a[j3-2];
x2i=a[j1-1]+a[j3-1];
x3r=a[j1-2]-a[j3-2];
x3i=a[j1-1]-a[j3-1];
a[j0-2]=x0r+x2r;
a[j0-1]=x0i+x2i;
a[j1-2]=x0r-x2r;
a[j1-1]=x0i-x2i;
x0r=x1r-x3i;
x0i=x1i+x3r;
a[j2-2]=wk1r*x0r-wk1i*x0i;
a[j2-1]=wk1r*x0i+wk1i*x0r;
x0r=x1r+x3i;
x0i=x1i-x3r;
a[j3-2]=wk3r*x0r+wk3i*x0i;
a[j3-1]=wk3r*x0i-wk3i*x0r;
x0r=a[j0]+a[j2];
x0i=a[j0+1]+a[j2+1];
x1r=a[j0]-a[j2];
x1i=a[j0+1]-a[j2+1];
x2r=a[j1]+a[j3];
x2i=a[j1+1]+a[j3+1];
x3r=a[j1]-a[j3];
x3i=a[j1+1]-a[j3+1];
a[j0]=x0r+x2r;
a[j0+1]=x0i+x2i;
a[j1]=x0r-x2r;
a[j1+1]=x0i-x2i;
x0r=x1r-x3i;
x0i=x1i+x3r;
a[j2]=wn4r*(x0r-x0i);
a[j2+1]=wn4r*(x0i+x0r);
x0r=x1r+x3i;
x0i=x1i-x3r;
a[j3]=-wn4r*(x0r+x0i);
a[j3+1]=-wn4r*(x0i-x0r);
x0r=a[j0+2]+a[j2+2];
x0i=a[j0+3]+a[j2+3];
x1r=a[j0+2]-a[j2+2];
x1i=a[j0+3]-a[j2+3];
x2r=a[j1+2]+a[j3+2];
x2i=a[j1+3]+a[j3+3];
x3r=a[j1+2]-a[j3+2];
x3i=a[j1+3]-a[j3+3];
a[j0+2]=x0r+x2r;
a[j0+3]=x0i+x2i;
a[j1+2]=x0r-x2r;
a[j1+3]=x0i-x2i;
x0r=x1r-x3i;
x0i=x1i+x3r;
a[j2+2]=wk1i*x0r-wk1r*x0i;
a[j2+3]=wk1i*x0i+wk1r*x0r;
x0r=x1r+x3i;
x0i=x1i-x3r;
a[j3+2]=wk3i*x0r+wk3r*x0i;
a[j3+3]=wk3i*x0i-wk3r*x0r;

}


template <class T> void cftb1st(int n,T *a,T *w){
int j,j0,j1,j2,j3,k,m,mh;
T wn4r,csc1,csc3,wk1r,wk1i,wk3r,wk3i,wd1r,wd1i,wd3r,wd3i;
T x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,y0r,y0i,y1r,y1i,y2r,y2i,y3r,y3i;
    
mh=n>>3;
m=2*mh;
j1=m;
j2=j1+m;
j3=j2+m;
x0r=a[0]+a[j2];
x0i=-a[1]-a[j2+1];
x1r=a[0]-a[j2];
x1i=-a[1]+a[j2+1];
x2r=a[j1]+a[j3];
x2i=a[j1+1]+a[j3+1];
x3r=a[j1]-a[j3];
x3i=a[j1+1]-a[j3+1];
a[0]=x0r+x2r;
a[1]=x0i-x2i;
a[j1]=x0r-x2r;
a[j1+1]=x0i+x2i;
a[j2]=x1r+x3i;
a[j2+1]=x1i+x3r;
a[j3]=x1r-x3i;
a[j3+1]=x1i-x3r;
wn4r=w[1];
csc1=w[2];
csc3=w[3];
wd1r=1;
wd1i=0;
wd3r=1;
wd3i=0;
k=0;
for(j=2;j<mh-2;j+=4){
   k+=4;
   wk1r=csc1*(wd1r+w[k]);
   wk1i=csc1*(wd1i+w[k+1]);
   wk3r=csc3*(wd3r+w[k+2]);
   wk3i=csc3*(wd3i+w[k+3]);
   wd1r=w[k];
   wd1i=w[k+1];
   wd3r=w[k+2];
   wd3i=w[k+3];
   j1=j+m;
   j2=j1+m;
   j3=j2+m;
   x0r=a[j]+a[j2];
   x0i=-a[j+1]-a[j2+1];
   x1r=a[j]-a[j2];
   x1i=-a[j+1]+a[j2+1];
   y0r=a[j+2]+a[j2+2];
   y0i=-a[j+3]-a[j2+3];
   y1r=a[j+2]-a[j2+2];
   y1i=-a[j+3]+a[j2+3];
   x2r=a[j1]+a[j3];
   x2i=a[j1+1]+a[j3+1];
   x3r=a[j1]-a[j3];
   x3i=a[j1+1]-a[j3+1];
   y2r=a[j1+2]+a[j3+2];
   y2i=a[j1+3]+a[j3+3];
   y3r=a[j1+2]-a[j3+2];
   y3i=a[j1+3]-a[j3+3];
   a[j]=x0r+x2r;
   a[j+1]=x0i-x2i;
   a[j+2]=y0r+y2r;
   a[j+3]=y0i-y2i;
   a[j1]=x0r-x2r;
   a[j1+1]=x0i+x2i;
   a[j1+2]=y0r-y2r;
   a[j1+3]=y0i+y2i;
   x0r=x1r+x3i;
   x0i=x1i+x3r;
   a[j2]=wk1r*x0r-wk1i*x0i;
   a[j2+1]=wk1r*x0i+wk1i*x0r;
   x0r=y1r+y3i;
   x0i=y1i+y3r;
   a[j2+2]=wd1r*x0r-wd1i*x0i;
   a[j2+3]=wd1r*x0i+wd1i*x0r;
   x0r=x1r-x3i;
   x0i=x1i-x3r;
   a[j3]=wk3r*x0r+wk3i*x0i;
   a[j3+1]=wk3r*x0i-wk3i*x0r;
   x0r=y1r-y3i;
   x0i=y1i-y3r;
   a[j3+2]=wd3r*x0r+wd3i*x0i;
   a[j3+3]=wd3r*x0i-wd3i*x0r;
   j0=m-j;
   j1=j0+m;
   j2=j1+m;
   j3=j2+m;
   x0r=a[j0]+a[j2];
   x0i=-a[j0+1]-a[j2+1];
   x1r=a[j0]-a[j2];
   x1i=-a[j0+1]+a[j2+1];
   y0r=a[j0-2]+a[j2-2];
   y0i=-a[j0-1]-a[j2-1];
   y1r=a[j0-2]-a[j2-2];
   y1i=-a[j0-1]+a[j2-1];
   x2r=a[j1]+a[j3];
   x2i=a[j1+1]+a[j3+1];
   x3r=a[j1]-a[j3];
   x3i=a[j1+1]-a[j3+1];
   y2r=a[j1-2]+a[j3-2];
   y2i=a[j1-1]+a[j3-1];
   y3r=a[j1-2]-a[j3-2];
   y3i=a[j1-1]-a[j3-1];
   a[j0]=x0r+x2r;
   a[j0+1]=x0i-x2i;
   a[j0-2]=y0r+y2r;
   a[j0-1]=y0i-y2i;
   a[j1]=x0r-x2r;
   a[j1+1]=x0i+x2i;
   a[j1-2]=y0r-y2r;
   a[j1-1]=y0i+y2i;
   x0r=x1r+x3i;
   x0i=x1i+x3r;
   a[j2]=wk1i*x0r-wk1r*x0i;
   a[j2+1]=wk1i*x0i+wk1r*x0r;
   x0r=y1r+y3i;
   x0i=y1i+y3r;
   a[j2-2]=wd1i*x0r-wd1r*x0i;
   a[j2-1]=wd1i*x0i+wd1r*x0r;
   x0r=x1r-x3i;
   x0i=x1i-x3r;
   a[j3]=wk3i*x0r+wk3r*x0i;
   a[j3+1]=wk3i*x0i-wk3r*x0r;
   x0r=y1r-y3i;
   x0i=y1i-y3r;
   a[j3-2]=wd3i*x0r+wd3r*x0i;
   a[j3-1]=wd3i*x0i-wd3r*x0r;
}
wk1r=csc1*(wd1r+wn4r);
wk1i=csc1*(wd1i+wn4r);
wk3r=csc3*(wd3r-wn4r);
wk3i=csc3*(wd3i-wn4r);
j0=mh;
j1=j0+m;
j2=j1+m;
j3=j2+m;
x0r=a[j0-2]+a[j2-2];
x0i=-a[j0-1]-a[j2-1];
x1r=a[j0-2]-a[j2-2];
x1i=-a[j0-1]+a[j2-1];
x2r=a[j1-2]+a[j3-2];
x2i=a[j1-1]+a[j3-1];
x3r=a[j1-2]-a[j3-2];
x3i=a[j1-1]-a[j3-1];
a[j0-2]=x0r+x2r;
a[j0-1]=x0i-x2i;
a[j1-2]=x0r-x2r;
a[j1-1]=x0i+x2i;
x0r=x1r+x3i;
x0i=x1i+x3r;
a[j2-2]=wk1r*x0r-wk1i*x0i;
a[j2-1]=wk1r*x0i+wk1i*x0r;
x0r=x1r-x3i;
x0i=x1i-x3r;
a[j3-2]=wk3r*x0r+wk3i*x0i;
a[j3-1]=wk3r*x0i-wk3i*x0r;
x0r=a[j0]+a[j2];
x0i=-a[j0+1]-a[j2+1];
x1r=a[j0]-a[j2];
x1i=-a[j0+1]+a[j2+1];
x2r=a[j1]+a[j3];
x2i=a[j1+1]+a[j3+1];
x3r=a[j1]-a[j3];
x3i=a[j1+1]-a[j3+1];
a[j0]=x0r+x2r;
a[j0+1]=x0i-x2i;
a[j1]=x0r-x2r;
a[j1+1]=x0i+x2i;
x0r=x1r+x3i;
x0i=x1i+x3r;
a[j2]=wn4r*(x0r-x0i);
a[j2+1]=wn4r*(x0i+x0r);
x0r=x1r-x3i;
x0i=x1i-x3r;
a[j3]=-wn4r*(x0r+x0i);
a[j3+1]=-wn4r*(x0i-x0r);
x0r=a[j0+2]+a[j2+2];
x0i=-a[j0+3]-a[j2+3];
x1r=a[j0+2]-a[j2+2];
x1i=-a[j0+3]+a[j2+3];
x2r=a[j1+2]+a[j3+2];
x2i=a[j1+3]+a[j3+3];
x3r=a[j1+2]-a[j3+2];
x3i=a[j1+3]-a[j3+3];
a[j0+2]=x0r+x2r;
a[j0+3]=x0i-x2i;
a[j1+2]=x0r-x2r;
a[j1+3]=x0i+x2i;
x0r=x1r+x3i;
x0i=x1i+x3r;
a[j2+2]=wk1i*x0r-wk1r*x0i;
a[j2+3]=wk1i*x0i+wk1r*x0r;
x0r=x1r-x3i;
x0i=x1i-x3r;
a[j3+2]=wk3i*x0r+wk3r*x0i;
a[j3+3]=wk3i*x0i-wk3r*x0r;

}

template <class T> void makect(int nc,int *ip,T *c){
int j,nch;
T delta;
    
ip[1]=nc;
if(nc>1){
  nch=nc>>1;
  delta=atan(1.0)/nch;
  c[0]=cos(delta*nch);
  c[nch]=0.5*c[0];
  for(j=1;j<nch;j++){
     c[j]=0.5*cos(delta*j);
     c[nc-j]=0.5*sin(delta*j);
  }
}

}


template <class T> void bitrv2(int n,int *ip,T *a){

int j,j1,k,k1,l,m,nh,nm;
T xr,xi,yr,yi;

m=1;
for(l=n>>2;l>8;l>>=2){
   m<<=1;
}
nh=n>>1;
nm=4*m;
if(l==8){
  for(k=0;k<m;k++){
     for(j=0;j<k;j++){
        j1=4*j+2*ip[m+k];
        k1=4*k+2*ip[m+j];
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nm;
        k1+=2*nm;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nm;
        k1-=nm;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nm;
        k1+=2*nm;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nh;
        k1+=2;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nm;
        k1-=2*nm;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nm;
        k1+=nm;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nm;
        k1-=2*nm;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=2;
        k1+=nh;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nm;
        k1+=2*nm;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nm;
        k1-=nm;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nm;
        k1+=2*nm;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nh;
        k1-=2;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nm;
        k1-=2*nm;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nm;
        k1+=nm;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nm;
        k1-=2*nm;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
     }
     k1=4*k+2*ip[m+k];
     j1=k1+2;
     k1+=nh;
     xr=a[j1];
     xi=a[j1+1];
     yr=a[k1];
     yi=a[k1+1];
     a[j1]=yr;
     a[j1+1]=yi;
     a[k1]=xr;
     a[k1+1]=xi;
     j1+=nm;
     k1+=2*nm;
     xr=a[j1];
     xi=a[j1+1];
     yr=a[k1];
     yi=a[k1+1];
     a[j1]=yr;
     a[j1+1]=yi;
     a[k1]=xr;
     a[k1+1]=xi;
     j1+=nm;
     k1-=nm;
     xr=a[j1];
     xi=a[j1+1];
     yr=a[k1];
     yi=a[k1+1];
     a[j1]=yr;
     a[j1+1]=yi;
     a[k1]=xr;
     a[k1+1]=xi;
     j1-=2;
     k1-=nh;
     xr=a[j1];
     xi=a[j1+1];
     yr=a[k1];
     yi=a[k1+1];
     a[j1]=yr;
     a[j1+1]=yi;
     a[k1]=xr;
     a[k1+1]=xi;
     j1+=nh+2;
     k1+=nh+2;
     xr=a[j1];
     xi=a[j1+1];
     yr=a[k1];
     yi=a[k1+1];
     a[j1]=yr;
     a[j1+1]=yi;
     a[k1]=xr;
     a[k1+1]=xi;
     j1-=nh-nm;
     k1+=2*nm-2;
     xr=a[j1];
     xi=a[j1+1];
     yr=a[k1];
     yi=a[k1+1];
     a[j1]=yr;
     a[j1+1]=yi;
     a[k1]=xr;
     a[k1+1]=xi;
  }
}else{
  for(k=0;k<m;k++){
     for(j=0;j<k;j++){
        j1=4*j+ip[m+k];
        k1=4*k+ip[m+j];
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nm;
        k1+=nm;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nh;
        k1+=2;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nm;
        k1-=nm;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=2;
        k1+=nh;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nm;
        k1+=nm;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nh;
        k1-=2;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nm;
        k1-=nm;
        xr=a[j1];
        xi=a[j1+1];
        yr=a[k1];
        yi=a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
     }
     k1=4*k+ip[m+k];
     j1=k1+2;
     k1+=nh;
     xr=a[j1];
     xi=a[j1+1];
     yr=a[k1];
     yi=a[k1+1];
     a[j1]=yr;
     a[j1+1]=yi;
     a[k1]=xr;
     a[k1+1]=xi;
     j1+=nm;
     k1+=nm;
     xr=a[j1];
     xi=a[j1+1];
     yr=a[k1];
     yi=a[k1+1];
     a[j1]=yr;
     a[j1+1]=yi;
     a[k1]=xr;
     a[k1+1]=xi;
  }
}

}


template <class T> void bitrv2conj(int n,int *ip,T *a){
int j,j1,k,k1,l,m,nh,nm;
T xr,xi,yr,yi;
    
m=1;
for(l=n>>2;l>8;l>>=2){
   m<<=1;
}
nh=n>>1;
nm=4*m;
if(l==8){
  for(k=0;k<m;k++){
     for(j=0;j<k;j++){
        j1=4*j+2*ip[m+k];
        k1=4*k+2*ip[m+j];
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nm;
        k1+=2*nm;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nm;
        k1-=nm;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nm;
        k1+=2*nm;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nh;
        k1+=2;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nm;
        k1-=2*nm;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nm;
        k1+=nm;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nm;
        k1-=2*nm;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=2;
        k1+=nh;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nm;
        k1+=2*nm;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nm;
        k1-=nm;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nm;
        k1+=2*nm;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nh;
        k1-=2;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nm;
        k1-=2*nm;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nm;
        k1+=nm;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nm;
        k1-=2*nm;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
     }
     k1=4*k+2*ip[m+k];
     j1=k1+2;
     k1+=nh;
     a[j1-1]=-a[j1-1];
     xr=a[j1];
     xi=-a[j1+1];
     yr=a[k1];
     yi=-a[k1+1];
     a[j1]=yr;
     a[j1+1]=yi;
     a[k1]=xr;
     a[k1+1]=xi;
     a[k1+3]=-a[k1+3];
     j1+=nm;
     k1+=2*nm;
     xr=a[j1];
     xi=-a[j1+1];
     yr=a[k1];
     yi=-a[k1+1];
     a[j1]=yr;
     a[j1+1]=yi;
     a[k1]=xr;
     a[k1+1]=xi;
     j1+=nm;
     k1-=nm;
     xr=a[j1];
     xi=-a[j1+1];
     yr=a[k1];
     yi=-a[k1+1];
     a[j1]=yr;
     a[j1+1]=yi;
     a[k1]=xr;
     a[k1+1]=xi;
     j1-=2;
     k1-=nh;
     xr=a[j1];
     xi=-a[j1+1];
     yr=a[k1];
     yi=-a[k1+1];
     a[j1]=yr;
     a[j1+1]=yi;
     a[k1]=xr;
     a[k1+1]=xi;
     j1+=nh+2;
     k1+=nh+2;
     xr=a[j1];
     xi=-a[j1+1];
     yr=a[k1];
     yi=-a[k1+1];
     a[j1]=yr;
     a[j1+1]=yi;
     a[k1]=xr;
     a[k1+1]=xi;
     j1-=nh-nm;
     k1+=2*nm-2;
     a[j1-1]=-a[j1-1];
     xr=a[j1];
     xi=-a[j1+1];
     yr=a[k1];
     yi=-a[k1+1];
     a[j1]=yr;
     a[j1+1]=yi;
     a[k1]=xr;
     a[k1+1]=xi;
     a[k1+3]=-a[k1+3];
  }
}else{
  for(k=0;k<m;k++){
     for(j=0;j<k;j++){
        j1=4*j+ip[m+k];
        k1=4*k+ip[m+j];
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nm;
        k1+=nm;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nh;
        k1+=2;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nm;
        k1-=nm;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=2;
        k1+=nh;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1+=nm;
        k1+=nm;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nh;
        k1-=2;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
        j1-=nm;
        k1-=nm;
        xr=a[j1];
        xi=-a[j1+1];
        yr=a[k1];
        yi=-a[k1+1];
        a[j1]=yr;
        a[j1+1]=yi;
        a[k1]=xr;
        a[k1+1]=xi;
     }
     k1=4*k+ip[m+k];
     j1=k1+2;
     k1+=nh;
     a[j1-1]=-a[j1-1];
     xr=a[j1];
     xi=-a[j1+1];
     yr=a[k1];
     yi=-a[k1+1];
     a[j1]=yr;
     a[j1+1]=yi;
     a[k1]=xr;
     a[k1+1]=xi;
     a[k1+3]=-a[k1+3];
     j1+=nm;
     k1+=nm;
     a[j1-1]=-a[j1-1];
     xr=a[j1];
     xi=-a[j1+1];
     yr=a[k1];
     yi=-a[k1+1];
     a[j1]=yr;
     a[j1+1]=yi;
     a[k1]=xr;
     a[k1+1]=xi;
     a[k1+3]=-a[k1+3];
  }
}

}


template <class T> void bitrv216(T *a){
T x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,
  x5r,x5i,x7r,x7i,x8r,x8i,x10r,x10i,
  x11r,x11i,x12r,x12i,x13r,x13i,x14r,x14i;
    
x1r=a[2];
x1i=a[3];
x2r=a[4];
x2i=a[5];
x3r=a[6];
x3i=a[7];
x4r=a[8];
x4i=a[9];
x5r=a[10];
x5i=a[11];
x7r=a[14];
x7i=a[15];
x8r=a[16];
x8i=a[17];
x10r=a[20];
x10i=a[21];
x11r=a[22];
x11i=a[23];
x12r=a[24];
x12i=a[25];
x13r=a[26];
x13i=a[27];
x14r=a[28];
x14i=a[29];
a[2]=x8r;
a[3]=x8i;
a[4]=x4r;
a[5]=x4i;
a[6]=x12r;
a[7]=x12i;
a[8]=x2r;
a[9]=x2i;
a[10]=x10r;
a[11]=x10i;
a[14]=x14r;
a[15]=x14i;
a[16]=x1r;
a[17]=x1i;
a[20]=x5r;
a[21]=x5i;
a[22]=x13r;
a[23]=x13i;
a[24]=x3r;
a[25]=x3i;
a[26]=x11r;
a[27]=x11i;
a[28]=x7r;
a[29]=x7i;
}


template <class T> void bitrv216neg(T *a){
T x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,
  x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,
  x9r,x9i,x10r,x10i,x11r,x11i,x12r,x12i,
  x13r,x13i,x14r,x14i,x15r,x15i;
    
x1r=a[2];
x1i=a[3];
x2r=a[4];
x2i=a[5];
x3r=a[6];
x3i=a[7];
x4r=a[8];
x4i=a[9];
x5r=a[10];
x5i=a[11];
x6r=a[12];
x6i=a[13];
x7r=a[14];
x7i=a[15];
x8r=a[16];
x8i=a[17];
x9r=a[18];
x9i=a[19];
x10r=a[20];
x10i=a[21];
x11r=a[22];
x11i=a[23];
x12r=a[24];
x12i=a[25];
x13r=a[26];
x13i=a[27];
x14r=a[28];
x14i=a[29];
x15r=a[30];
x15i=a[31];
a[2]=x15r;
a[3]=x15i;
a[4]=x7r;
a[5]=x7i;
a[6]=x11r;
a[7]=x11i;
a[8]=x3r;
a[9]=x3i;
a[10]=x13r;
a[11]=x13i;
a[12]=x5r;
a[13]=x5i;
a[14]=x9r;
a[15]=x9i;
a[16]=x1r;
a[17]=x1i;
a[18]=x14r;
a[19]=x14i;
a[20]=x6r;
a[21]=x6i;
a[22]=x10r;
a[23]=x10i;
a[24]=x2r;
a[25]=x2i;
a[26]=x12r;
a[27]=x12i;
a[28]=x4r;
a[29]=x4i;
a[30]=x8r;
a[31]=x8i;
}


template <class T> void bitrv208(T *a){
T x1r,x1i,x3r,x3i,x4r,x4i,x6r,x6i;
    
x1r=a[2];
x1i=a[3];
x3r=a[6];
x3i=a[7];
x4r=a[8];
x4i=a[9];
x6r=a[12];
x6i=a[13];
a[2]=x4r;
a[3]=x4i;
a[6]=x6r;
a[7]=x6i;
a[8]=x1r;
a[9]=x1i;
a[12]=x3r;
a[13]=x3i;
}


template <class T> void bitrv208neg(T *a){
T x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,
  x5r,x5i,x6r,x6i,x7r,x7i;
    
x1r=a[2];
x1i=a[3];
x2r=a[4];
x2i=a[5];
x3r=a[6];
x3i=a[7];
x4r=a[8];
x4i=a[9];
x5r=a[10];
x5i=a[11];
x6r=a[12];
x6i=a[13];
x7r=a[14];
x7i=a[15];
a[2]=x7r;
a[3]=x7i;
a[4]=x3r;
a[5]=x3i;
a[6]=x5r;
a[7]=x5i;
a[8]=x1r;
a[9]=x1i;
a[10]=x6r;
a[11]=x6i;
a[12]=x2r;
a[13]=x2i;
a[14]=x4r;
a[15]=x4i;
}


template <class T> int cfttree(int n,int j,int k,T *a,int nw,T *w){

int i,isplt,m;
    
if((k & 3)!=0){
  isplt=k & 1;
  if(isplt !=0){
    cftmdl1(n,&a[j-n],&w[nw-(n>>1)]);
  }else{
    cftmdl2(n,&a[j-n],&w[nw-n]);
  }
}else{
  m=n;
  for(i=k;(i & 3)==0;i>>=2){
     m<<=2;
  }
  isplt=i & 1;
  if(isplt!=0){
    while(m>128){
         cftmdl1(m,&a[j-m],&w[nw-(m>>1)]);
         m>>=2;
    }
  }else{
    while(m>128){
         cftmdl2(m,&a[j-m],&w[nw-m]);
         m>>=2;
    }
  }
}
return(isplt);
}

template <class T> void cftleaf(int n,int isplt,T *a,int nw,T *w){
    
if(n==512){
  cftmdl1(128,a,&w[nw-64]);
  cftf161(a,&w[nw-8]);
  cftf162(&a[32],&w[nw-32]);
  cftf161(&a[64],&w[nw-8]);
  cftf161(&a[96],&w[nw-8]);
  cftmdl2(128,&a[128],&w[nw-128]);
  cftf161(&a[128],&w[nw-8]);
  cftf162(&a[160],&w[nw-32]);
  cftf161(&a[192],&w[nw-8]);
  cftf162(&a[224],&w[nw-32]);
  cftmdl1(128,&a[256],&w[nw-64]);
  cftf161(&a[256],&w[nw-8]);
  cftf162(&a[288],&w[nw-32]);
  cftf161(&a[320],&w[nw-8]);
  cftf161(&a[352],&w[nw-8]);
  if(isplt!=0){
    cftmdl1(128,&a[384],&w[nw-64]);
    cftf161(&a[480],&w[nw-8]);
  }else{
    cftmdl2(128,&a[384],&w[nw-128]);
    cftf162(&a[480],&w[nw-32]);
  }
  cftf161(&a[384],&w[nw-8]);
  cftf162(&a[416],&w[nw-32]);
  cftf161(&a[448],&w[nw-8]);
}else{
  cftmdl1(64,a,&w[nw-32]);
  cftf081(a,&w[nw-8]);
  cftf082(&a[16],&w[nw-8]);
  cftf081(&a[32],&w[nw-8]);
  cftf081(&a[48],&w[nw-8]);
  cftmdl2(64,&a[64],&w[nw-64]);
  cftf081(&a[64],&w[nw-8]);
  cftf082(&a[80],&w[nw-8]);
  cftf081(&a[96],&w[nw-8]);
  cftf082(&a[112],&w[nw-8]);
  cftmdl1(64,&a[128],&w[nw-32]);
  cftf081(&a[128],&w[nw-8]);
  cftf082(&a[144],&w[nw-8]);
  cftf081(&a[160],&w[nw-8]);
  cftf081(&a[176],&w[nw-8]);
  if(isplt!=0){
    cftmdl1(64,&a[192],&w[nw-32]);
    cftf081(&a[240],&w[nw-8]);
  }else{
    cftmdl2(64,&a[192],&w[nw-64]);
    cftf082(&a[240],&w[nw-8]);
  }
  cftf081(&a[192],&w[nw-8]);
  cftf082(&a[208],&w[nw-8]);
  cftf081(&a[224],&w[nw-8]);
}
}





template <class T> void cftfx41(int n,T *a,int nw,T *w){
    
if(n==128){
  cftf161(a,&w[nw-8]);
  cftf162(&a[32],&w[nw-32]);
  cftf161(&a[64],&w[nw-8]);
  cftf161(&a[96],&w[nw-8]);
}else{
  cftf081(a,&w[nw-8]);
  cftf082(&a[16],&w[nw-8]);
  cftf081(&a[32],&w[nw-8]);
  cftf081(&a[48],&w[nw-8]);
}
}


template <class T> void cftrec4(int n,T *a,int nw,T *w){

int isplt,j,k,m;

m=n;
while(m>512){
     m>>=2;
     cftmdl1(m,&a[n-m],&w[nw-(m>>1)]);
}
cftleaf(m,1,&a[n-m],nw,w);
k=0;
for(j=n-m;j>0;j-=m){
   k++;
   isplt=cfttree(m,j,k,a,nw,w);
   cftleaf(m,isplt,&a[j-m],nw,w);
}
}



template <class T> void cftbsub(int n,T *a,int *ip,int nw,T *w){
    
if(n>8){
  if(n>32){
    cftb1st(n,a,&w[nw-(n>>2)]);
    if(n>512){
      cftrec4(n,a,nw,w);
    }else if(n>128){
      cftleaf(n,1,a,nw,w);
    }else{
      cftfx41(n,a,nw,w);
    }
    bitrv2conj(n,ip,a);
  }else if(n==32){
    cftf161(a,&w[nw-8]);
    bitrv216neg(a);
  }else{
    cftf081(a,w);
    bitrv208neg(a);
  }
}else if(n==8){
  cftb040(a);
}else if(n==4){
  cftx020(a);
}
}

template <class T> void cftfsub(int n,T *a,int *ip,int nw,T *w){
    
if(n>8){
  if(n>32){
    cftf1st(n,a,&w[nw-(n>>2)]);
    if(n>512){
      cftrec4(n,a,nw,w);
    }else if(n>128){
      cftleaf(n,1,a,nw,w);
    }else{
      cftfx41(n,a,nw,w);
    }
    bitrv2(n,ip,a);
  }else if(n==32){
    cftf161(a,&w[nw-8]);
    bitrv216(a);
  }else{
    cftf081(a,w);
    bitrv208(a);
  }
}else if(n==8){
  cftf040(a);
}else if(n==4){
  cftx020(a);
}

}


template <class T> void rdft(int n,int isgn,T *a,int *ip,T *w){
int nw,nc;
T xi;
    
nw=ip[0];
if(n>(nw<<2)){
  nw=n>>2;
  makewt(nw,ip,w);
}
nc=ip[1];
if(n>(nc<<2)){
  nc=n>>2;
  makect(nc,ip,w+nw);
}
if(isgn>=0){
  if(n>4){
    cftfsub(n,a,ip,nw,w);
    rftfsub(n,a,nc,w+nw);
  }else if(n==4){
    cftfsub(n,a,ip,nw,w);
  }
  xi=a[0]-a[1];
  a[0]+=a[1];
  a[1]=xi;
}else{
  a[1]=0.5*(a[0]-a[1]);
  a[0]-=a[1];
  if(n>4){
    rftbsub(n,a,nc,w+nw);
    cftbsub(n,a,ip,nw,w);
  }else if(n==4){
    cftbsub(n,a,ip,nw,w);
  }
}

}


template <class T> void ddct(int n,int isgn,T *a,int *ip,T *w){
int j,nw,nc;
T xr;
    
nw=ip[0];
if(n>(nw<<2)){
  nw=n>>2;
  makewt(nw,ip,w);
}
nc=ip[1];
if(n>nc){
  nc=n;
  makect(nc,ip,w+nw);
}
if(isgn<0){
  xr=a[n-1];
  for(j=n-2;j>=2;j-=2){
     a[j+1]=a[j]-a[j-1];
     a[j]+=a[j-1];
  }
  a[1]=a[0]-xr;
  a[0]+=xr;
  if(n>4){
    rftbsub(n,a,nc,w+nw);
    cftbsub(n,a,ip,nw,w);
  }else if(n==4){
    cftbsub(n,a,ip,nw,w);
  }
}
dctsub(n,a,nc,w+nw);
if(isgn>=0){
  if(n>4){
    cftfsub(n,a,ip,nw,w);
    rftfsub(n,a,nc,w+nw);
  }else if(n==4){
    cftfsub(n,a,ip,nw,w);
  }
  xr=a[0]-a[1];
  a[0]+=a[1];
  for(j=2;j<n;j+=2){
     a[j-1]=a[j]-a[j+1];
     a[j]+=a[j+1];
  }
  a[n-1]=xr;
}

}


template <class T> void ddst(int n,int isgn,T *a,int *ip,T *w){
int j,nw,nc;
T xr;
    
nw=ip[0];
if(n>(nw<<2)){
  nw=n>>2;
  makewt(nw,ip,w);
}
nc=ip[1];
if(n>nc){
  nc=n;
  makect(nc,ip,w+nw);
}
if(isgn<0){
  xr=a[n-1];
  for(j=n-2;j >=2;j-=2){
     a[j+1]=-a[j]-a[j-1];
     a[j]-=a[j-1];
  }
   a[1]=a[0]+xr;
   a[0]-=xr;
   if(n>4){
     rftbsub(n,a,nc,w+nw);
     cftbsub(n,a,ip,nw,w);
   }else if(n==4){
     cftbsub(n,a,ip,nw,w);
   }
}
dstsub(n,a,nc,w+nw);
if(isgn>=0){
  if(n>4){
    cftfsub(n,a,ip,nw,w);
    rftfsub(n,a,nc,w+nw);
  }else if(n==4){
    cftfsub(n,a,ip,nw,w);
  }
  xr=a[0]-a[1];
  a[0]+=a[1];
  for(j=2;j<n;j+=2){
     a[j-1]=-a[j]-a[j+1];
     a[j]-=a[j+1];
  }
  a[n-1]=-xr;
}

}


template <class T> void dfct(int n,T *a,T *t,int *ip,T *w){
int j,k,l,m,mh,nw,nc;
T xr,xi,yr,yi;
    
nw=ip[0];
if(n>(nw<<3)){
  nw=n>>3;
  makewt(nw,ip,w);
}
nc=ip[1];
if(n>(nc<<1)){
  nc=n>>1;
  makect(nc,ip,w+nw);
}
m=n>>1;
yi=a[m];
xi=a[0]+a[n];
a[0]-=a[n];
t[0]=xi-yi;
t[m]=xi+yi;
if(n>2){
  mh=m>>1;
  for(j=1;j<mh;j++){
     k=m-j;
     xr=a[j]-a[n-j];
     xi=a[j]+a[n-j];
     yr=a[k]-a[n-k];
     yi=a[k]+a[n-k];
     a[j]=xr;
     a[k]=yr;
     t[j]=xi-yi;
     t[k]=xi+yi;
  }
  t[mh]=a[mh]+a[n-mh];
  a[mh]-=a[n-mh];
  dctsub(m,a,nc,w+nw);
  if(m>4){
    cftfsub(m,a,ip,nw,w);
    rftfsub(m,a,nc,w+nw);
  }else if(m==4){
    cftfsub(m,a,ip,nw,w);
  }
  a[n-1]=a[0]-a[1];
  a[1]=a[0]+a[1];
  for(j=m-2;j >=2;j-=2){
     a[2*j+1]=a[j]+a[j+1];
     a[2*j-1]=a[j]-a[j+1];
  }
  l=2;
  m=mh;
  while(m>=2){
       dctsub(m,t,nc,w+nw);
       if(m >4){
         cftfsub(m,t,ip,nw,w);
         rftfsub(m,t,nc,w+nw);
       }else if(m==4){
         cftfsub(m,t,ip,nw,w);
       }
       a[n-l]=t[0]-t[1];
       a[l]=t[0]+t[1];
       k=0;
       for(j=2;j<m;j+=2){
          k+=l<<2;
          a[k-l]=t[j]-t[j+1];
          a[k+l]=t[j]+t[j+1];
       }
       l<<=1;
       mh=m>>1;
       for(j=0;j<mh;j++){
          k=m-j;
          t[j]=t[m+k]-t[m+j];
          t[k]=t[m+k]+t[m+j];
       }
       t[mh]=t[m+mh];
       m=mh;
  }
  a[l]=t[0];
  a[n]=t[2]-t[1];
  a[0]=t[2]+t[1];
}else{
  a[1]=a[0];
  a[2]=t[0];
  a[0]=t[1];
}

}


template <class T> void dfst(int n,T *a,T *t,int *ip,T *w){
int j,k,l,m,mh,nw,nc;
T xr,xi,yr,yi;
    
nw=ip[0];
if(n>(nw<<3)){
  nw=n>>3;
  makewt(nw,ip,w);
}
nc=ip[1];
if(n>(nc<<1)){
  nc=n>>1;
  makect(nc,ip,w+nw);
}
if(n>2){
  m=n>>1;
  mh=m>>1;
  for(j=1;j<mh;j++){
     k=m-j;
     xr=a[j]+a[n-j];
     xi=a[j]-a[n-j];
     yr=a[k]+a[n-k];
     yi=a[k]-a[n-k];
     a[j]=xr;
     a[k]=yr;
     t[j]=xi+yi;
     t[k]=xi-yi;
  }
  t[0]=a[mh]-a[n-mh];
  a[mh]+=a[n-mh];
  a[0]=a[m];
  dstsub(m,a,nc,w+nw);
  if(m>4){
    cftfsub(m,a,ip,nw,w);
    rftfsub(m,a,nc,w+nw);
  }else if(m==4){
    cftfsub(m,a,ip,nw,w);
  }
  a[n-1]=a[1]-a[0];
  a[1]=a[0]+a[1];
  for(j=m-2;j>=2;j-=2){
     a[2*j+1]=a[j]-a[j+1];
     a[2*j-1]=-a[j]-a[j+1];
  }
  l=2;
  m=mh;
  while(m>=2){
       dstsub(m,t,nc,w+nw);
       if(m>4){
         cftfsub(m,t,ip,nw,w);
         rftfsub(m,t,nc,w+nw);
       }else if(m==4){
         cftfsub(m,t,ip,nw,w);
       }
        a[n-l]=t[1]-t[0];
        a[l]=t[0]+t[1];
        k=0;
        for(j=2;j<m;j+=2){
           k+=l<<2;
           a[k-l]=-t[j]-t[j+1];
           a[k+l]=t[j]-t[j+1];
        }
        l<<=1;
        mh=m>>1;
        for(j=1;j<mh;j++){
           k=m-j;
           t[j]=t[m+k]+t[m+j];
           t[k]=t[m+k]-t[m+j];
         }
         t[0]=t[m+mh];
         m=mh;
  }
  a[l]=t[0];
}
a[0]=0;

}

template <class T> void cdft(int n,int isgn,T *a,int *ip,T *w){

int nw;
    
nw=ip[0];
if(n>(nw<<2)){
  nw=n>>2;
  makewt(nw,ip,w);
}
if(isgn>=0){
  cftfsub(n,a,ip,nw,w);
}else{
  cftbsub(n,a,ip,nw,w);
}

}




#endif