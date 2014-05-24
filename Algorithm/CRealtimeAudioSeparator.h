

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



#include "CAudioModalizer.h"

#ifndef CREALTIMEAUDIOSEPARATOR_HH
#define CREALTIMEAUDIOSEPARATOR_HH



template <class T> inline T nmf_error(CRingBufferMatrix<T>& original,CRingBufferMatrix<T>& approx,T eps=1.0e-12){
const int col=original._col;
const int row=original._row;
T err=0.0;

for(int j=row-1;j<row;j++){
   for(int i=0;i<col;i++){
      if(approx[i][j]<eps) continue;
	  if(original[j][i]<eps) continue;
	  err+=square(original[j][i]-approx[i][j]);
   }
}

return(err);
}



template <class T> class CRealtimeAudioSeparator{

private:

CDenseMatrix<T> _dictionary; // teacher spectrum data
CRingBufferMatrix<T> _activation;
CRingBufferMatrix<T> _original;
CRingBufferMatrix<T> _approx;
CVector<T> _coef;
CDenseMatrix<T> _weight;
CVector<T> _buffers;
CVector<T> _damping;
T _threthold;
T _decay;
unsigned int _max_iter;
unsigned int _halfWindowSize;
unsigned int _num_frames;
unsigned int _num_timbres;


public:

CVector<int> _active;

CRealtimeAudioSeparator(){}

CRealtimeAudioSeparator(unsigned int max_iter,unsigned int windowSize,unsigned int num_timbres,unsigned int num_frames,T threthold){
_init(max_iter,windowSize,num_timbres,num_timbres,threthold);
}


void _init(unsigned int max_iter,unsigned int windowSize,unsigned int num_timbres,unsigned int num_frames,T threthold){
_dictionary._resize(num_timbres,windowSize/2);
_activation._resize(num_frames,num_timbres);
_activation._zero();
_original._resize(num_frames,windowSize/2);
_original._zero();
_approx._resize(windowSize/2,num_frames);
_approx._zero();
_coef._resize(windowSize/2);
_active._resize(num_timbres);
_active=0;
_threthold=threthold;
_halfWindowSize=windowSize/2;
_num_frames=num_frames;
_num_timbres=num_timbres;
_max_iter=max_iter;
_weight._resize(_num_timbres,windowSize/2);
_weight=0;
_decay=0.999;
_damping._resize(_num_timbres);
}

void _registerTeacher(int timbre,T* teacher,T damp){
_dictionary._copy_vec(timbre,teacher);
_active[timbre]=1;
_damping[timbre]=damp;
}

void _unregisterTeacher(int timbre){
_dictionary._zeroRow(timbre);
_active[timbre]=0;
}


void _calcWeight(){

const int HalfWindowSize=_halfWindowSize;

_buffers._resize(_num_timbres);
_buffers=1;

T maxt=0;
for(int t=0;t<_num_timbres;++t){
   if(_active[t]){
     if(_damping[t]>maxt) maxt=_damping[t];
   }
}
for(int t=0;t<_num_timbres;++t){
   if(_active[t]){
     _damping[t]=0.999*_damping[t]/maxt;
   }
}

for(int w=0;w<HalfWindowSize;++w){
   T sum=0;
   for(int t=0;t<_num_timbres;++t){
       if(_active[t]) sum+=_dictionary[t][w];
   }
   
   if(sum>0){
     for(int t=0;t<_num_timbres;++t){
        if(_active[t]) _weight[t][w]=_dictionary[t][w]/sum;
     }
   }
}

}


void _get_current_result(const unsigned int timbre,CVector2d<T>* source,CVector2d<T>* result){

const unsigned int WindowSize=_halfWindowSize*2;
const unsigned int HalfWindowSize=_halfWindowSize;

if(_active[timbre]==0){
  for(int i=0;i<WindowSize;i++) result[i]=0;
  return;
}

T rate;

result[0]=0;
if(_approx[HalfWindowSize-1][_num_frames-1]<1.0e-12){
  result[HalfWindowSize]=0.0;
}else{
  rate=sqrt(_dictionary[timbre][HalfWindowSize-1]*_activation[_num_frames-1][timbre]/_approx[HalfWindowSize-1][_num_frames-1]);
  result[HalfWindowSize]=rate*source[HalfWindowSize];
}

for(int i=1;i<HalfWindowSize;i++){
   if(_approx[i][_num_frames-1]<1.0e-12){
     result[i]=0.0;
	 result[WindowSize-i]=0.0;
	 continue;
   }
   rate=sqrt(_dictionary[timbre][i-1]*_activation[_num_frames-1][timbre]/_approx[i-1][_num_frames-1]);
   result[i]=rate*source[i];
   result[i]._conjugate(result[WindowSize-i]);
}

}



void _process(T* source){

const unsigned int max_iter=_max_iter;
const unsigned int HalfWindowSize=_halfWindowSize;
const unsigned int num_frames=_num_frames;
const unsigned int num_timbres=_num_timbres;
const T decay=_decay;
const T thr=1/(T)num_timbres;

_original._push(source);
_activation._advance();

T maxv=0;
T power=0;
for(int i=0;i<HalfWindowSize;++i){
   power+=source[i];
   if(source[i]>maxv) maxv=source[i];
}

if(power<=0){
  for(int t=0;t<num_timbres;++t){
     _activation[num_frames-1][t]=0;
  }
  for(int i=0;i<HalfWindowSize;++i){
     _approx[i][num_frames-1]=0;
  }
  return;
}


for(int t=0;t<num_timbres;++t){
   if(_active[t]){
     
     if(_activation[num_frames-2][t]>0.01) _activation[num_frames-1][t]=_activation[num_frames-2][t];
     else _activation[num_frames-1][t]=0.01;
     
     T sum=0;
     for(int w=0;w<HalfWindowSize;++w){
        sum+=source[w]*_weight[t][w];
     }
     
     if(sum>power*thr*_buffers[t]) _activation[num_frames-1][t]=sum;
     else _activation[num_frames-1][t]*=decay;
     
     if(sum>power*thr*_buffers[t]){
       if(sum>power*thr*2.7) sum=power*thr*2.7;
       _buffers[t]=sum/(power*thr);
     }else{
       _buffers[t]*=_damping[t];
       if(_buffers[t]<1.0) _buffers[t]=1.0;
     }
     
     for(int w=0;w<HalfWindowSize;++w){
        if(_dictionary[t][w]>0.6){
          if(source[w]/maxv<0.1){
            _activation[num_frames-1][t]*=decay;
            break;
          }
        }
     }
     
   }
}



T total=0;
for(int j=0;j<num_timbres;++j){
   for(int i=0;i<HalfWindowSize;++i){
      total+=_dictionary[j][i]*_activation[num_frames-1][j];
   }
}
if(total>power){
  for(int j=0;j<num_timbres;++j){
     _activation[num_frames-1][j]*=power/total;
  }
}

for(int i=0;i<HalfWindowSize;++i){
   T sum=0.0;
   for(int j=0;j<num_timbres;++j){
      if(_active[j]) sum+=_dictionary[j][i]*_activation[num_frames-1][j];
   }
   _approx[i][num_frames-1]=sum;
}


T beta=0.9;
const T preset[10]={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
int pid=9;
T prev=FLT_MAX;

for(int iter=0;iter<max_iter;++iter){

   for(int w=0;w<HalfWindowSize;++w){
      if(_approx[w][num_frames-1]>1.0e-12){
	    _coef[w]=pow(_approx[w][num_frames-1],beta-2.0)*_original[num_frames-1][w];
	  }else{ 
	    _coef[w]=0.0;
	  }
   }
   

   for(int t=0;t<num_timbres;++t){
      if(_activation[num_frames-1][t]<=0) continue;
      T a=0.0;
      T b=0.0;
	  if(_active[t]){
	    for(int w=0;w<HalfWindowSize;++w){
		   if(_approx[w][num_frames-1]>1.0e-12){
             if(_weight[t][w]>=thr){
               a+=_dictionary[t][w]*_coef[w];
	           b+=_dictionary[t][w]*pow(_approx[w][num_frames-1],beta-1.0);
             }
		   }
	    }
	    if(b>0){
          _activation[num_frames-1][t]*=a/b;
        }else{
          _activation[num_frames-1][t]=0;
        }
	  }else{
	    _activation[num_frames-1][t]=0;
	  }
   }
      
   for(int t=0;t<num_timbres;++t){
      if(_activation[num_frames-1][t]<=0) continue;
      T sum=0.0;
      for(int i=0;i<HalfWindowSize;++i){
         sum+=_dictionary[t][i]*_activation[num_frames-1][t];
      }
      if(sum>power){
        _activation[num_frames-1][t]*=power/sum;
      }
   }
   
   for(int w=0;w<HalfWindowSize;++w){
      T sum=0.0;
      for(int t=0;t<num_timbres;++t){
         if(_active[t]){
           sum+=_dictionary[t][w]*_activation[num_frames-1][t];
         }
      }
      _approx[w][num_frames-1]=sum;
   }
   
   T err=nmf_error(_original,_approx);
   if(err<_threthold) break;
   if(err>prev){
     pid=(pid+3)%10;
     beta=preset[pid];
   }
   prev=err;
}

}



//first :: window size
//second :: num timbres
//third :: threthold
//4th :: sequence data
/*       format
num_timbres(int) num_events(int) max_time(sec)
timbre_idx0(int) time0(sec)
           .
           .
           .
timbre_idxN(int) timeN(sec)
*/
//5th::output file path
//6~6+num_timbres :: individual timbre wavefile path
void _test(std::vector<std::string>& arg){

int windowsize=atoi(arg[0].c_str());
int num_timbres=atoi(arg[1].c_str());
T max_time;
int step=windowsize/2;
char outputpath[512];

std::ifstream seq(arg[3].c_str());


int timb;
int start;
T time;
int num_events;
seq>>num_events;
seq>>max_time;

CWaveData<T>* waves=new CWaveData<T>[num_timbres];
int length=(int)(waves[0]._sampleRate*max_time+0.5);
const int num_frames=length/step-windowsize/step-1;

CVector<T> source(num_timbres*length);
source=0;
CVector<T> mixture(length);
mixture=0;

for(int i=0;i<num_timbres;i++){
   waves[i]._open(arg[5+i].c_str());
   waves[i]._cache_from_file();
   waves[i]._stereo_to_mono();
}

for(int j=0;j<num_events;j++){
   seq>>timb;
   seq>>time;
   start=(int)(waves[timb]._sampleRate*time+0.5);
   for(int i=0;i<waves[timb]._length;i++){
      source[length*timb+start+i]+=waves[timb]._data[i];
      mixture[start+i]+=waves[timb]._data[i];
   }
}
seq.close();

CWaveData<T> savewav;
for(int i=0;i<num_timbres;i++){
   savewav._set_wave(length,&source[length*i]);
   sprintf(outputpath,"%s/source%d.wav",arg[4].c_str(),i);
   savewav._save(outputpath,0,length);
}

savewav._set_wave(length,mixture._getPointer());
sprintf(outputpath,"%s/mix.wav",arg[4].c_str());
savewav._save(outputpath,0,length);


_init(10,windowsize,num_timbres,8,atof(arg[2].c_str()));


CVector<T> power(windowsize/2);
CFFT_Analyzer<T> fft;
fft._init(windowsize);
fft._setWindow(FFT_WINDOW_TYPE_SINBELL);

for(int i=0;i<num_timbres;i++){
   unsigned int num_sp=calc_average_power_spectrum(fft,waves[i]._data,power._getPointer(),0,waves[i]._length,step);
   
   sprintf(outputpath,"%s/spectrum%d.res",arg[4].c_str(),i);
   std::ofstream ofs(outputpath);

   for(int j=0;j<windowsize/2;j++){
      power[j]*=power[j];
      ofs<<power[j]<<std::endl;
   }

   _registerTeacher(i,power._getPointer(),num_sp);
}
delete [] waves;

_calcWeight();

CVector<T> result(num_timbres*length);
result=0;
CVector<T> fragment(windowsize);
CVector< CVector2d<T> > spectrum(windowsize);
CVector< CVector2d<T> > res_spectrum(windowsize);


for(int f=0;f<num_frames;f++){

   fft._process(&mixture[step*f],&(spectrum[0][0]));
   for(int j=0;j<windowsize/2;j++){
      power[j]=spectrum[j+1]._normSquare();
   }
   _process(power._getPointer());
   
   for(int t=0;t<num_timbres;t++){
      _get_current_result(t,spectrum._getPointer(),res_spectrum._getPointer());
     
      fft._inverse(&res_spectrum[0][0],fragment._getPointer());

      for(int j=0;j<windowsize;j++){
         result[length*t+step*f+j]+=fragment[j];
      }
      
   }
}


double err_av,err;
for(int t=0;t<num_timbres;t++){
   err_av=0.0;
   sprintf(outputpath,"%s/output%d.res",arg[4].c_str(),t);
   std::ofstream ofs(outputpath);
   for(int l=0;l<length;l++){
      err=(result[length*t+l]-source[length*t+l])*(result[length*t+l]-source[length*t+l]);
      err_av+=err;
      ofs<<err<<std::endl;
   }
   err_av/=length;
   std::cout<<"separation error - timbre:"<<t<<" "<<"error:"<<err_av<<std::endl;
   savewav._set_wave(length,&result[length*t]);
   sprintf(outputpath,"%s/result%d.wav",arg[4].c_str(),t);
   savewav._save(outputpath,0,length);
}


}



};



#endif
