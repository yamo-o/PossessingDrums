
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


#include "./Algorithm/PossessingDrums.h"



#define BUFFER_SIZE 256 // audio buffer size


//example for testing sound source separation
void realtime_sound_source_separation_test(){

CRealtimeAudioSeparator<float>* sep=new CRealtimeAudioSeparator<float>();
std::vector<std::string> testarg;
testarg.push_back(std::string("256")); // window size
testarg.push_back(std::string("5")); // quantity of timbre
testarg.push_back(std::string("0.01")); // threshold
testarg.push_back(std::string("./testdata/seq.sq")); // sequence file
testarg.push_back(std::string("./testdata")); // directory for save
testarg.push_back(std::string("./testdata/pianoA4.wav")); // timbre teacher wav file
testarg.push_back(std::string("./testdata/pianoF2.wav"));
testarg.push_back(std::string("./testdata/hihat.wav"));
testarg.push_back(std::string("./testdata/tom.wav"));
testarg.push_back(std::string("./testdata/guiterB3.wav"));
sep->_test(testarg);
delete sep;

}



// offline example
int main(int argc,char *argv[]){


PossessingDrums<float> pd;

std::cout<<"initialize"<<std::endl;

//Initialize.
//1:quantity of timbres, 2:fft window size, 3:audio buffer size, 4:sampling rate.
pd._init(5,BUFFER_SIZE,BUFFER_SIZE,44100); //5 polys


std::cout<<"assigning timbres"<<std::endl;

// Assigning a timbres to a input signal
//1:timbre id, 2:driving signal(input), 3:timbre to assign(output)
pd._setPair(0,"./testdata/pianoA4.wav","./testdata/pianoE3.wav");
pd._setPair(1,"./testdata/pianoF2.wav","./testdata/pianoE3.wav");
pd._setPair(2,"./testdata/hihat.wav","./testdata/pianoE3.wav");
pd._setPair(3,"./testdata/tom.wav","./testdata/pianoE3.wav");
pd._setPair(4,"./testdata/guiterB3.wav","./testdata/pianoE3.wav");


//pre-processing for realtime sound source separation
pd._ready();


CWaveData<float> wav;
wav._open("./testdata/mix.wav"); // input signal
wav._cache_from_file(); // load onto RAM
wav._stereo_to_mono(); // convert to monaural signal
int num_frames=wav._length/BUFFER_SIZE+1;

CWaveData<float> result;

float frame[BUFFER_SIZE];
float resultBuf[BUFFER_SIZE];

std::cout<<"start"<<std::endl;

// play
for(int f=0;f<num_frames;++f){

   if(f!=num_frames-1) memcpy(frame,&wav[BUFFER_SIZE*f],BUFFER_SIZE*sizeof(float));
   else{
     memset(frame,0,BUFFER_SIZE*sizeof(float));
     int j=0;
     for(int i=f*BUFFER_SIZE;i<wav._length;++i,++j){
        frame[j]=wav[i];
     }
   }
   
   memset(resultBuf,0,BUFFER_SIZE*sizeof(float));
   
   pd._record(BUFFER_SIZE,frame); // input
   pd._render(BUFFER_SIZE,resultBuf); // output
   
   result._record(BUFFER_SIZE,resultBuf); // recording
   
}


std::cout<<"finished"<<std::endl;


result._render_Recdata(); // cleanup after recording
result._save("./testdata/result.wav",0,result._length); // save final result


return(0);
}




