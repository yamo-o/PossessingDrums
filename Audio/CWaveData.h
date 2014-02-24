

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

#ifndef CWAVDATA_HH
#define CWAVDATA_HH

#include "CTypes.h"


typedef enum CWAVE_FORMAT{
CWAVE_FORMAT_SINT8,
CWAVE_FORMAT_SINT16,
CWAVE_FORMAT_SINT24,
CWAVE_FORMAT_SINT32,
CWAVE_FORMAT_FLOAT32,
CWAVE_FORMAT_FLOAT64
} CWAVE_FORMAT;


typedef struct RECORD_FRAME{
float* _buf;
RECORD_FRAME *_prev,*_next;
} RECORD_FRAME;


typedef struct CWAVHDR{
char _riff[4];
int _file_size;
char _wave[4];
char _fmt[4];
int _chunk_size;
short int _format_tag;
short int _num_chans;
int _sample_rate;
int _bytes_per_sec;
short int _bytes_per_samp;
short int _bits_per_samp;
char _data[4];
int _data_length;
} CWAVHDR;

static inline void swap32(unsigned int* dwData){
*dwData=((*dwData >> 24) & 0x000000FF) | 
	   ((*dwData >> 8)  & 0x0000FF00) | 
	   ((*dwData << 8)  & 0x00FF0000) | 
	   ((*dwData << 24) & 0xFF000000);
}

static inline void swap64(unsigned long int* dwData){
*dwData=((*dwData >> 48) & 0x00000000000000FF) |
	   ((*dwData >> 32) & 0x000000000000FF00) | 
	   ((*dwData >> 24) & 0x0000000000FF0000) | 
	   ((*dwData >> 8)  & 0x00000000FF000000) | 
	   ((*dwData << 8)  & 0x000000FF00000000) | 
	   ((*dwData << 24) & 0x0000FF0000000000) | 
	   ((*dwData << 32) & 0x00FF000000000000) | 
	   ((*dwData << 48) & 0xFF00000000000000);
}    

static inline void swap16(unsigned short* wData){
*wData=((*wData >> 8) & 0x00FF) | 
	  ((*wData << 8) & 0xFF00);
}

static inline void swap16Buffer(unsigned short *pData,unsigned int dwNumWords){
unsigned long i;
for(i=0;i<dwNumWords;i++){
   swap16(&pData[i]);
}
}


class CWaveData{

public:

FILE* _fd;
bool _opened;
bool _big_endian;
unsigned int _dataOffset;
CWAVE_FORMAT _dataType;
unsigned int _format_tag;
unsigned int _current_pos;
bool _on_memory;
unsigned int _channels;
unsigned int _length;
unsigned int _sampleRate;
int _beginpoint,_endpoint;
float* _data;
float* _cache;
float _cache_size;
bool _loop;
unsigned int _rec_frameSize;
RECORD_FRAME *_rec_frame,*_tail_frame;
int _num_frames;


CWaveData(){
_on_memory=false;
_opened=false;
_fd=0;
_big_endian=false;
_tail_frame=_rec_frame=NULL;
_rec_frameSize=44100*8;
_num_frames=0;
_channels=1;
_sampleRate=44100;
_length=0;
_data=NULL;
_cache=NULL;
_cache_size=0;
_loop=false;
_current_pos=0;
_beginpoint=_endpoint=0;
}

~CWaveData(){
if(_fd>0) _close();
if(_on_memory) free(_data);
if(_num_frames>0){
  RECORD_FRAME* tmp=_rec_frame;
  while(tmp){
     RECORD_FRAME* tmpp=tmp->_next;
	 free(tmp->_buf);
	 delete tmp;
	 tmp=tmpp;
  }
  _rec_frame=_tail_frame=NULL;
  _num_frames=0;
}
}


inline float& operator[](int i){ return(_data[i]); }
inline const float& operator[](int i) const { return(_data[i]); }


void _open(const char* path){
if(_fd>0) return;
_fd=fopen(path,"rb");

if(_fd<=0) return;

char header[12];
if(fread(&header,4,3,_fd)!=3) return;
if(!strncmp(header,"RIFF",4) && !strncmp(&header[8],"WAVE",4)) _getWavInfo();
  
}

void _close(){
if(_fd>0) fclose(_fd);
_fd=0;
_opened=false;
}

void _resize(unsigned int size){
float* tmp=(float*)malloc(_channels*size*sizeof(float));
if(_data){
  if(size<_length) _length=size;
  memcpy(tmp,_data,_channels*_length*sizeof(float));
  free(_data);
}
_data=tmp;
_length=size;
}


void _getWavInfo(){
char fid[4];
int chunkSize;
if(fread(&fid,4,1,_fd)!=1) return;

while(strncmp(fid,"fmt ",4)){
     if(fread(&chunkSize,4,1,_fd)!=1) return;
	 if(_big_endian){
	   swap32((unsigned int*)&chunkSize);
	 }
	 if(fseek(_fd,chunkSize,SEEK_CUR)==-1) return;
     if(fread(&fid,4,1,_fd)!=1) return;
}


if(fread(&chunkSize,4,1,_fd)!=1) return;
if(fread(&_format_tag,2,1,_fd)!=1) return;

if(_big_endian){
  swap16((unsigned short*)&_format_tag);
  swap32((unsigned int*)&chunkSize);
}
if(_format_tag==0xFFFE){
  _dataOffset=ftell(_fd);
  if(fseek(_fd,14,SEEK_CUR)==-1) return;
  unsigned short extSize;
  if(fread(&extSize,2,1,_fd)!=1) return;

  if(_big_endian){
    swap16(&extSize);
  }
  if(extSize==0) return;

  if(fseek(_fd,6,SEEK_CUR)==-1) return;

  if(fread(&_format_tag,2,1,_fd)!=1) return;

  if(_big_endian){
    swap16((unsigned short*)&_format_tag);
  }
  if(fseek(_fd,_dataOffset,SEEK_SET)==-1) return;
}
if(_format_tag!=1 && _format_tag!=3){
  //return;
}


short int temp;
if(fread(&temp,2,1,_fd)!=1) return;

if(_big_endian){
  swap16((unsigned short*)&temp);
}
_channels=(unsigned int)temp;

if(fread(&_sampleRate,4,1,_fd)!=1) return;

if(_big_endian){
  swap32(&_sampleRate);
}

_dataType=(CWAVE_FORMAT)0;
if(fseek(_fd,6,SEEK_CUR)==-1) return;

if(fread(&temp,2,1,_fd)!=1) return;

if(_big_endian){
  swap16((unsigned short*)&temp);
}
if(_format_tag==1){
  if(temp==8)
	_dataType=CWAVE_FORMAT_SINT8;
  else if(temp==16)
	_dataType=CWAVE_FORMAT_SINT16;
  else if(temp==32)
	_dataType=CWAVE_FORMAT_SINT32;
}else if(_format_tag==3){
  if(temp==3)
	_dataType=CWAVE_FORMAT_FLOAT32;
  else if(temp==64)
    _dataType=CWAVE_FORMAT_FLOAT64;
}
if(_dataType==0){
  _dataType=CWAVE_FORMAT_SINT16;
}

if(fseek(_fd,chunkSize-16,SEEK_CUR)==-1) return;
if(fread(&fid,4,1,_fd)!=1) return;

while(strncmp(fid,"data",4)){
     if(fread(&chunkSize,4,1,_fd)!=1) return;
     if(_big_endian){
       swap32((unsigned int*)&chunkSize);
     }
     chunkSize+=chunkSize%2;
     if(fseek(_fd,chunkSize,SEEK_CUR)==-1) return;
     if(fread(&fid,4,1,_fd)!=1) return;
}

int bytes;
if(fread(&bytes,4,1,_fd)!=1) return;
if(_big_endian){
  swap32((unsigned int*)&bytes);
}
_length=8*bytes/temp/_channels;

_dataOffset=ftell(_fd);

_opened=true;
_current_pos=0;
}

void _rewind(){
if(_opened){
  fseek(_fd,_dataOffset,SEEK_SET);
}
_current_pos=0;
}

void _seek(unsigned int pos){
if(pos>_length-1) pos=_length-1;
if(_on_memory){
  _current_pos=pos;
}else if(_opened){
  int offset=pos*_channels;
  if(_dataType==CWAVE_FORMAT_SINT16) fseek(_fd,_dataOffset+2*offset,SEEK_SET);
  else if(_dataType==CWAVE_FORMAT_SINT32) fseek(_fd,_dataOffset+4*offset,SEEK_SET); 
  else if(_dataType==CWAVE_FORMAT_FLOAT32) fseek(_fd,_dataOffset+4*offset,SEEK_SET);
  else if(_dataType==CWAVE_FORMAT_FLOAT64) fseek(_fd,_dataOffset+8*offset,SEEK_SET);
  else if(_dataType==CWAVE_FORMAT_SINT8) fseek(_fd,_dataOffset+offset,SEEK_SET);
  else if(_dataType==CWAVE_FORMAT_SINT24) fseek(_fd,_dataOffset+3*offset,SEEK_SET);
  _current_pos=pos;
}
}

int _read(unsigned int size,float* buf){

int readSize;

if(_current_pos+size>=_length) readSize=_length-_current_pos;
int nSamples=readSize*_channels;
float max_val;

if(_on_memory){
  memcpy(buf,&_data[_channels*_current_pos],nSamples*sizeof(float));
  _current_pos+=readSize;
  if(_loop){
    if(readSize<size){
      _current_pos=0;
	  readSize+=_read(size-readSize,&buf[_channels*readSize]);
    }
  }
  return(readSize);
}

if(_dataType==CWAVE_FORMAT_SINT16){
  short int* ptr=(short int*)buf;
  fread(ptr,2*nSamples,1,_fd);
  if(_big_endian){
	for(int i=nSamples-1;i>=0;i--) swap16((unsigned short*)ptr++);
	ptr=(short int*)buf;
  }
  max_val=1.0f/32768.0f;
  for(int i=nSamples-1;i>=0;i--){
     ptr[i]*=max_val;
     buf[i]=ptr[i];
  }
}else if(_dataType==CWAVE_FORMAT_SINT32){
  int* ptr=(int*)buf;
  fread(ptr,4*nSamples,1,_fd);
  if(_big_endian){
	for(int i=nSamples-1;i>=0;i--) swap32((unsigned int*)ptr++);
	ptr=(int*)buf;
  }
  max_val=1.0f/2147483648.0f;
  for(int i=nSamples-1;i>=0;i--){ 
     ptr[i]*=max_val;
	 buf[i]=ptr[i];
  }
}else if(_dataType==CWAVE_FORMAT_FLOAT32){
  fread(buf,4*nSamples,1,_fd);
  if(_big_endian){
    float* ptr=buf;
	for(int i=nSamples-1;i>=0;i--) swap32((unsigned int*)ptr++);
  }
}else if(_dataType==CWAVE_FORMAT_FLOAT64){
  double *ptr=(double*)buf;
  fread(ptr,8*nSamples,1,_fd);
  if(_big_endian){
    for(int i=nSamples-1;i>=0;i--) swap64((unsigned long int*)ptr++);
  }
}else if(_dataType==CWAVE_FORMAT_SINT8){
  unsigned char* ptr=(unsigned char*)buf;
  fread(ptr,nSamples,1,_fd);
  max_val=1.0f/128.0f;
  for(int i=nSamples-1;i>=0;i--){ 
     buf[i]=(ptr[i]-128)*max_val;
  }
}else if(_dataType==CWAVE_FORMAT_SINT24){
  unsigned int tmp;
  max_val=1.0f/8388608.0f;
  for(int i=0;i<nSamples;i++){
	 fread(&tmp,3,1,_fd);
	 tmp>>=8;
	 if(_big_endian) swap32((unsigned int*)&tmp);
	 buf[i]=tmp*max_val;
  }
}

_current_pos+=readSize;
if(_loop){
  if(readSize<size){
	_rewind();
	readSize+=_read(size-readSize,&buf[_channels*readSize]);
  }
}
return(readSize);
}


void _cache_from_file(){
if(!_opened) return;
if(_data) free(_data);

_data=(float*)malloc(_channels*_length*sizeof(float));
memset(_data,0,_channels*_length*sizeof(float));
_rewind();
short int* ptr=(short int*)_data;
//_read(_length,_data);
fread(ptr,_channels*_length*sizeof(short int),1,_fd);

float max_val=1.0f/32768.0f;
for(int i=(_channels*_length-1);i>=0;i--){
   _data[i]=(float)(ptr[i]);
   _data[i]*=max_val;
}
_endpoint=_length;
_close();

_on_memory=true;

}


void _stereo_to_mono(){
if(_channels==1) return;
if(_on_memory){
  for(int i=0;i<_length;i++){
     _data[i]=_data[2*i];
  }
}
_channels=1;
}

void _set_wave(int length,float* data){
if(_data) free(_data);
_data=(float*)malloc(length*sizeof(float));
memcpy(_data,data,length*sizeof(float));
_channels=1;
_length=length;
_on_memory=true;
}

void _save(const char *fileName,int start,int length){
char name[8192];
strncpy(name,fileName,8192);
if(strstr(name,".wav")==NULL) strcat(name,".wav");
FILE* fd=fopen(name,"wb");

CWAVHDR hdr={"RIF",44,"WAV","fmt",16,1,1,_sampleRate,0,2,16,"dat",0};
hdr._riff[3]='F';
hdr._wave[3]='E';
hdr._fmt[3]=' ';
hdr._data[3]='a';
hdr._num_chans=_channels;
hdr._bits_per_samp=16;
hdr._bytes_per_samp=_channels* hdr._bits_per_samp/8;
hdr._bytes_per_sec=hdr._sample_rate*hdr._bytes_per_samp;

fwrite(&hdr,4,11,fd);

short int data;
for(int i=0;i<length;i++){
   data=_data[start+i]*32768.0f;
   fwrite(&data,2,1,fd);
}

int bytes_per_sample=2;
int bytes=length*_channels*bytes_per_sample;

fseek(fd,40,SEEK_SET);
fwrite(&bytes,4,1,fd);

bytes=length*_channels*bytes_per_sample+44;

fseek(fd,4,SEEK_SET);
fwrite(&bytes,4,1,fd);
fclose(fd);
}


void _normalize(){
float maxf=0.0f;
for(int i=0;i<_length;i++){
   if(maxf<sqrt(_data[i]*_data[i])) maxf=sqrt(_data[i]*_data[i]);
}
maxf=1.0f/maxf;
for(int i=0;i<_length;i++){
   _data[i]*=maxf;
}
}

void _ready_cache(int begin,int size){
if(!_opened) return;
if(begin>=_length) return;

if((begin+size)>=_length) size=_length-begin;

if(_cache){
  if(size>_cache_size){
    free(_cache);
	_cache=(float*)malloc(_channels*size*sizeof(float));
  }
}else _cache=(float*)malloc(_channels*_length*sizeof(float));

memset(_cache,0,_channels*size*sizeof(float));

_seek(begin);

_read(size,_cache);

_close();

}


void _copy(unsigned int begin,unsigned int end,CWaveData* data){
data->_fd=_fd;
data->_opened=true;
data->_big_endian=_big_endian;
data->_dataOffset=_dataOffset;
data->_dataType=_dataType;
data->_format_tag=_format_tag;
data->_current_pos=0;
data->_on_memory=true;
data->_channels=_channels;
data->_length=begin-end;
data->_sampleRate=_sampleRate;
if(data->_data) free(data->_data);
data->_data=(float*)malloc(_channels*data->_length*sizeof(float));
if(_on_memory){
  memcpy(data->_data,&_data[_channels*begin],_channels*(begin-end)*sizeof(float));
}else{
  _seek(begin);
  _read(begin-end,data->_data);
}
}


void _rec_start(){}
void _rec_end(){ _render_Recdata(); }

void _record(int size,float* buf){
if(_tail_frame){
  if((_current_pos+size)>=_rec_frameSize){
	int remain=_rec_frameSize-_current_pos;
	memcpy(&_tail_frame->_buf[_current_pos],buf,_channels*remain*sizeof(float));
	RECORD_FRAME* tmp=new RECORD_FRAME;
	tmp->_buf=(float*)malloc(_rec_frameSize*sizeof(float));
	memset(tmp->_buf,0,_rec_frameSize*sizeof(float));
	tmp->_next=NULL;
	tmp->_prev=_tail_frame;
	_tail_frame->_next=tmp;
	_tail_frame=tmp;
	if((size-remain)>0) memcpy(_tail_frame->_buf,&buf[remain],_channels*(size-remain)*sizeof(float));
	_current_pos=size-remain;
	_num_frames++;
  }else{
    memcpy(&_tail_frame->_buf[_current_pos],buf,_channels*size*sizeof(float));
	_current_pos+=size;
  }
}else{
  _rec_frame=new RECORD_FRAME;
  _rec_frame->_buf=(float*)malloc(_rec_frameSize*sizeof(float));
  memset(_rec_frame->_buf,0,_rec_frameSize*sizeof(float));
  _rec_frame->_prev=_rec_frame->_next=NULL;
  _tail_frame=_rec_frame;
  memcpy(_rec_frame->_buf,buf,_channels*size*sizeof(float));
  _current_pos=size;
  _num_frames=1;
}
}

void _render_Recdata(){
if(_num_frames==0) return;
if(_data) free(_data);
_length=(_num_frames-1)*_rec_frameSize+_current_pos;
_data=(float*)malloc(_channels*_length*sizeof(float));

int idx=0;;
RECORD_FRAME* tmp=_rec_frame;
while(tmp){
     if(tmp==_tail_frame) break;
	 memcpy(&_data[_channels*idx*_rec_frameSize],tmp->_buf,_channels*_rec_frameSize*sizeof(float));
	 idx++;
     tmp=tmp->_next;
}
memcpy(&_data[_channels*idx*_rec_frameSize],tmp->_buf,_channels*_current_pos*sizeof(float));

tmp=_rec_frame;
while(tmp){
     RECORD_FRAME* tmpp=tmp->_next;
	 free(tmp->_buf);
	 delete tmp;
	 tmp=tmpp;
}
_rec_frame=_tail_frame=NULL;
_num_frames=0;


}



};




#endif
