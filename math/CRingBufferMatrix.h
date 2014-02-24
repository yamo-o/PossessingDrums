

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


#include "CDenseMatrix.h"

#ifndef CRINGBUFFERMAT_HH
#define CRINGBUFFERMAT_HH


template <class T> class CRingBufferMatrix:public CDenseMatrix<T>{
public:
unsigned int _head,_tail;

CRingBufferMatrix():CDenseMatrix<T>(){}

CRingBufferMatrix(int row,int col):CDenseMatrix<T>(row,col){
_head=0;
_tail=row-1;
}

virtual inline T &operator()(int i,int j){
int idx=_head+i;
if(idx>=this->_row) idx-=this->_row;
return(this->_value[idx][j]);
}

virtual inline const T &operator()(int i,int j) const {
int idx=_head+i;
if(idx>=this->_row) idx-=this->_row;
return(this->_value[idx][j]);
}

virtual inline T* operator[](int i){
int idx=_head+i;
if(idx>=this->_row) idx-=this->_row;
return(this->_value[idx]);
}

virtual inline const T* operator[](int i) const {
int idx=_head+i;
if(idx>=this->_row) idx-=this->_row;
return(this->_value[idx]);
}

void _resize(unsigned int row){
_head=0;
_tail=row-1;
CDenseMatrix<T>::_resize(row);
}

void _resize(unsigned int row,unsigned int col){
_head=0;
_tail=row-1;
CDenseMatrix<T>::_resize(row,col);
}

inline void _copy_vec(int row,T* vec){
int idx=_head+row;
if(idx>=this->_row) idx-=this->_row;
memcpy(this->_value[idx],vec,this->_col*sizeof(T));
}

inline void _zeroRow(int row){
int idx=_head+row;
if(idx>=this->_row) idx-=this->_row;
memset(this->_value[idx],0,this->_col*sizeof(T));
}

inline void _advance(){
++_head;
++_tail;
if(_tail>=this->_row){
  _tail=0;
}else if(_head>=this->_row){
  _head=0;
}
}

inline void _push(T* vec){
_advance();
memcpy(this->_value[_tail],vec,this->_col*sizeof(T));
}

};


#endif
