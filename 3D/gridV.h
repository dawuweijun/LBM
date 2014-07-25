/*********************************************************************************************
 * Copyright (c) <2014>, <Santosh Ansumali@JNCASR>                                           *
 *  All rights reserved.                                                                     *
 *   Redistribution and use in source and binary forms, with or without modification, are    *
 *   permitted provided that the following conditions are met:                               *
 *                                                                                           *
 *    1. Redistributions of source code must retain the above copyright notice, this list of *
 *       conditions and the following disclaimer.                                            *
 *    2. Redistributions in binary form must reproduce the above copyright notice, this list *
 *       of conditions and the following disclaimer in the documentation and/or other        *
 *       materials provided with the distribution.                                           *
 *    3. Neither the name of the <JNCASR> nor the names of its contributors may be used to   *
 *       endorse or promote products derived from this software without specific prior       *
 *       written permission.                                                                 *
 *                                                                                           *
 *       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND     *
 *       ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED       *
 *       WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  *
 *       IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,    *
 *       INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,      *
 *       BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,       *
 *       DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF     *
 *       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE     *
 *       OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED   *
 *       OF THE POSSIBILITY OF SUCH DAMAGE.                                                  *
 *                                                                                           *
 *       Suggestions:          ansumali@jncasr.ac.in                                         *
 *       Bugs:                 ansumali@jncasr.ac.in                                         *
 *                                                                                           *
 *********************************************************************************************/
#ifndef _GRID_V_H_
#define _GRID_V_H_
#include<iostream>
using namespace std;
 
typedef double myReal;
template<int N>
class gridV{
public:
gridV(int _n1, int _n2, int _numPad1, int _numPad2){
	n1= _n1;
	n2 = _n2;
	n3= 1;
	nB1= _numPad1;
	nB2=_numPad2;
	nB3= 0;
	m1= _n1+2*_numPad1;
	m2 = _n2+2*_numPad2;
	m3=  1;
	iB1 = _numPad1;
	iB2 = _numPad2;
	iB3 = 1;
	iE1= _n1+_numPad1-1;
	iE2 = _n2+_numPad2-1;
	iE3=   1;
	size = m1*m2*m3;
	data = new myReal[size*N];
}
gridV(int _n1, int _n2,int _n3, int _numPad1, int _numPad2, int _numPad3){
	n1= _n1;
	n2 = _n2;
	n3= _n3;
	nB1= _numPad1;
	nB2= _numPad2;
	nB3= _numPad3;		
	m1= _n1+2*_numPad1;
	m2 = _n2+2*_numPad2;
	m3 = _n3+2*_numPad3;
	
	iB1 = _numPad1;
	iB2 = _numPad2;
	iB3 = _numPad3;
	iE1= _n1+_numPad1-1;
	iE2 = _n2+_numPad2-1;
	iE3 = _n3+_numPad3-1;
	size = m1*m2*m3;
	data = new myReal[size*N];
}

myReal operator() (int i, int j, int k,int dv) const {
	return data[k*m1*m2*N+j*m1*N+i*N+dv];
}
myReal& operator() (int i, int j, int k,int dv)  
{
return data[k*m1*m2*N+j*m1*N+i*N+dv];
}
//Use as 2D array;
 myReal operator() (int i, int j, int dv) const
{
return data[j*m1*N+i*N+dv];
}
myReal& operator() (int i, int j, int dv)  
{
return data[j*m1*N+i*N+dv];
}
//use as 1D Array
myReal operator[] (int i) const
{
return data[i];
}
myReal& operator[] (int i)  
{
return data[i];
}
~gridV(){}
//USer Specified Sizes
int n1;
int n2;
int n3;

  
// Num of Ghost Nodes at any one end
int nB1;
int nB2;
int nB3;
//Total Sizes
int m1;
int m2;
int m3;
//Beginning of actual nodes - physical nodes
int iB1;
int iB2;
int iB3;
//End point for actual nodes
int iE1;
int iE2;
int iE3;

int size;
myReal *data; 
};

//Class to hold node type and cell type
class gridtype {
	public:
		gridtype(int _n1, int _n2,int _n3, int _numPad1, int _numPad2, int _numPad3) {
			m1= _n1+2*_numPad1;
			m2 = _n2+2*_numPad2;
			m3 = _n3+2*_numPad3;
			nB1= _numPad1;
			nB2= _numPad2;
			nB3= _numPad3;
			iB1 = _numPad1;
			iB2 = _numPad2;
			iB3 = _numPad3;
			iE1= _n1+_numPad1-1;
			iE2 = _n2+_numPad2-1;
			iE3 = _n3+_numPad3-1;
			size = m1*m2*m3;
			data = new int[size];
 			for (int i = 0; i < size; i++) data[i] = 1;
		}
		int n1,n2,n3,m1,m2,m3,size;
		int iE1,iE2,iE3,iB1,iB2,iB3;
		int nB1,nB2,nB3;
		int *data;
		//Overload the () operator
		int operator() (int i, int j, int k) const {
			return data[k*m1*m2 + j*m1 + i];
		}
		int& operator() (int i, int j, int k) {
			return data[k*m1*m2 + j*m1 + i];
		}
};

//				
class arrayRT {	
//				
	public:
		arrayRT(int _n1, int _n2,int _n3,int _numPad1, int _numPad2, int _numPad3) {
			m1= _n1+2*_numPad1;
			m2 = _n2+2*_numPad2;
			m3 = _n3+2*_numPad3;
			size = m1*m2*m3;
			data = new myReal[size];
 			for (int i = 0; i < size; i++) data[i] = 1.0;
		}
		int m1,m2,m3,size;
		myReal *data;
		myReal operator() (int i, int j, int k) const {
			return data[k*m1*m2 + j*m1 + i];
		}
		myReal& operator() (int i, int j, int k) {
			return data[k*m1*m2 + j*m1 + i];
		}
};
 
/***************************************************************************************************************
 * function to  print dual grid strucure                                                                            * 
 *                                                                                                             *
 * @Author: Santosh Ansumali                                                                                   *
 * Last Change: 05/02/2014                                                                                     *
 * *************************************************************************************************************/
template<int N>
void printNodeCell(gridV<N> node, gridV<N> cell, std::ostream &file1, int k) {
	for (int i3 = 1; i3 < node.n3; i3++) {
		for (int i2 = 1; i2 < node.n2; i2++) {
			for(int i1=1; i1<= node.n1; i1++) {
				file1 << cell(i1,i2,i3,k) << " ";
			}
			file1 << std::endl;
			for(int i1=1; i1<= node.n1; i1++) {
				file1 << node(i1,i2,i3,k) << " ";
			}
			file1 << std::endl;
		}
		file1 << std::endl;
	}
}


template<int N>
std::istream& operator>> (std::istream& os, gridV<N> node) { 
	for (int i3 = node.n3; i3 >=1; i3--) {	
		for(int i2=node.n2; i2>= 1; i2--) { 
			for(int i1=1; i1<= node.n1; i1++) {
				for(int dv =0; dv<N;dv++)   
					os>>node(i1,i2,i3,dv);
			}
		}
	}
   	return os;
}

   	
template<int N>
std::ostream& operator<< (std::ostream& os, gridV<N> node) { 
	for (int i3 = node.n3; i3 >=1; i3--) {	
		for(int i2=node.n2; i2>= 1; i2--) { 
			for(int i1=1; i1<= node.n1; i1++) {
				for(int dv =0; dv<N;dv++)     
					os<<node(i1,i2,i3,dv);
			}
		}
	}
   	return os;
}

 

/***************************************************************************************************************
 * function to  print   grid strucure                                                                            * 
 *                                                                                                             *
 * @Author: Santosh Ansumali                                                                                   *
 * Last Change: 05/02/2014                                                                                     *
 * *************************************************************************************************************/
 template<int N>
void printNode(gridV<N> node,std::ostream &file1, int k) {
	for (int i3 = 1; i3 < node.n3; i3++) {
		for (int i2 = 1; i2 < node.n2; i2++) {
			for(int i1=1; i1<= node.n1; i1++) {
				file1 << node(i1,i2,i3,k) << " ";
			}
			file1 << std::endl;
		}
		file1 << std::endl;
	}
}

template<int N>
void printNode(gridV<N> node,std::ostream &file1) {
	for (int i3 = 1; i3 < node.n3; i3++) {
		for (int i2 = 1; i2 < node.n2; i2++) {
			for(int i1=1; i1<= node.n1; i1++) {
				for (int k = 0; k < N; k++)
					file1 << node(i1,i2,i3,k) << " ";
			}
			file1 << std::endl;
		}
		file1 << std::endl;
	}
} 									


#endif
