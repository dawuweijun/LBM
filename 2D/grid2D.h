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
#ifndef _GRID2D_H_
#define _GRID2D_H_
#include<iostream>
enum marker{SOLID, FLUID};
 
 typedef  double myReal;
template<int N>
class grid2D{
public:
grid2D(int _nX, int _nY){
size = (_nX+2)*(_nY+2);
nX = _nX;
nY = _nY;
data = new myReal[size*N];
}
 
myReal operator() (int i, int j, int k) const
{
return data[i*(nY+2)*N+j*N+k];
}

myReal& operator() (int i, int j, int k)
{
return  (data[i*(nY+2)*N+j*N+k]);
}

~grid2D(){}
int size;
int nX;
int nY;
//marker type;
myReal *data; 
};

class GridType {
	public:
		int *type,size,nX,nY;
		GridType(int _nX, int _nY) {
			size = (_nX+2)*(_nY+2);
			nX = _nX;
			nY = _nY;
			type = new int[size];
		}
		int operator() (int i, int j) const {
			return type[i*(nY+2) + j]	;
		}
		int& operator() (int i, int j) {
			return type[i*(nY+2) + j];
		}
		
};


/***************************************************************************************************************
 * function to  set Periodic BC                                                                                * 
 *                                                                                                             *
 * @Author: Santosh Ansumali                                                                                   *
 * Last Change: 05/02/2014                                                                                     *
 * *************************************************************************************************************/
template<int N>
void periodic(grid2D<N> node){
int i2BDummy= 0;
int i2EDummy= node.nY+1;
// Set periodic in Y first 
for (int i1=0; i1< node.nX+2;  i1++)
{
	for(int k=0;k<N;k++) 
	node(i1,i2BDummy,k)  = node(i1,node.nY,k); 
	for(int k=0;k<N;k++) 
	node(i1,i2EDummy,k)  = node(i1,1,k); 
}
// Now set periodic in X
int i1BDummy= 0;
int i1EDummy= node.nY+1; 
for (int  i2=0; i2<  node.nY+2;  i2++)
{
	for(int k=0;k<N;k++) 
	node(i1BDummy,i2,k)  = node(node.nX,i2,k); 
	for(int k=0;k<N;k++) 
	node(i1EDummy,i2,k)  = node(1,i2,k); 
}
}


/***************************************************************************************************************
 * function to  set Periodic BC   in X                                                                         * 
 *                                                                                                             *
 * @Author: Santosh Ansumali                                                                                   *
 * Last Change: 05/02/2014                                                                                     *
 * *************************************************************************************************************/
template<int N>
void periodicX(grid2D<N> node){ 
// Now set periodic in X
int i1BDummy= 0;
int i1EDummy= node.nY+1; 
for (int  i2=0; i2<  node.nY+2;  i2++)
{
	for(int k=0;k<N;k++) 
	node(i1BDummy,i2,k)  = node(node.nX,i2,k); 
	for(int k=0;k<N;k++) 
	node(i1EDummy,i2,k)  = node(1,i2,k); 
}
}



/***************************************************************************************************************
 * function to  set Periodic BC in Y                                                                           * 
 *                                                                                                             *
 * @Author: Santosh Ansumali                                                                                   *
 * Last Change: 05/02/2014                                                                                     *
 * *************************************************************************************************************/
template<int N>
void periodicY(grid2D<N> node){
int i2BDummy= 0;
int i2EDummy= node.nY+1;
// Set periodic in Y first 
for (int i1=0; i1< node.nX+2;  i1++)
{
	for(int k=0;k<N;k++) 
	node(i1,i2BDummy,k)  = node(i1,node.nY,k); 
	for(int k=0;k<N;k++) 
	node(i1,i2EDummy,k)  = node(i1,1,k); 
} 
}

/***************************************************************************************************************
 * function to  get grid Laplacian using CD2                                                                   * 
 *                                                                                                             *
 * @Author: Santosh Ansumali                                                                                   *
 * Last Change: 05/02/2014                                                                                     *
 * *************************************************************************************************************/
template<int N>
void Laplacian_CD2(grid2D<N> node, grid2D<N> nodeLap, myReal dxInv){
for (int  i2=1; i2<=  node.nY;  i2++)
{ 
for (int i1=0; i1<= node.nX;  i1++)
{
  for(int k=0;k<N; k++)
  nodeLap(i1,i2,k) = dxInv*dxInv*(node(i1+1,i2,k)+node(i1-1,i2,k)+node(i1,i2+1,k)+node(i1,i2-1,k)-4.0*node(i1,i2,k)); 
} 
}
}

template<int N> 
void prepareBC(grid2D<N> &node)
{
   int nY= node.nY;
   int nX= node.nX;
   
//for(int i1=node.nX; i1>= 1; i1--)
//  {
//    for(int dv =0;dv<N;dv++)
//      {
//	node(i1,nY+1, dv) = node(i1,nY, dv); 
//	  node(i1,0, dv) = node(i1,1, dv);
//	  }
//	}
for(int i2=node.nY; i2>= 1; i2--)
	{
	  for(int dv =0;dv<N;dv++)
	  { 
	  node(nX+1,i2, dv) = node(nX,i2, dv); 
	  node(0,i2, dv) = node(1,i2, dv);
	  }
	}

//Corners
for (int dv = 0; dv<N; dv++) {
	node(nX+1,0,dv) = node(nX,1,dv);
	node(nX+1,nY+1,dv) = node(nX,nY,dv);
	node(0,0,dv) = node(1,1,dv);
	node(0,nY+1,dv) = node(1,nY,dv);
}
}


template<int N> 
void prepareBCY(grid2D<N> &node, grid2D<N> &cell)
{
   int nY= node.nY;
  for(int i1=node.nX; i1>= 1; i1--)
	{
	  for(int dv =0;dv<N;dv++){
	  cell(i1,nY+1, dv) = cell(i1,nY, dv); 
	  node(i1,nY+1, dv) = node(i1,nY, dv); 
	  }
	}
	 for(int i1=node.nX; i1>= 1; i1--)
	{
	  for(int dv =0;dv<N;dv++){
	  cell(i1,0, dv) = cell(i1,1, dv); 
	  node(i1,0, dv) = node(i1,1, dv); 
	  }
	}
}
/***************************************************************************************************************
 * function to  print dual grid strucure                                                                            * 
 *                                                                                                             *
 * @Author: Santosh Ansumali                                                                                   *
 * Last Change: 05/02/2014                                                                                     *
 * *************************************************************************************************************/
 template<int N>
void printNodeCell(grid2D<N> node, grid2D<N>  cell,std::ostream &file1, int k)
{

for(int i2=node.nY; i2>= 1; i2--)
{
file1 << " ";
	for(int i1=1; i1<= node.nX; i1++)
	{
	file1 << cell(i1,i2,k) << " ";
	}
	
file1 << std::endl;
	for(int i1=1; i1<= node.nX; i1++)
	{
	file1 << node(i1,i2,k) << " ";
	}

file1 << std::endl;
}

}


template<int N>
	std::istream& 
	operator>> (std::istream& os, grid2D<N> nodeField){ 
	 for(int i2=nodeField.nY; i2>= 1; i2--)
{ 
  
	for(int i1=1; i1<= nodeField.nX; i1++)
	{
	  for(int dv =0; dv<N;dv++)   
	os>>nodeField(i1,i2,dv);
	}
}
   	return os;
   	}

   	
   	template<int N>
	std::ostream& 
	operator<< (std::ostream& os, grid2D<N> nodeField){ 
	 for(int i2=nodeField.nY; i2>= 1; i2--)
{ 
  
	for(int i1=1; i1<= nodeField.nX; i1++)
	{
	  for(int dv =0; dv<N;dv++)   
	os<<nodeField(i1,i2,dv);
	}
}
   	return os;
   	}
	
template<int N>
void interpolateCellFromNode(grid2D<N> nodeField, grid2D<N>  cellField)
{

for(int i2=nodeField.nY; i2>= 1; i2--)
{ 
  
	for(int i1=1; i1<= nodeField.nX; i1++)
	{
	  for(int dv =0; dv<N;dv++)
	  cellField(i1,i2,dv) =0.25*(nodeField(i1,i2,dv) +nodeField(i1+1,i2,dv)+nodeField(i1,i2+1,dv)+nodeField(i1+1,i2+1,dv));
	}
	
  
}

}


#endif
