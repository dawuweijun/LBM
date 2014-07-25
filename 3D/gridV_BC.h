/*********************************************************************************************
 * Copyright (c) <2014>, <Santosh Ansumali@JNCASR>                                           *
 *  All rights reserved.                                                                     *
 *   Redistrnode.iBution and use in source and binary forms, with or without modification, are    *
 *   permitted provided that the following conditions are met:                               *
 *                                                                                           *
 *    1. Redistrnode.iButions of source code must retain the above copyright notice, this list of *
 *       conditions and the following disclaimer.                                            *
 *    2. Redistrnode.iButions in binary form must reproduce the above copyright notice, this list *
 *       of conditions and the following disclaimer in the documentation and/or other        *
 *       materials provided with the distrnode.iBution.                                           *
 *    3. Neither the name of the <JNCASR> nor the names of its contrnode.iButors may be used to   *
 *       endorse or promote products derived from this software without specific prior       *
 *       written permission.                                                                 *
 *                                                                                           *
 *       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRnode.iBUTORS "AS IS" AND     *
 *       ANY EXPRESS OR IMPLnode.iED WARRANTnode.iES, INCLUDING, BUT NOT LIMITED TO, THE IMPLnode.iED       *
 *       WARRANTnode.iES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  *
 *       IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRnode.iBUTORS BE LIABLE FOR ANY DIRECT,    *
 *       INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,      *
 *       BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,       *
 *       DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF     *
 *       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE     *
 *       OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED   *
 *       OF THE POSSnode.iBILITY OF SUCH DAMAGE.                                                  *
 *                                                                                           *
 *       Suggestions:          ansumali@jncasr.ac.in                                         *
 *       Bugs:                 ansumali@jncasr.ac.in                                         *
 *                                                                                           *
 *********************************************************************************************/
#ifndef _GRIDV_BC_H_
#define  _GRIDV_BC_H_
#include<iostream>
#include "gridV.h"
#include "globals3D.h"
#include <unistd.h>
 

/***************************************************************************************************************
 * function to  set periodic BC in single Dir                                                               * 
 *                                                                                                             *
 * @Author: Santosh Ansumali                                                                                   *
 * Last Change: 05/02/2014                                                                                     *
 * *************************************************************************************************************/
template<int N>
void periodicNode1(gridV<N> &node){
  // Now set periodic in   Left end of 1
	for (int i3 = 0; i3 < node.m3; i3++) {
		for (int  i2=0; i2< node.m2;  i2++) {
			for(int i1BDummy=0; i1BDummy<node.nB1; i1BDummy++)
				for(int k=0;k<N;k++) 
					node(i1BDummy,i2,i3,k) = node(node.iE1+i1BDummy,i2,i3,k); 
//					node(0,i2,i3,k) = node(node.iE1,i2,i3,k);					
		}
	}
//Right end of 1
	for (int i3 = 0; i3 < node.m3; i3++) {
		for (int  i2=0; i2<node.m2;  i2++) {
			int left= node.iB1;
			for(int i1EDummy=node.iE1+1; i1EDummy<node.m1; i1EDummy++,left++)
				for(int k=0;k<N;k++) 
					node(i1EDummy,i2,i3,k)  = node(left,i2,i3,k); 
//					node(node.iE1+1,i2,i3,k)  = node(node.iB1,i2,i3,k); 
				
			  
		}
	}
 
}
template<int N>
void periodicNode2(gridV<N> &node){
	for (int i3 = 0; i3 < node.m3; i3++) {
		for (int  i1=0; i1< node.m1;  i1++) {
			for(int i2BDummy=0; i2BDummy<node.nB2; i2BDummy++)
				for(int k=0;k<N;k++) 
					node(i1,i2BDummy,i3,k)  = node(i1,node.iE2+i2BDummy,i3,k);
//					node(i1,0,i3,k)  = node(i1,node.iE2,i3,k);  //Debug
		}
	}
	for (int i3 = 0; i3 < node.m3; i3++) {
		for (int  i1=0; i1<node.m1;  i1++) {
			int left= node.iB2;
			for(int i2EDummy=node.iE2+1; i2EDummy<node.m2; i2EDummy++,left++)
				for(int k=0;k<N;k++) 
					node(i1,i2EDummy,i3,k)  = node(i1,left,i3,k); 
//					node(i1,node.iE2+1,i3,k)  = node(i1,1,i3,k);
		}
	}

 
}

template<int N>
void periodicNode3(gridV<N> &node){
	for (int i1 = 0; i1 < node.m1; i1++) {
		for (int  i2=0; i2< node.m2;  i2++) {
			for(int i3BDummy=0; i3BDummy<node.nB3; i3BDummy++)
				for(int k=0;k<N;k++) 
					node(i1,i2,i3BDummy,k)  = node(i1,i2,node.iE3+i3BDummy,k); 
//					node(i1,i2,0,k)  = node(i1,i2,node.iE3,k);
		}
	}

	for (int i1 = 0; i1 < node.m1; i1++) {
		for (int  i2=0; i2< node.m2;  i2++) {
			int left= node.iB3;
			for(int i3EDummy=node.iE3+1; i3EDummy<node.m3; i3EDummy++,left++)
			  	for(int k=0;k<N;k++) 
					node(i1,i2,i3EDummy,k)  = node(i1,i2,left,k); 
//					node(i1,i2,node.iE3+1,k)  = node(i1,i2,1,k);
			  
		}
	}
 
}


template<int N>
void periodicNode12(gridV<N> &node){
   periodicNode1( node);
    periodicNode2( node);
 
}

template<int N>
void periodicNode13(gridV<N> &node){
  periodicNode1( node);
    periodicNode3( node);
}

 

template<int N>
void periodicNode23(gridV<N> &node){
  periodicNode2( node);
    periodicNode3( node);
}

template<int N>
void periodicNode123(gridV<N> &node){
  periodicNode3( node);
	 periodicNode2( node);
	 periodicNode1( node); 
}

//Now set periodic for all the 7 distributions
void periodic1(lbgrid& node) {
	periodicNode1(node.SCP);
	periodicNode1(node.SCM);
	periodicNode1(node.FCC12);
	periodicNode1(node.FCC13);
	periodicNode1(node.FCC23);
	periodicNode1(node.BCCP);
	periodicNode1(node.BCCM);
}
void periodic2(lbgrid& node) {
	periodicNode2(node.SCP);
	periodicNode2(node.SCM);
	periodicNode2(node.FCC12);
	periodicNode2(node.FCC13);
	periodicNode2(node.FCC23);
	periodicNode2(node.BCCP);
	periodicNode2(node.BCCM);
}
void periodic3(lbgrid& node) {
	periodicNode3(node.SCP);
	periodicNode3(node.SCM);
	periodicNode3(node.FCC12);
	periodicNode3(node.FCC13);
	periodicNode3(node.FCC23);
	periodicNode3(node.BCCP);
	periodicNode3(node.BCCM);
}
void periodic12(lbgrid& node) {
	periodicNode12(node.SCP);
	periodicNode12(node.SCM);
	periodicNode12(node.FCC12);
	periodicNode12(node.FCC13);
	periodicNode12(node.FCC23);
	periodicNode12(node.BCCP);
	periodicNode12(node.BCCM);
}
void periodic23(lbgrid& node) {
	periodicNode23(node.SCP);
	periodicNode23(node.SCM);
	periodicNode23(node.FCC12);
	periodicNode23(node.FCC13);
	periodicNode23(node.FCC23);
	periodicNode23(node.BCCP);
	periodicNode23(node.BCCM);
}
void periodic13(lbgrid& node) {
	periodicNode13(node.SCP);
	periodicNode13(node.SCM);
	periodicNode13(node.FCC12);
	periodicNode13(node.FCC13);
	periodicNode13(node.FCC23);
	periodicNode13(node.BCCP);
	periodicNode13(node.BCCM);
}
void periodic123(lbgrid& node) {
	periodicNode123(node.SCP);
	periodicNode123(node.SCM);
	periodicNode123(node.FCC12);
	periodicNode123(node.FCC13);
	periodicNode123(node.FCC23);
	periodicNode123(node.BCCP);
	periodicNode123(node.BCCM);
}

//prepare BC
template<int N>
void prepareBCNode1(gridV<N> &node) {  
	int n1 = node.n1;
	for (int i3 = 0; i3 < node.m3; i3++) {
		for (int  i2=0; i2< node.m2;  i2++) {
			for(int iDummy = 0; iDummy < node.nB1; iDummy++) {
				for (int k = 0; k < N; k++) {
					node(iDummy,i2,i3,k) = node(node.iB1-iDummy,i2,i3,k);  //left end
					node(node.iE1+iDummy+1,i2,i3,k) = node(node.iE1+iDummy,i2,i3,k); //right end
				}
			}
		}
	}
	
}

template<int N>
void prepareBCNode2(gridV<N> &node) {  
	int n2 = node.n2;
	for (int i3 = 0; i3 < node.m3; i3++) {
		for (int  i1=0; i1< node.m1;  i1++) {
			for(int iDummy = 0; iDummy < node.nB2; iDummy++) {
				for (int k = 0; k < N; k++) {
					node(i1,iDummy,i3,k) = node(i1,node.iB2-iDummy,i3,k);
					node(i1,node.iE2+iDummy+1,i3,k) = node(i1,node.iE2+iDummy,i3,k);
				}
			}
		}
	}
	
}

template<int N>
void prepareBCNode3(gridV<N> &node) {  
	for (int i2 = 0; i2 < node.m2; i2++) {
		for (int  i1=0; i1< node.m1;  i1++) {
			for(int iDummy = 0; iDummy < node.nB3; iDummy++) {
				for (int k = 0; k < N; k++) {
					node(i1,i2,iDummy,k) = node(i1,i2,node.iB3-iDummy,k);
					node(i1,i2,node.iE3+iDummy+1,k) = node(i1,i2,node.iE3+iDummy,k);
				}
			}
		}
	}
	
}

template<int N>
void prepareBCNode12(gridV<N> &node) { 
	prepareBCNode1(node);
	prepareBCNode2(node);
}

template<int N>
void prepareBCNode13(gridV<N> &node) { 
	prepareBCNode1(node);
	prepareBCNode3(node);
}

template<int N>
void prepareBCNode23(gridV<N> &node) { 
	prepareBCNode3(node);
	prepareBCNode2(node);
}

template<int N>
void prepareBCNode123(gridV<N> &node) { 
	prepareBCNode3(node);
	prepareBCNode2(node);
	prepareBCNode1(node);
}

void prepareBC1(lbgrid &node) {
	prepareBCNode1(node.SCP);
	prepareBCNode1(node.SCM);
	prepareBCNode1(node.FCC12);
	prepareBCNode1(node.FCC13);
	prepareBCNode1(node.FCC23);
	prepareBCNode1(node.BCCP);
	prepareBCNode1(node.BCCM);
}
void prepareBC2(lbgrid& node) {
	prepareBCNode2(node.SCP);
	prepareBCNode2(node.SCM);
	prepareBCNode2(node.FCC12);
	prepareBCNode2(node.FCC13);
	prepareBCNode2(node.FCC23);
	prepareBCNode2(node.BCCP);
	prepareBCNode2(node.BCCM);
}
void prepareBC3(lbgrid& node) {
	prepareBCNode3(node.SCP);
	prepareBCNode3(node.SCM);
	prepareBCNode3(node.FCC12);
	prepareBCNode3(node.FCC13);
	prepareBCNode3(node.FCC23);
	prepareBCNode3(node.BCCP);
	prepareBCNode3(node.BCCM);
}
void prepareBC12(lbgrid& node) {
	prepareBCNode12(node.SCP);
	prepareBCNode12(node.SCM);
	prepareBCNode12(node.FCC12);
	prepareBCNode12(node.FCC13);
	prepareBCNode12(node.FCC23);
	prepareBCNode12(node.BCCP);
	prepareBCNode12(node.BCCM);
}
void prepareBC23(lbgrid& node) {
	prepareBCNode23(node.SCP);
	prepareBCNode23(node.SCM);
	prepareBCNode23(node.FCC12);
	prepareBCNode23(node.FCC13);
	prepareBCNode23(node.FCC23);
	prepareBCNode23(node.BCCP);
	prepareBCNode23(node.BCCM);
}
void prepareBC13(lbgrid& node) {
	prepareBCNode13(node.SCP);
	prepareBCNode13(node.SCM);
	prepareBCNode13(node.FCC12);
	prepareBCNode13(node.FCC13);
	prepareBCNode13(node.FCC23);
	prepareBCNode13(node.BCCP);
	prepareBCNode13(node.BCCM);
}
void prepareBC123(lbgrid& node) {
	prepareBCNode123(node.SCP);
	prepareBCNode123(node.SCM);
	prepareBCNode123(node.FCC12);
	prepareBCNode123(node.FCC13);
	prepareBCNode123(node.FCC23);
	prepareBCNode123(node.BCCP);
	prepareBCNode123(node.BCCM);
}

#endif
