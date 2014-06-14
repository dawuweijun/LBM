#ifndef _GRIDV_PACK_H_
#define  _GRIDV_PACK_H_
#include<iostream>
#include<stdlib.h>

#include "gridV.h"
#include "gridV_Probe.h"
#include "globals3D.h"

//Grid spit in 1-direction: Pack in 23 plane
//Allocate size
myReal* allocateForPacking(int n1, int n2, int n3) {
	myReal *temp;
	//sort in descending order
	if (n1<n2) swap(n1,n2);
	if (n2<n3) swap(n2,n3);
	if (n1<n2) swap(n1,n2);
	temp = (double*) malloc(n1*n2*27);
}
	
	
template<int N> void packToPlane23(lbgrid &node,myReal *temp1P, myReal *temp1M) {
//	int size23 = node.SCP.n2*node.SCP.n3*(2+2+2+2+1);
	int size = node.SCP.n2*node.SCP.n3;
//	myReal temp1P[size23];
//	myReal temp1M[size23];
	int index = 0;
	int stride = 1;
	//Pack at the end only the positive distributions in i1 direction
	//Pack BCC
	//Pack BCCP
	probe23(node.BCCP, temp1P, node.BCCP.iE1, DV_P1_P1_P1, stride, index);
	index += size;
	probe23(node.BCCP, temp1P, node.BCCP.iE1, DV_P1_M1_P1, stride, index);
	index += size;
	//Pack BCCM
	probe23(node.BCCM, temp1P, node.BCCM.iE1, DV_P1_P1_M1, stride, index);
	index += size;
	probe23(node.BCCM, temp1P, node.BCCM.iE1, DV_P1_M1_M1, stride, index);
	index += size;
	//Pack SC
	//Pack SCP
	probe23(node.SCP, temp1P, node.SCP.iE1, DV_P1_ZERO_ZERO, stride, index);
	index += size;
	//Pack FCC
	//Pack FCC12
	probe23(node.FCC12, temp1P, node.FCC12.iE1, DV_P1_P1_ZERO, stride, index);
	index += size;
	probe23(node.FCC12, temp1P, node.FCC12.iE1, DV_P1_M1_ZERO, stride, index);
	index += size;
	//Pack FCC13
	probe23(node.FCC13, temp1P, node.FCC13.iE1, DV_P1_ZERO_P1, stride, index);
	index += size;
	probe23(node.FCC13, temp1P, node.FCC13.iE1, DV_P1_ZERO_M1, stride, index);
	index += size;
	//Pack FCC23
	//No need: ZERO in i1 direction 
	
	//Pact at the beginning only the negative distributions in i1 directions
	index = 0;
	//Pack BCC
	//Pack BCCP
	probe23(node.BCCP, temp1M, node.BCCP.iB1, DV_M1_P1_P1, stride, index);
	index += size;
	probe23(node.BCCP, temp1M, node.BCCP.iB1, DV_M1_M1_P1, stride, index);
	index += size;
	//Pack BCCM
	probe23(node.BCCM, temp1M, node.BCCM.iB1, DV_M1_P1_M1, stride, index);
	index += size;
	probe23(node.BCCM, temp1M, node.BCCM.iB1, DV_M1_M1_M1, stride, index);
	index += size;
	//Pack SC
	//Pack SCM
	probe23(node.SCM, temp1M, node.SCM.iB1, DV_M1_ZERO_ZERO, stride, index);
	index += size;
	//Pack FCC
	//Pack FCC12
	probe23(node.FCC12, temp1M, node.FCC12.iB1, DV_M1_P1_ZERO, stride, index);
	index += size;
	probe23(node.FCC12, temp1M, node.FCC12.iB1, DV_M1_M1_ZERO, stride, index);
	index += size;
	//Pack FCC13
	probe23(node.FCC13, temp1M, node.FCC13.iB1, DV_M1_ZERO_P1, stride, index);
	index += size;
	probe23(node.FCC13, temp1M, node.FCC13.iB1, DV_M1_ZERO_M1, stride, index);
	index += size;
	//Pack FCC23
	//No need: ZERO in i1 direction 
}



//Grid spit in 2-direction: Pack in 13 plane
template<int N> void packToPlane13(lbgrid &node, myReal *temp2P, myReal *temp2M) {
	//int size13 = node.SCP.n1*node.SCP.n3*(2+2+2+2+1);
	int size = node.SCP.n1*node.SCP.n3;
//	myReal temp2P[size13];
//	myReal temp2M[size13];
	int index = 0;
	int stride = 1;
	//Pack at the end only the positive distributions in i2 direction
	//Pack BCC
	//Pack BCCP
	probe13(node.BCCP, temp2P, node.BCCP.iE2, DV_P1_P1_P1, stride, index);
	index += size;
	probe13(node.BCCP, temp2P, node.BCCP.iE2, DV_M1_P1_P1, stride, index);
	index += size;
	//Pack BCCM
	probe13(node.BCCM, temp2P, node.BCCM.iE2, DV_P1_P1_M1, stride, index);
	index += size;
	probe13(node.BCCM, temp2P, node.BCCM.iE2, DV_M1_P1_M1, stride, index);
	index += size;
	//Pack SC
	//Pack SCP
	probe13(node.SCP, temp2P, node.SCP.iE2, DV_ZERO_P1_ZERO, stride, index);
	index += size;
	//Pack FCC
	//Pack FCC12
	probe13(node.FCC12, temp2P, node.FCC12.iE2, DV_P1_P1_ZERO, stride, index);
	index += size;
	probe13(node.FCC12, temp2P, node.FCC12.iE2, DV_M1_P1_ZERO, stride, index);
	index += size;
	//Pack FCC13
	//No need:Zero in i2 direction
	//Pack FCC23
	probe13(node.FCC23, temp2P, node.FCC23.iE2, DV_ZERO_P1_P1, stride, index);
	index += size;
	probe13(node.FCC23, temp2P, node.FCC23.iE2, DV_ZERO_P1_M1, stride, index);
	index += size; 
	
	//Pact at the beginning only the negative distributions in i2 directions
	index = 0;
	//Pack BCC
	//Pack BCCP
	probe13(node.BCCP, temp2M, node.BCCP.iB2, DV_P1_M1_P1, stride, index);
	index += size;
	probe13(node.BCCP, temp2M, node.BCCP.iB2, DV_M1_M1_P1, stride, index);
	index += size;
	//Pack BCCM
	probe13(node.BCCM, temp2M, node.BCCM.iB2, DV_P1_M1_M1, stride, index);
	index += size;
	probe13(node.BCCM, temp2M, node.BCCM.iB2, DV_M1_M1_M1, stride, index);
	index += size;
	//Pack SC
	//Pack SCM
	probe13(node.SCM, temp2M, node.SCM.iB2, DV_ZERO_M1_ZERO, stride, index);
	index += size;
	//Pack FCC
	//Pack FCC12
	probe13(node.FCC12, temp2M, node.FCC12.iB2, DV_P1_M1_ZERO, stride, index);
	index += size;
	probe13(node.FCC12, temp2M, node.FCC12.iB2, DV_M1_M1_ZERO, stride, index);
	index += size;
	//Pack FCC13
	//No need:Zero in i2 direction
	//Pack FCC23
	probe13(node.FCC23, temp2M, node.FCC23.iB2, DV_P1_M1_ZERO, stride, index);
	index += size;
	probe13(node.FCC23, temp2M, node.FCC23.iB2, DV_M1_M1_ZERO, stride, index);
	index += size;
}

//Grid spit in 3-direction: Pack in 12 plane
template<int N> void packToPlane12(lbgrid &node,myReal *temp3P, myReal *temp3M) {
	//int size12 = node.SCP.n1*node.SCP.n2*(2+2+2+2+1);
	int size = node.SCP.n1*node.SCP.n2;
//	myReal temp3P[size12];
//	myReal temp3M[size12];
	int index = 0;
	int stride = 1;
	//Pack at the end only the positive distributions in i3 direction
	//Pack BCC
	//Pack BCCP
	probe12(node.BCCP, temp3P, node.BCCP.iE3, DV_P1_P1_P1, stride, index);
	index += size;
	probe12(node.BCCP, temp3P, node.BCCP.iE3, DV_M1_P1_P1, stride, index);
	index += size;
	probe12(node.BCCP, temp3P, node.BCCP.iE3, DV_P1_M1_P1, stride, index);
	index += size;
	probe12(node.BCCP, temp3P, node.BCCP.iE3, DV_M1_M1_P1, stride, index);
	index += size;
	//Pack SC
	//Pack SCP
	probe12(node.SCP, temp3P, node.SCP.iE3, DV_ZERO_ZERO_P1, stride, index);
	index += size;
	//Pack FCC
	probe12(node.FCC23, temp3P, node.FCC23.iE3, DV_ZERO_P1_P1, stride, index);
	index += size;
	probe12(node.FCC23, temp3P, node.FCC23.iE3, DV_ZERO_M1_P1, stride, index);
	index += size; 
	probe12(node.FCC13, temp3P, node.FCC13.iE3, DV_P1_ZERO_P1, stride, index);
	index += size;
	probe12(node.FCC13, temp3P, node.FCC13.iE3, DV_M1_ZERO_P1, stride, index);
	index += size; 
	
	//Pact at the beginning only the negative distributions in i3 directions
	index = 0;
	//Pack BCC
	probe12(node.BCCM, temp3M, node.BCCM.iE3, DV_M1_M1_M1, stride, index);
	index += size;
	probe12(node.BCCM, temp3M, node.BCCM.iE3, DV_P1_M1_M1, stride, index);
	index += size;
	probe12(node.BCCM, temp3M, node.BCCM.iE3, DV_M1_P1_M1, stride, index);
	index += size;
	probe12(node.BCCM, temp3M, node.BCCM.iE3, DV_P1_P1_M1, stride, index);
	index += size;
	//Pack SC
	//Pack SCM
	probe12(node.SCP, temp3M, node.SCP.iE3, DV_ZERO_ZERO_M1, stride, index);
	index += size;
	//Pack FCC
	probe12(node.FCC23, temp3M, node.FCC23.iE3, DV_ZERO_P1_M1, stride, index);
	index += size;
	probe12(node.FCC23, temp3M, node.FCC23.iE3, DV_ZERO_M1_M1, stride, index);
	index += size; 
	probe12(node.FCC13, temp3M, node.FCC13.iE3, DV_P1_ZERO_M1, stride, index);
	index += size;
	probe12(node.FCC13, temp3M, node.FCC13.iE3, DV_M1_ZERO_M1, stride, index);
	index += size; 
}

//....................................................................................................................................
//PACK FROM PLANES:UNPACK
//USE REVERSE PROBE
template<int N> void unpackFromPlane23(lbgrid &node,myReal *temp1P, myReal *temp1M) {
	int size23 = node.SCP.n2*node.SCP.n3*(2+2+2+2+1);
	int size = node.SCP.n2*node.SCP.n3;
//	myReal temp1P[size23];
//	myReal temp1M[size23];
	int index = 0;
	int stride = 1;
	//Unpack at the end only the positive distributions in i1 direction
	//Unpack BCC
	//Unpack BCCP
	reverseProbe23(node.BCCP, temp1P, node.BCCP.iE1, DV_P1_P1_P1, stride, index);
	index += size;
	reverseProbe23(node.BCCP, temp1P, node.BCCP.iE1, DV_P1_M1_P1, stride, index);
	index += size;
	//Unpack BCCM
	reverseProbe23(node.BCCM, temp1P, node.BCCM.iE1, DV_P1_P1_M1, stride, index);
	index += size;
	reverseProbe23(node.BCCM, temp1P, node.BCCM.iE1, DV_P1_M1_M1, stride, index);
	index += size;
	//Unpack SC
	//Unpack SCP
	reverseProbe23(node.SCP, temp1P, node.SCP.iE1, DV_P1_ZERO_ZERO, stride, index);
	index += size;
	//Unpack FCC
	//Unpack FCC12
	reverseProbe23(node.FCC12, temp1P, node.FCC12.iE1, DV_P1_P1_ZERO, stride, index);
	index += size;
	reverseProbe23(node.FCC12, temp1P, node.FCC12.iE1, DV_P1_M1_ZERO, stride, index);
	index += size;
	//Unpack FCC13
	reverseProbe23(node.FCC13, temp1P, node.FCC13.iE1, DV_P1_ZERO_P1, stride, index);
	index += size;
	reverseProbe23(node.FCC13, temp1P, node.FCC13.iE1, DV_P1_ZERO_M1, stride, index);
	index += size;
	//Unpack FCC23
	//No need: ZERO in i1 direction 
	
	//Pact at the beginning only the negative distributions in i1 directions
	index = 0;
	//Unpack BCC
	//Unpack BCCP
	reverseProbe23(node.BCCP, temp1M, node.BCCP.iB1, DV_M1_P1_P1, stride, index);
	index += size;
	reverseProbe23(node.BCCP, temp1M, node.BCCP.iB1, DV_M1_M1_P1, stride, index);
	index += size;
	//Unpack BCCM
	reverseProbe23(node.BCCM, temp1M, node.BCCM.iB1, DV_M1_P1_M1, stride, index);
	index += size;
	reverseProbe23(node.BCCM, temp1M, node.BCCM.iB1, DV_M1_M1_M1, stride, index);
	index += size;
	//Unpack SC
	//Unpack SCM
	reverseProbe23(node.SCM, temp1M, node.SCM.iB1, DV_M1_ZERO_ZERO, stride, index);
	index += size;
	//Unpack FCC
	//Unpack FCC12
	reverseProbe23(node.FCC12, temp1M, node.FCC12.iB1, DV_M1_P1_ZERO, stride, index);
	index += size;
	reverseProbe23(node.FCC12, temp1M, node.FCC12.iB1, DV_M1_M1_ZERO, stride, index);
	index += size;
	//Unpack FCC13
	reverseProbe23(node.FCC13, temp1M, node.FCC13.iB1, DV_M1_ZERO_P1, stride, index);
	index += size;
	reverseProbe23(node.FCC13, temp1M, node.FCC13.iB1, DV_M1_ZERO_M1, stride, index);
	index += size;
	//Unpack FCC23
	//No need: ZERO in i1 direction 
}



//Grid spit in 2-direction: Unpack in 13 plane
template<int N> void unpackFromPlane13(lbgrid &node,myReal *temp2P, myReal *temp2M) {
	int size13 = node.SCP.n1*node.SCP.n3*(2+2+2+2+1);
	int size = node.SCP.n1*node.SCP.n3;
//	myReal temp2P[size13];
//	myReal temp2M[size13];
	int index = 0;
	int stride = 1;
	//Unpack at the end only the positive distributions in i2 direction
	//Unpack BCC
	//Unpack BCCP
	reverseProbe13(node.BCCP, temp2P, node.BCCP.iE2, DV_P1_P1_P1, stride, index);
	index += size;
	reverseProbe13(node.BCCP, temp2P, node.BCCP.iE2, DV_M1_P1_P1, stride, index);
	index += size;
	//Unpack BCCM
	reverseProbe13(node.BCCM, temp2P, node.BCCM.iE2, DV_P1_P1_M1, stride, index);
	index += size;
	reverseProbe13(node.BCCM, temp2P, node.BCCM.iE2, DV_M1_P1_M1, stride, index);
	index += size;
	//Unpack SC
	//Unpack SCP
	reverseProbe13(node.SCP, temp2P, node.SCP.iE2, DV_ZERO_P1_ZERO, stride, index);
	index += size;
	//Unpack FCC
	//Unpack FCC12
	reverseProbe13(node.FCC12, temp2P, node.FCC12.iE2, DV_P1_P1_ZERO, stride, index);
	index += size;
	reverseProbe13(node.FCC12, temp2P, node.FCC12.iE2, DV_M1_P1_ZERO, stride, index);
	index += size;
	//Unpack FCC13
	//No need:Zero in i2 direction
	//Unpack FCC23
	reverseProbe13(node.FCC23, temp2P, node.FCC23.iE2, DV_ZERO_P1_P1, stride, index);
	index += size;
	reverseProbe13(node.FCC23, temp2P, node.FCC23.iE2, DV_ZERO_P1_M1, stride, index);
	index += size; 
	
	//Pact at the beginning only the negative distributions in i2 directions
	index = 0;
	//Unpack BCC
	//Unpack BCCP
	reverseProbe13(node.BCCP, temp2M, node.BCCP.iB2, DV_P1_M1_P1, stride, index);
	index += size;
	reverseProbe13(node.BCCP, temp2M, node.BCCP.iB2, DV_M1_M1_P1, stride, index);
	index += size;
	//Unpack BCCM
	reverseProbe13(node.BCCM, temp2M, node.BCCM.iB2, DV_P1_M1_M1, stride, index);
	index += size;
	reverseProbe13(node.BCCM, temp2M, node.BCCM.iB2, DV_M1_M1_M1, stride, index);
	index += size;
	//Unpack SC
	//Unpack SCM
	reverseProbe13(node.SCM, temp2M, node.SCM.iB2, DV_ZERO_M1_ZERO, stride, index);
	index += size;
	//Unpack FCC
	//Unpack FCC12
	reverseProbe13(node.FCC12, temp2M, node.FCC12.iB2, DV_P1_M1_ZERO, stride, index);
	index += size;
	reverseProbe13(node.FCC12, temp2M, node.FCC12.iB2, DV_M1_M1_ZERO, stride, index);
	index += size;
	//Unpack FCC13
	//No need:Zero in i2 direction
	//Unpack FCC23
	reverseProbe13(node.FCC23, temp2M, node.FCC23.iB2, DV_P1_M1_ZERO, stride, index);
	index += size;
	reverseProbe13(node.FCC23, temp2M, node.FCC23.iB2, DV_M1_M1_ZERO, stride, index);
	index += size;
}

//Grid spit in 3-direction: Unpack in 12 plane
template<int N> void unpackFromPlane12(lbgrid &node, myReal *temp3P, myReal *temp3M) {
//	int size12 = node.SCP.n1*node.SCP.n2*(2+2+2+2+1);
	int size = node.SCP.n1*node.SCP.n2;
//	myReal temp3P[size12];
//	myReal temp3M[size12];
	int index = 0;
	int stride = 1;
	//Unpack at the end only the positive distributions in i3 direction
	//Unpack BCC
	//Unpack BCCP
	reverseProbe12(node.BCCP, temp3P, node.BCCP.iE3, DV_P1_P1_P1, stride, index);
	index += size;
	reverseProbe12(node.BCCP, temp3P, node.BCCP.iE3, DV_M1_P1_P1, stride, index);
	index += size;
	reverseProbe12(node.BCCP, temp3P, node.BCCP.iE3, DV_P1_M1_P1, stride, index);
	index += size;
	reverseProbe12(node.BCCP, temp3P, node.BCCP.iE3, DV_M1_M1_P1, stride, index);
	index += size;
	//Unpack SC
	//Unpack SCP
	reverseProbe12(node.SCP, temp3P, node.SCP.iE3, DV_ZERO_ZERO_P1, stride, index);
	index += size;
	//Unpack FCC
	reverseProbe12(node.FCC23, temp3P, node.FCC23.iE3, DV_ZERO_P1_P1, stride, index);
	index += size;
	reverseProbe12(node.FCC23, temp3P, node.FCC23.iE3, DV_ZERO_M1_P1, stride, index);
	index += size; 
	reverseProbe12(node.FCC13, temp3P, node.FCC13.iE3, DV_P1_ZERO_P1, stride, index);
	index += size;
	reverseProbe12(node.FCC13, temp3P, node.FCC13.iE3, DV_M1_ZERO_P1, stride, index);
	index += size; 
	
	//Pact at the beginning only the negative distributions in i3 directions
	index = 0;
	//Unpack BCC
	reverseProbe12(node.BCCM, temp3M, node.BCCM.iE3, DV_M1_M1_M1, stride, index);
	index += size;
	reverseProbe12(node.BCCM, temp3M, node.BCCM.iE3, DV_P1_M1_M1, stride, index);
	index += size;
	reverseProbe12(node.BCCM, temp3M, node.BCCM.iE3, DV_M1_P1_M1, stride, index);
	index += size;
	reverseProbe12(node.BCCM, temp3M, node.BCCM.iE3, DV_P1_P1_M1, stride, index);
	index += size;
	//Unpack SC
	//Unpack SCM
	reverseProbe12(node.SCP, temp3M, node.SCP.iE3, DV_ZERO_ZERO_M1, stride, index);
	index += size;
	//Unpack FCC
	reverseProbe12(node.FCC23, temp3M, node.FCC23.iE3, DV_ZERO_P1_M1, stride, index);
	index += size;
	reverseProbe12(node.FCC23, temp3M, node.FCC23.iE3, DV_ZERO_M1_M1, stride, index);
	index += size; 
	reverseProbe12(node.FCC13, temp3M, node.FCC13.iE3, DV_P1_ZERO_M1, stride, index);
	index += size;
	reverseProbe12(node.FCC13, temp3M, node.FCC13.iE3, DV_M1_ZERO_M1, stride, index);
	index += size; 
}
	
#endif	
	
