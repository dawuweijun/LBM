#ifndef _GRIDV_PACK_H_
#define  _GRIDV_PACK_H_
#include<iostream>
#include<stdlib.h>

#include "gridV.h"
#include "gridV_Probe.h"
#include "globals3D.h"

//Grid spit in 1-direction: Pack in 23 plane
//Allocate size
myReal* allocateForPacking(int m1, int m2, int m3) {
	myReal *temp;
	//sort in descending order
	if (m1<m2) swap(m1,m2);
	if (m2<m3) swap(m2,m3);
	if (m1<m2) swap(m1,m2);
	// temp = (double*) malloc(m1*m2*18);
	temp = new double[m1*m2*18];
	return temp;
}
	
	
void packTo1P(lbgrid &node,myReal *temp1P,int &index) {
//	int size23 = node.SCP.m2*node.SCP.m3*(2+2+2+2+1);
	int size = node.SCP.m2*node.SCP.m3;
//	myReal temp1P[size23];
//	myReal temp1M[size23];
	int stride = 1;
	//Pack at the end only the positive distributions in i1 direction
	//Pack BCC
	//Pack BCCP
	probe23(node.BCCP, temp1P, node.BCCP.iE1, DV_P1_P1_P1, stride, index);
	index += size;
	probe23(node.BCCP, temp1P, node.BCCP.iE1, DV_P1_M1_P1, stride, index);
	index += size;
	//////Pack BCCM
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
	//////Pack FCC13
	probe23(node.FCC13, temp1P, node.FCC13.iE1, DV_P1_ZERO_P1, stride, index);
	index += size;
	probe23(node.FCC13, temp1P, node.FCC13.iE1, DV_P1_ZERO_M1, stride, index);
	index += size;
	//Pack FCC23
	//No need: ZERO in i1 direction 
}
void packTo1M(lbgrid &node,myReal *temp1M, int &index) {	
	//Pact at the beginning only the negative distributions in i1 directions
	//Pack BCC
	//Pack BCCP
	int size = node.SCP.m2*node.SCP.m3;
	int stride = 1;
	probe23(node.BCCP, temp1M, node.BCCP.iB1, DV_M1_P1_P1, stride, index);
	index += size;
	probe23(node.BCCP, temp1M, node.BCCP.iB1, DV_M1_M1_P1, stride, index);
	index += size;
	///////////////Pack BCCM
	probe23(node.BCCM, temp1M, node.BCCM.iB1, DV_M1_P1_M1, stride, index);
	index += size;
	probe23(node.BCCM, temp1M, node.BCCM.iB1, DV_M1_M1_M1, stride, index);
	index += size;
	//Pack SC
	//Pack SCM
	probe23(node.SCM, temp1M, node.SCM.iB1, DV_M1_ZERO_ZERO, stride, index);
	index += size;
	// Pack FCC
	//Pack FCC12
	probe23(node.FCC12, temp1M, node.FCC12.iB1, DV_M1_P1_ZERO, stride, index);
	index += size;
	probe23(node.FCC12, temp1M, node.FCC12.iB1, DV_M1_M1_ZERO, stride, index);
	index += size;
	//////////////Pack FCC13
	probe23(node.FCC13, temp1M, node.FCC13.iB1, DV_M1_ZERO_P1, stride, index);
	index += size;
	probe23(node.FCC13, temp1M, node.FCC13.iB1, DV_M1_ZERO_M1, stride, index);
	index += size;
	//Pack FCC23
	//No need: ZERO in i1 direction 
}

void packMSgTo1Neb(lbgrid &node, lbgrid &cell, myReal *temp1Psend, myReal *temp1Msend) {
		int index = 0;
		packTo1P(node,temp1Psend,index);
		packTo1P(cell,temp1Psend,index);
		index = 0;
		packTo1M(node,temp1Msend,index);
		packTo1M(cell,temp1Msend,index);
}

//Grid spit in 2-direction: Pack in 13 plane
void packTo2P(lbgrid &node, myReal *temp2P, int &index) {
	//int size13 = node.SCP.m1*node.SCP.m3*(2+2+2+2+1);
	int size = node.SCP.m1*node.SCP.m3;
//	myReal temp2P[size13];
//	myReal temp2M[size13];
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
}
void packTo2M(lbgrid &node, myReal *temp2M, int &index) {
	//Pact at the beginning only the negative distributions in i2 directions
	int size = node.SCP.m1*node.SCP.m3;
	int stride = 1;
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
	probe13(node.FCC23, temp2M, node.FCC23.iB2, DV_ZERO_M1_P1, stride, index);
	index += size;
	probe13(node.FCC23, temp2M, node.FCC23.iB2, DV_ZERO_M1_M1, stride, index);
	index += size;
}

void packMSgTo2Neb(lbgrid &node, lbgrid &cell, myReal *temp2Psend, myReal *temp2Msend) {
		int index = 0;
		packTo2P(node,temp2Psend,index);
		packTo2P(cell,temp2Psend,index);
		index = 0;
		packTo2M(node,temp2Msend,index);
		packTo2M(cell,temp2Msend,index);
}

//Grid spit in 3-direction: Pack in 12 plane
void packTo3P(lbgrid &node,myReal *temp3P, int &index) {
	//int size12 = node.SCP.m1*node.SCP.m2*(2+2+2+2+1);
	int size = node.SCP.m1*node.SCP.m2;
//	myReal temp3P[size12];
//	myReal temp3M[size12];

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
}
void packTo3M(lbgrid &node,myReal *temp3M, int &index) {	
	//Pact at the beginning only the negative distributions in i3 directions
	int size = node.SCP.m1*node.SCP.m2;
	int stride = 1;
	//Pack BCC
	probe12(node.BCCM, temp3M, node.BCCM.iB3, DV_M1_M1_M1, stride, index);
	index += size;
	probe12(node.BCCM, temp3M, node.BCCM.iB3, DV_P1_M1_M1, stride, index);
	index += size;
	probe12(node.BCCM, temp3M, node.BCCM.iB3, DV_M1_P1_M1, stride, index);
	index += size;
	probe12(node.BCCM, temp3M, node.BCCM.iB3, DV_P1_P1_M1, stride, index);
	index += size;
	//Pack SC
	//Pack SCM
	probe12(node.SCM, temp3M, node.SCP.iB3, DV_ZERO_ZERO_M1, stride, index);
	index += size;
	//Pack FCC
	probe12(node.FCC23, temp3M, node.FCC23.iB3, DV_ZERO_P1_M1, stride, index);
	index += size;
	probe12(node.FCC23, temp3M, node.FCC23.iB3, DV_ZERO_M1_M1, stride, index);
	index += size; 
	probe12(node.FCC13, temp3M, node.FCC13.iB3, DV_P1_ZERO_M1, stride, index);
	index += size;
	probe12(node.FCC13, temp3M, node.FCC13.iB3, DV_M1_ZERO_M1, stride, index);
	index += size; 
}

void packMSgTo3Neb(lbgrid &node, lbgrid &cell, myReal *temp3Psend, myReal *temp3Msend) {
		int index = 0;
		packTo3P(node,temp3Psend,index);
		packTo3P(cell,temp3Psend,index);
		index = 0;
		packTo3M(node,temp3Msend,index);
		packTo3M(cell,temp3Msend,index);
}


//....................................................................................................................................
//PACK FROM PLANES:UNPACK
//USE REVERSE PROBE
void unpackFrom1P(lbgrid &node,myReal *temp1P, int &index) {
	//int size23 = node.SCP.m2*node.SCP.m3*(2+2+2+2+1);
	int size = node.SCP.m2*node.SCP.m3;
//	myReal temp1P[size23];
//	myReal temp1M[size23];

	int stride = 1;
	//Unpack at the end only the positive distributions in i1 direction
	//Unpack BCC
	//Unpack BCCP
	reverseProbe23(node.BCCP, temp1P, node.BCCP.iB1-1, DV_P1_P1_P1, stride, index);
	index += size;
	reverseProbe23(node.BCCP, temp1P, node.BCCP.iB1-1, DV_P1_M1_P1, stride, index);
	index += size;
	/////////Unpack BCCM
	reverseProbe23(node.BCCM, temp1P, node.BCCM.iB1-1, DV_P1_P1_M1, stride, index);
	index += size;
	reverseProbe23(node.BCCM, temp1P, node.BCCM.iB1-1, DV_P1_M1_M1, stride, index);
	index += size;
	//Unpack SC
	//Unpack SCP
	reverseProbe23(node.SCP, temp1P, node.SCP.iB1-1, DV_P1_ZERO_ZERO, stride, index);
	index += size;
	//Unpack FCC
	//Unpack FCC12
	reverseProbe23(node.FCC12, temp1P, node.FCC12.iB1-1, DV_P1_P1_ZERO, stride, index);
	index += size;
	reverseProbe23(node.FCC12, temp1P, node.FCC12.iB1-1, DV_P1_M1_ZERO, stride, index);
	index += size;
	//////////Unpack FCC13
	reverseProbe23(node.FCC13, temp1P, node.FCC13.iB1-1, DV_P1_ZERO_P1, stride, index);
	index += size;
	reverseProbe23(node.FCC13, temp1P, node.FCC13.iB1-1, DV_P1_ZERO_M1, stride, index);
	index += size;
	//Unpack FCC23
	//No need: ZERO in i1 direction 
}
void unpackFrom1M(lbgrid &node,myReal *temp1M, int &index) {
	//Pact at the beginning only the negative distributions in i1 directions
	int size = node.SCP.m2*node.SCP.m3;
	int stride = 1;
	//Unpack BCC
	//Unpack BCCP
	reverseProbe23(node.BCCP, temp1M, node.BCCP.iE1+1, DV_M1_P1_P1, stride, index);
	index += size;
	reverseProbe23(node.BCCP, temp1M, node.BCCP.iE1+1, DV_M1_M1_P1, stride, index);
	index += size;
	/////Unpack BCCM
	reverseProbe23(node.BCCM, temp1M, node.BCCM.iE1+1, DV_M1_P1_M1, stride, index);
	index += size;
	reverseProbe23(node.BCCM, temp1M, node.BCCM.iE1+1, DV_M1_M1_M1, stride, index);
	index += size;
	//Unpack SC
	//Unpack SCM
	reverseProbe23(node.SCM, temp1M, node.SCM.iE1+1, DV_M1_ZERO_ZERO, stride, index);
	index += size;
	//Unpack FCC
	//Unpack FCC12
	reverseProbe23(node.FCC12, temp1M, node.FCC12.iE1+1, DV_M1_P1_ZERO, stride, index);
	index += size;
	reverseProbe23(node.FCC12, temp1M, node.FCC12.iE1+1, DV_M1_M1_ZERO, stride, index);
	index += size;
	//////Unpack FCC13
	reverseProbe23(node.FCC13, temp1M, node.FCC13.iE1+1, DV_M1_ZERO_P1, stride, index);
	index += size;
	reverseProbe23(node.FCC13, temp1M, node.FCC13.iE1+1, DV_M1_ZERO_M1, stride, index);
	index += size;
	//Unpack FCC23
	//No need: ZERO in i1 direction 
}

void receiveMSgFrom1Neb(lbgrid &node, lbgrid &cell, myReal *temp1Prec, myReal *temp1Mrec) {
	int index = 0;
	unpackFrom1P(node,temp1Prec,index);
	unpackFrom1P(cell,temp1Prec,index);
	index = 0;
	unpackFrom1M(node,temp1Mrec,index);
	unpackFrom1M(cell,temp1Mrec,index);
}

//Grid spit in 2-direction: Unpack in 13 plane
void unpackFrom2P(lbgrid &node,myReal *temp2P, int &index) {
	// int size13 = node.SCP.m1*node.SCP.m3*(2+2+2+2+1);
	int size = node.SCP.m1*node.SCP.m3;
//	myReal temp2P[size13];
//	myReal temp2M[size13];

	int stride = 1;
	//Unpack at the end only the positive distributions in i2 direction
	//Unpack BCC
	//Unpack BCCP
	reverseProbe13(node.BCCP, temp2P, node.BCCP.iB2-1, DV_P1_P1_P1, stride, index);
	index += size;
	reverseProbe13(node.BCCP, temp2P, node.BCCP.iB2-1, DV_M1_P1_P1, stride, index);
	index += size;
	//Unpack BCCM
	reverseProbe13(node.BCCM, temp2P, node.BCCM.iB2-1, DV_P1_P1_M1, stride, index);
	index += size;
	reverseProbe13(node.BCCM, temp2P, node.BCCM.iB2-1, DV_M1_P1_M1, stride, index);
	index += size;
	//Unpack SC
	//Unpack SCP
	reverseProbe13(node.SCP, temp2P, node.SCP.iB2-1, DV_ZERO_P1_ZERO, stride, index);
	index += size;
	//Unpack FCC
	//Unpack FCC12
	reverseProbe13(node.FCC12, temp2P, node.FCC12.iB2-1, DV_P1_P1_ZERO, stride, index);
	index += size;
	reverseProbe13(node.FCC12, temp2P, node.FCC12.iB2-1, DV_M1_P1_ZERO, stride, index);
	index += size;
	//Unpack FCC13
	//No need:Zero in i2 direction
	//Unpack FCC23
	reverseProbe13(node.FCC23, temp2P, node.FCC23.iB2-1, DV_ZERO_P1_P1, stride, index);
	index += size;
	reverseProbe13(node.FCC23, temp2P, node.FCC23.iB2-1, DV_ZERO_P1_M1, stride, index);
	index += size; 
}
void unpackFrom2M(lbgrid &node,myReal *temp2M, int &index) {
	//Pact at the beginning only the negative distributions in i2 directions
	int size = node.SCP.m1*node.SCP.m3;
	int stride = 1;
	//Unpack BCC
	//Unpack BCCP
	reverseProbe13(node.BCCP, temp2M, node.BCCP.iE2+1, DV_P1_M1_P1, stride, index);
	index += size;
	reverseProbe13(node.BCCP, temp2M, node.BCCP.iE2+1, DV_M1_M1_P1, stride, index);
	index += size;
	//Unpack BCCM
	reverseProbe13(node.BCCM, temp2M, node.BCCM.iE2+1, DV_P1_M1_M1, stride, index);
	index += size;
	reverseProbe13(node.BCCM, temp2M, node.BCCM.iE2+1, DV_M1_M1_M1, stride, index);
	index += size;
	//Unpack SC
	//Unpack SCM
	reverseProbe13(node.SCM, temp2M, node.SCM.iE2+1, DV_ZERO_M1_ZERO, stride, index);
	index += size;
	//Unpack FCC
	//Unpack FCC12
	reverseProbe13(node.FCC12, temp2M, node.FCC12.iE2+1, DV_P1_M1_ZERO, stride, index);
	index += size;
	reverseProbe13(node.FCC12, temp2M, node.FCC12.iE2+1, DV_M1_M1_ZERO, stride, index);
	index += size;
	//Unpack FCC13
	//No need:Zero in i2 direction
	//Unpack FCC23
	reverseProbe13(node.FCC23, temp2M, node.FCC23.iE2+1, DV_ZERO_M1_P1, stride, index);
	index += size;
	reverseProbe13(node.FCC23, temp2M, node.FCC23.iE2+1, DV_ZERO_M1_M1, stride, index);
	index += size;
}

void receiveMSgFrom2Neb(lbgrid &node, lbgrid &cell, myReal *temp2Prec, myReal *temp2Mrec) {
	int index = 0;
	unpackFrom2P(node,temp2Prec,index);
	unpackFrom2P(cell,temp2Prec,index);
	index = 0;
	unpackFrom2M(node,temp2Mrec,index);
	unpackFrom2M(cell,temp2Mrec,index);
}

//Grid spit in 3-direction: Unpack in 12 plane
void unpackFrom3P(lbgrid &node, myReal *temp3P, int &index) {
//	int size12 = node.SCP.m1*node.SCP.m2*(2+2+2+2+1);
	int size = node.SCP.m1*node.SCP.m2;
//	myReal temp3P[size12];
//	myReal temp3M[size12];

	int stride = 1;
	//Unpack at the end only the positive distributions in i3 direction
	//Unpack BCC
	//Unpack BCCP
	reverseProbe12(node.BCCP, temp3P, node.BCCP.iB3-1, DV_P1_P1_P1, stride, index);
	index += size;
	reverseProbe12(node.BCCP, temp3P, node.BCCP.iB3-1, DV_M1_P1_P1, stride, index);
	index += size;
	reverseProbe12(node.BCCP, temp3P, node.BCCP.iB3-1, DV_P1_M1_P1, stride, index);
	index += size;
	reverseProbe12(node.BCCP, temp3P, node.BCCP.iB3-1, DV_M1_M1_P1, stride, index);
	index += size;
	//Unpack SC
	//Unpack SCP
	reverseProbe12(node.SCP, temp3P, node.SCP.iB3-1, DV_ZERO_ZERO_P1, stride, index);
	index += size;
	//Unpack FCC
	reverseProbe12(node.FCC23, temp3P, node.FCC23.iB3-1, DV_ZERO_P1_P1, stride, index);
	index += size;
	reverseProbe12(node.FCC23, temp3P, node.FCC23.iB3-1, DV_ZERO_M1_P1, stride, index);
	index += size; 
	reverseProbe12(node.FCC13, temp3P, node.FCC13.iB3-1, DV_P1_ZERO_P1, stride, index);
	index += size;
	reverseProbe12(node.FCC13, temp3P, node.FCC13.iB3-1, DV_M1_ZERO_P1, stride, index);
	index += size; 
}
void unpackFrom3M(lbgrid &node, myReal *temp3M, int &index) {
	//Pact at the beginning only the negative distributions in i3 directions
	int size = node.SCP.m1*node.SCP.m2;
	int stride = 1;
	//Unpack BCC
	
	reverseProbe12(node.BCCM, temp3M, node.BCCM.iE3+1, DV_M1_M1_M1, stride, index);
	index += size;
	reverseProbe12(node.BCCM, temp3M, node.BCCM.iE3+1, DV_P1_M1_M1, stride, index);
	index += size;
	reverseProbe12(node.BCCM, temp3M, node.BCCM.iE3+1, DV_M1_P1_M1, stride, index);
	index += size;
	reverseProbe12(node.BCCM, temp3M, node.BCCM.iE3+1, DV_P1_P1_M1, stride, index);
	index += size;
	
	//Unpack SC
	//Unpack SCM
	reverseProbe12(node.SCM, temp3M, node.SCP.iE3+1, DV_ZERO_ZERO_M1, stride, index);
	
	index += size;
	//Unpack FCC
	reverseProbe12(node.FCC23, temp3M, node.FCC23.iE3+1, DV_ZERO_P1_M1, stride, index);
	index += size;
	
	reverseProbe12(node.FCC23, temp3M, node.FCC23.iE3+1, DV_ZERO_M1_M1, stride, index);
	index += size; 
	reverseProbe12(node.FCC13, temp3M, node.FCC13.iE3+1, DV_P1_ZERO_M1, stride, index);
	index += size;
	reverseProbe12(node.FCC13, temp3M, node.FCC13.iE3+1, DV_M1_ZERO_M1, stride, index);
	index += size;
}

void receiveMSgFrom3Neb(lbgrid &node, lbgrid &cell, myReal *temp3Prec, myReal *temp3Mrec) {
		
	int index = 0;
	unpackFrom3P(node,temp3Prec,index);
	unpackFrom3P(cell,temp3Prec,index);
	index = 0;
	unpackFrom3M(node,temp3Mrec,index);
	unpackFrom3M(cell,temp3Mrec,index);
}
	
#endif	
	
