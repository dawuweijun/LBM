#ifndef _GRIDV_PROBE_H_
#define  _GRIDV_PROBE_H_
#include<iostream>

#include "gridV.h"

//Probe: Get data at any plane i and put it in a temporary 1D array
template<int N>
void probe12(gridV<N> &node, myReal *tempData, int i3, int k, int stride, int index) { 
	
	for (int i1 = node.iB1-1; i1 <= node.iE1+1; i1+=stride) {
		for (int i2 = node.iB2-1; i2 <= node.iE2+1; i2+=stride) {
			tempData[index] = node(i1,i2,i3,k);
			index++;
		}
	}
}

template<int N>
void probe13(gridV<N> &node, myReal *tempData, int i2, int k, int stride, int index) {
	
	for (int i1 = node.iB1-1; i1 <= node.iE1+1; i1+=stride) {
		for (int i3 = node.iB3-1; i3 <= node.iE3+1; i3+=stride) {
			tempData[index] = node(i1,i2,i3,k);
			index++;
		}
	}
}

template<int N>
void probe23(gridV<N> &node, myReal *tempData, int i1, int k, int stride, int index) {
	
	for (int i2 = node.iB2-1; i2 <= node.iE2+1; i2+=stride) {
		for (int i3 = node.iB3-1; i3 <= node.iE3+1; i3+=stride) {
			tempData[index] = node(i1,i2,i3,k);
			index++;
		}
	}
}
//Probe over

//Get my plane from a temporary array - Reverse Probe
template<int N>
void reverseProbe12(gridV<N> &node, myReal *tempData, int i3, int k, int stride, int index) {
	for (int i1 = node.iB1-1; i1 <= node.iE1+1; i1+=stride) {
		for (int i2 = node.iB2-1; i2 <= node.iE2+1; i2+=stride) {
			node(i1,i2,i3,k) = tempData[index];
			index++;
		}
	}
}

template<int N>
void reverseProbe13(gridV<N> &node, myReal *tempData, int i2, int k, int stride, int index) {
	
	for (int i1 = node.iB1-1; i1 <= node.iE1+1; i1+=stride) {
		for (int i3 = node.iB3-1; i3 <= node.iE3+1; i3+=stride) {
			node(i1,i2,i3,k) = tempData[index];
			index++;
			
		}
	}
}

template<int N>
void reverseProbe23(gridV<N> &node, myReal *tempData, int i1, int k, int stride, int index) {
	
	for (int i2 = node.iB2-1; i2 <= node.iE2+1; i2+=stride) {
		for (int i3 = node.iB3-1; i3 <= node.iE3+1; i3+=stride) {
			node(i1,i2,i3,k) = tempData[index];
			index++;
		}
	}
}

//Reverse probe over
#endif
