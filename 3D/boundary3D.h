#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "gridV.h"
#include "globals3D.h"
#include "moments3D.h"

//Class to hold node type and cell type
class gridtype {
	public:
		gridtype(int _n1, int _n2, int _n3): n1(_n1), n2(_n2), n3(_n3), size(n1*n2*n3) {
			data = new myReal[size];
		}
		int n1,n2,n3,size;
		myReal *data;
		//Overload the () operator
		myReal operator() (int i, int j, int k) const {
			return data[k*n1*n2 + j*n1 + i];
		}
		myReal& operator() (int i, int j, int k) {
			return data[k*n1*n2 + j*n1 + i];
		}
};

void setTypeSphere(gridtype nodeType, gridtype cellType) {
	for (int i3 = 0; i3 <= nodeType.n3; i3++) {
		for (int i2 = 0; i2 <= nodeType.n2; i2++) {
			for (int i1 = 0; i1 <= nodeType.n1; i1++) {
				if (pow(i1,2) + pow(i2,2) + pow(i3,2) <= pow(radius,2)) 
					nodeType(i1,i2,i3) = FLUID;
				else
					nodeType(i1,i2,i3) = SOLID;
			}
		}
	}
}

int getFCCReverse(int k) {
	int k_opp;
	switch(k) {
		case 0: k_opp = 3; break;
		case 1: k_opp = 2; break;
		case 2: k_opp = 1; break;
		case 3: k_opp = 0; break;
		default: cout<<"Looping is wrong in FCC Boundary settings"<<endl;
	}
	return k_opp;
}

void setSCAndFCC(lbgrid &node, gridtype nodeType, int i1, int i2, int i3) {
	int k_local, k_local_opp;
	//SC swap
	//SCP:1-3
	for (int k = 1; k <= 3; k++)
		if (nodeType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) 
			swap(node.SCM(i1,i2,i3,k-1), node.SCP(i1+c_x[k],i2+c_y[k],i3+c_z[k],k-1));
	//SCM:4-6
	for (int k = 4; k <= 6; k++)
		if (nodeType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) 
			swap(node.SCP(i1,i2,i3,k-4), node.SCM(i1+c_x[k],i2+c_y[k],i3+c_z[k],k-4));
	//FCC
	//FCC 12:7-10
	for (int k = 7; k <=10; k++) {
		k_local = k-7;
		if (nodeType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) { 
			k_local_opp = getFCCReverse(k_local);
			swap(node.FCC12(i1,i2,i3,k_local_opp), node.FCC12(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local));
		}
	}
	//FCC 13:11-14
	for (int k = 11; k <=14; k++) {
		k_local = k-11;
		if (nodeType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) { 
			k_local_opp = getFCCReverse(k_local);
			swap(node.FCC13(i1,i2,i3,k_local_opp), node.FCC13(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local));
		}
	}
	//FCC 23:15-19
	for (int k = 15; k <=19; k++) {
		k_local = k-15;
		if (nodeType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) { 
			k_local_opp = getFCCReverse(k_local);
			swap(node.FCC23(i1,i2,i3,k_local_opp), node.FCC23(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local));
		}
	}
	
}

void setBCC(lbgrid &node, lbgrid &cell, gridtype nodeType, gridtype cellType, int i1, int i2, int i3) {
	//BCCP:20-23
	for (int k = 20; k <=23; k++) 
		if (nodeType(i1,i2,i3) == FLUID && cellType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) 
			swap(node.BCCM(i1,i2,i3,k-20), cell.BCCP(i1+c_x[k],i2+c_y[k],i3+c_z[k],k-20));	
	//BCCM:24-27
	for (int k = 24; k <=27; k++) 
		if (nodeType(i1,i2,i3) == FLUID && cellType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) 
			swap(node.BCCP(i1,i2,i3,k-24), cell.BCCM(i1+c_x[k],i2+c_y[k],i3+c_z[k],k-24));
}

void setInternalBoundary(lbgrid &node, lbgrid &cell, gridtype nodeType, gridtype cellType) {
	int n1, n2, n3; 
	n1 = node.SCP.n1;
	n2 = node.SCP.n2;
	n3 = node.SCP.n3;
	
	for (int i3 = 0; i3 <= n3; i3++) {
		for (int i2 = 0; i2 <= n2; i2++) {
			for (int i1 = 0; i1 <= n1; i1++) {
				//SC and FCC: bounce back between similar grid points: No node to cell advection
				setSCAndFCC(node,nodeType,i1,i2,i3);  
				setSCAndFCC(cell,cellType,i1,i2,i3);
				//BCC boundary: bounce back between node and cell
				setBCC(node,cell,nodeType,cellType,i1,i2,i3); //Boundary btw node and cell
				setBCC(cell,node,cellType,nodeType,i1,i2,i3); //Boundary btw cell and node
			}
		}
	}			
}

#endif
