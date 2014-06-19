#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include<math.h>
#include<cstdlib>
#include "gridV.h"
#include "globals3D.h"
#include "moments3D.h"

//Class to hold node type and cell type
class gridtype {
	public:
		gridtype(int _n1, int _n2,int _n3, int _numPad1, int _numPad2, int _numPad3) {
			m1= _n1+2*_numPad1;
			m2 = _n2+2*_numPad2;
			m3 = _n3+2*_numPad3;
			size = m1*m2*m3;
			data = new myReal[size];
		}
		int n1,n2,n3,m1,m2,m3,size;
		myReal *data;
		//Overload the () operator
		myReal operator() (int i, int j, int k) const {
			return data[k*m1*m2 + j*m1 + i];
		}
		myReal& operator() (int i, int j, int k) {
			return data[k*m1*m2 + j*m1 + i];
		}
};

void setTypeSphere(gridtype nodeType, gridtype cellType) {
	for (int i3 = 0; i3 < nodeType.m3; i3++) {
		for (int i2 = 0; i2 < nodeType.m2; i2++) {
			for (int i1 = 0; i1 < nodeType.m1; i1++) {
				if (pow(i1-center_x,2) + pow(i2-center_y,2) + pow(i3-center_z,2) <= pow(radius,2)) 
					nodeType(i1,i2,i3) = FLUID;
				else
					nodeType(i1,i2,i3) = SOLID;
			}
		}
	}
	for (int i3 = 0; i3 < cellType.m3; i3++) {
		for (int i2 = 0; i2 < cellType.m2; i2++) {
			for (int i1 = 0; i1 < cellType.m1; i1++) {
				if (pow(i1+0.5-center_x,2) + pow(i2+0.5-center_y,2) + pow(i3+0.5-center_z,2) <= pow(radius,2))  
					cellType(i1,i2,i3) = FLUID;
				else
					cellType(i1,i2,i3) = SOLID;
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
	for (int k = 15; k <=18; k++) {
		k_local = k-15;
		if (nodeType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) { 
			k_local_opp = getFCCReverse(k_local);
			swap(node.FCC23(i1,i2,i3,k_local_opp), node.FCC23(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local));
		}
	}
	
}

void setBCC(lbgrid &node, lbgrid &cell, gridtype nodeType, gridtype cellType, int i1, int i2, int i3) {
	//Btw node and cell
	//BCCP:20-23
	int cx, cy, cz;
	for (int k = 19; k <=22; k++) {
		cx = (c_x[k]>0?0:2*c_x[k]); //Replica grid is displaced by 0.5 units in each direction
		cy = (c_y[k]>0?0:2*c_y[k]);
		cz = (c_z[k]>0?0:2*c_z[k]);
		if (nodeType(i1,i2,i3) == FLUID && cellType(i1-cx,i2-cy,i3-cz) == SOLID) {
			cout<<"Sphere encountered"<<endl;
			swap(node.BCCM(i1,i2,i3,k-19), cell.BCCP(i1-cx,i2-cy,i3-cz,k-19));
		}	
	}
	//BCCM:24-27
	for (int k = 23; k <=26; k++) {
		cx = (c_x[k]>0?0:2*c_x[k]); 
		cy = (c_y[k]>0?0:2*c_y[k]);
		cz = (c_z[k]>0?0:2*c_z[k]);
		if (nodeType(i1,i2,i3) == FLUID && cellType(i1+cx,i2+cy,i3+cz) == SOLID) 
			swap(node.BCCP(i1,i2,i3,k-23), cell.BCCM(i1+cx,i2+cy,i3+cz,k-23));
	}
	//Btw cell and node
	for (int k = 19; k <=22; k++) {
		cx = (c_x[k]<0?0:2*c_x[k]); //Replica grid is displaced by 0.5 units in each direction
		cy = (c_y[k]<0?0:2*c_y[k]);
		cz = (c_z[k]<0?0:2*c_z[k]);
		if (cellType(i1,i2,i3) == FLUID && nodeType(i1-cx,i2-cy,i3-cz) == SOLID) 
			swap(cell.BCCM(i1,i2,i3,k-19), node.BCCP(i1-cx,i2-cy,i3-cz,k-19));	
	}
	//BCCM:24-27
	for (int k = 23; k <=26; k++) {
		cx = (c_x[k]<0?0:2*c_x[k]); 
		cy = (c_y[k]<0?0:2*c_y[k]);
		cz = (c_z[k]<0?0:2*c_z[k]);
		if (cellType(i1,i2,i3) == FLUID && nodeType(i1+cx,i2+cy,i3+cz) == SOLID) 
			swap(cell.BCCP(i1,i2,i3,k-23), node.BCCM(i1+cx,i2+cy,i3+cz,k-23));
	}
}

void setInternalBoundary(lbgrid &node, lbgrid &cell, gridtype nodeType, gridtype cellType) {
	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
		for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
			for (int i1 = node.SCP.iB1; i1 <= node.SCP.iE1; i1++) {
				//SC and FCC: bounce back between similar grid points: No node to cell advection
				setSCAndFCC(node,nodeType,i1,i2,i3);  
				setSCAndFCC(cell,cellType,i1,i2,i3);
				//BCC boundary: bounce back between node and cell
				setBCC(node,cell,nodeType,cellType,i1,i2,i3); //Boundary btw node and cell
			}
		}
	}			
}

void setInletOutlet(lbgrid &node, lbgrid &cell) {
	myReal feq[27],temp[27],rho,ux,uy,uz,theta,R;
	//inlet
	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
		for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
			copyFromGrid(node.SCP.iB1,i2,i3,node,temp);
			getMoments(temp,rho,ux,uy,uz,theta,R);
			getFeq(feq,rho,u_inlet,0.0,0.0,theta); 
			copyToGrid(node.SCP.iB1-node.SCP.nB1,i2,i3,node,feq);
			copyToGrid(cell.SCP.iB1-cell.SCP.nB1,i2,i3,cell,feq);
		}
	}
	//outlet
	for (int i3 = cell.SCP.iB3; i3 <= cell.SCP.iE3; i3++) {
		for (int i2 = cell.SCP.iB2; i2 <= cell.SCP.iE2; i2++) {
			copyFromGrid(cell.SCP.iE1,i2,i3,cell,temp);
			getMoments(temp,rho,ux,uy,uz,theta,R);
			getFeq(feq,rho,ux,uy,uz,theta); 
			copyToGrid(node.SCP.iE1+node.SCP.nB1,i2,i3,node,feq);
			copyToGrid(cell.SCP.iE1+cell.SCP.nB1,i2,i3,cell,feq);
		}
	}
}


#endif
