#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include<math.h>
#include<fstream>
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
			data = new myReal[size];
			for (int i = 0; i < size; i++) data[i] = FLUID;
		}
		int n1,n2,n3,m1,m2,m3,size;
		int iE1,iE2,iE3,iB1,iB2,iB3;
		int nB1,nB2,nB3;
		myReal *data;
		//Overload the () operator
		myReal operator() (int i, int j, int k) const {
			return data[k*m1*m2 + j*m1 + i];
		}
		myReal& operator() (int i, int j, int k) {
			return data[k*m1*m2 + j*m1 + i];
		}
};

void setTypeWall(gridtype nodeType, gridtype cellType) {
	for (int i3 = 3; i3 < nodeType.m3-3; i3++) {
		for (int i1 = 3; i1 < nodeType.m1-3; i1++) {
			for (int index = cellType.iE2-1;  index < cellType.m2; index++) 
				cellType(i1,index,i3) = nodeType(i1,index,i3) = SOLID; //Solid at the top-cells
			for (int index = nodeType.iB2+1;  index >= 0; index--) 
				cellType(i1,index,i3) = nodeType(i1,index,i3) = SOLID; //Solid at the bottom-nodes
		}
	}
//	ofstream file;
//	file.open("type.txt");
//	for (int i2 = 0; i2 < nodeType.m2; i2++) {
//		for (int i1 = 0; i1 < nodeType.m1; i1++) {
//			file<<nodeType(i1,i2,NZ/2);
//		}
//		file<<endl;
//	}
//	file.close();
}
void setTypeSphere(gridtype nodeType, gridtype cellType) {
	for (int i3 = 0; i3 < nodeType.m3; i3++) {
		for (int i2 = 0; i2 < nodeType.m2; i2++) {
			for (int i1 = 0; i1 < nodeType.m1; i1++) {
				if (pow(i1-center_x,2) + pow(i2-center_y,2) + pow(i3-center_z,2) <= pow(radius,2)) 
					nodeType(i1,i2,i3) = SOLID;
				else
					nodeType(i1,i2,i3) = FLUID;
			}
		}
	}
	for (int i3 = 0; i3 < cellType.m3; i3++) {
		for (int i2 = 0; i2 < cellType.m2; i2++) {
			for (int i1 = 0; i1 < cellType.m1; i1++) {
				if (pow(i1+0.5-center_x,2) + pow(i2+0.5-center_y,2) + pow(i3+0.5-center_z,2) <= pow(radius,2))  
					cellType(i1,i2,i3) = SOLID;
				else
					cellType(i1,i2,i3) = FLUID;
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

void setSCAndFCC(lbgrid &node, gridtype &nodeType, int i1, int i2, int i3) {
	int k_local, k_local_opp;
	//SC swap
	//SCP:1-3
	for (int k = 1; k <= 3; k++)
		if (nodeType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) {
			swap(node.SCM(i1,i2,i3,k-1), node.SCP(i1+c_x[k],i2+c_y[k],i3+c_z[k],k-1));
			forceX += (node.SCM(i1,i2,i3,k-1)-node.SCP(i1+c_x[k],i2+c_y[k],i3+c_z[k],k-1))*c_x[k];
		}
	//SCM:4-6
	for (int k = 4; k <= 6; k++)
		if (nodeType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) {
			swap(node.SCP(i1,i2,i3,k-4), node.SCM(i1+c_x[k],i2+c_y[k],i3+c_z[k],k-4));
			forceX += (node.SCP(i1,i2,i3,k-4)-node.SCM(i1+c_x[k],i2+c_y[k],i3+c_z[k],k-4))*c_x[k];
		}
	//FCC
	//FCC 12:7-10
	for (int k = 7; k <=10; k++) {
		k_local = k-7;
		if (nodeType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) { 
			k_local_opp = getFCCReverse(k_local);
			swap(node.FCC12(i1,i2,i3,k_local_opp), node.FCC12(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local));
			forceX += (node.FCC12(i1,i2,i3,k_local_opp)-node.FCC12(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local))*c_x[k];
		}
	}
	//FCC 13:11-14
	for (int k = 11; k <=14; k++) {
		k_local = k-11;
		if (nodeType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) { 
			k_local_opp = getFCCReverse(k_local);
			swap(node.FCC13(i1,i2,i3,k_local_opp), node.FCC13(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local));
			forceX += (node.FCC13(i1,i2,i3,k_local_opp)-node.FCC13(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local))*c_x[k];
		}
	}
	//FCC 23:15-19
	for (int k = 15; k <=18; k++) {
		k_local = k-15;
		if (nodeType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) { 
			k_local_opp = getFCCReverse(k_local);
			swap(node.FCC23(i1,i2,i3,k_local_opp), node.FCC23(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local));
			forceX += (node.FCC23(i1,i2,i3,k_local_opp)-node.FCC23(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local))*c_x[k];
		}
	}
	
}

void setBCC(lbgrid &node, lbgrid &cell, gridtype &nodeType, gridtype &cellType, int i1, int i2, int i3) {
	//Btw node and cell
	//BCCP:20-23
	int cx, cy, cz;
	for (int k = 19; k <=22; k++) {
		cx = (c_x[k]>0?0:2*c_x[k]); //Replica grid is displaced by 0.5 units in each direction
		cy = (c_y[k]>0?0:2*c_y[k]);
		cz = (c_z[k]>0?0:2*c_z[k]);
		if (nodeType(i1,i2,i3) == FLUID && cellType(i1+cx,i2+cy,i3+cz) == SOLID) {
			//cout<<"CHECK"<<i1<<i2<<i3<<endl;
			swap(node.BCCM(i1,i2,i3,k-19), cell.BCCP(i1+cx,i2+cy,i3+cz,k-19));	
			forceX += (node.BCCM(i1,i2,i3,k-19)-cell.BCCP(i1+cx,i2+cy,i3+cz,k-19))*c_x[k];
		}	
	}
	//BCCM:24-27
	for (int k = 23; k <=26; k++) {
		cx = (c_x[k]>0?0:2*c_x[k]); 
		cy = (c_y[k]>0?0:2*c_y[k]);
		cz = (c_z[k]>0?0:2*c_z[k]);
		if (nodeType(i1,i2,i3) == FLUID && cellType(i1+cx,i2+cy,i3+cz) == SOLID) {
			swap(node.BCCP(i1,i2,i3,k-23), cell.BCCM(i1+cx,i2+cy,i3+cz,k-23));
			forceX += (node.BCCP(i1,i2,i3,k-23)-cell.BCCM(i1+cx,i2+cy,i3+cz,k-23))*c_x[k];
		}
	}
	//Btw cell and node
	//BCCP
	for (int k = 19; k <=22; k++) {
		cx = (c_x[k]<0?0:2*c_x[k]); //Replica grid is displaced by 0.5 units in each direction
		cy = (c_y[k]<0?0:2*c_y[k]);
		cz = (c_z[k]<0?0:2*c_z[k]);
		if (cellType(i1,i2,i3) == FLUID && nodeType(i1+cx,i2+cy,i3+cz) == SOLID) {
			swap(cell.BCCM(i1,i2,i3,k-19), node.BCCP(i1+cx,i2+cy,i3+cz,k-19));
			forceX += (cell.BCCM(i1,i2,i3,k-19), node.BCCP(i1+cx,i2+cy,i3+cz,k-19))*c_x[k];
		}
	}
	//BCCM:24-27
	for (int k = 23; k <=26; k++) {
		cx = (c_x[k]<0?0:2*c_x[k]); 
		cy = (c_y[k]<0?0:2*c_y[k]);
		cz = (c_z[k]<0?0:2*c_z[k]);
		if (cellType(i1,i2,i3) == FLUID && nodeType(i1+cx,i2+cy,i3+cz) == SOLID) {
			swap(cell.BCCP(i1,i2,i3,k-23), node.BCCM(i1+cx,i2+cy,i3+cz,k-23));
			forceX += (cell.BCCP(i1,i2,i3,k-23), node.BCCM(i1+cx,i2+cy,i3+cz,k-23))*c_x[k];
		}
	}
}

void setInternalBoundary(lbgrid &node, lbgrid &cell, gridtype nodeType, gridtype cellType) {
	forceX = 0.0;
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
	//cout<<forceX/(0.5*u_inlet*u_inlet*3.14*radius*radius)<<endl;		
}

void setInletOutlet(lbgrid &node, lbgrid &cell) {
	myReal feq[27],temp[27],rho,ux,uy,uz,theta,R,pxx,pyy,pzz,pxy,pxz,pyz;
	//inlet
	for (int i3 = node.SCP.iB3-1; i3 < node.SCP.m3; i3++) {
		for (int i2 = node.SCP.iB2-1; i2 < node.SCP.m2; i2++) {
			copyFromGrid(node.SCP.iB1,i2,i3,node,temp);
			getMoments(temp,rho,ux,uy,uz,theta,R);
			getP(temp,pxx,pyy,pzz,pxy,pxz,pyz);
			getGrad(feq,1.0,u_inlet,uy,uz,pxx,pyy,pzz,pxy,pxz,pyz);
//			getFeq(feq,rho,u_inlet,uy,uz,theta); 
			for (int nBi = 0 ; nBi < node.SCP.iB1; nBi++) {
				copyToGrid(nBi,i2,i3,node,feq);
				copyToGrid(nBi,i2,i3,cell,feq);
			}
		}
	}
	//outlet
	for (int i3 = node.SCP.iB3-1; i3 < node.SCP.m3; i3++) {
		for (int i2 = node.SCP.iB2-1; i2 < node.SCP.m2; i2++) {
			copyFromGrid(cell.SCP.iE1,i2,i3,cell,temp);
			getMoments(temp,rho,ux,uy,uz,theta,R);
			getP(temp,pxx,pyy,pzz,pxy,pxz,pyz);
			getGrad(feq,rho,ux,uy,uz,pxx,pyy,pzz,pxy,pxz,pyz);
//			getFeq(feq,rho,ux,uy,uz,theta); 
			for (int nBi = node.SCP.iE1+1 ; nBi < node.SCP.m1; nBi++) {
				copyToGrid(nBi,i2,i3,node,feq);
				copyToGrid(nBi,i2,i3,cell,feq);
			}
		}
	}
}

void setInletOutletCorrection(lbgrid &node, lbgrid &cell) {
	int iB1 = node.SCP.iB1;
	int iE1 = node.SCP.iE1;
	myReal temp[27];
	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
		for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
			//FCC
			node.FCC12(iB1,i2,i3,DV_P1_P1_ZERO) = node.FCC12(0,i2,i3,DV_P1_P1_ZERO);
			node.FCC12(iB1,i2,i3,DV_P1_M1_ZERO) = node.FCC12(0,i2,i3,DV_P1_M1_ZERO);
			node.FCC13(iB1,i2,i3,DV_P1_ZERO_P1) = node.FCC13(0,i2,i3,DV_P1_ZERO_P1);
			node.FCC13(iB1,i2,i3,DV_P1_ZERO_M1) = node.FCC13(0,i2,i3,DV_P1_ZERO_M1);
			//BCC
			node.BCCP(iB1,i2,i3,DV_P1_P1_P1) = cell.BCCP(0,i2,i3,DV_P1_P1_P1);
			node.BCCP(iB1,i2,i3,DV_P1_M1_P1) = cell.BCCP(0,i2,i3,DV_P1_M1_P1);
			node.BCCM(iB1,i2,i3,DV_P1_P1_M1) = cell.BCCM(0,i2,i3,DV_P1_P1_M1);
			node.BCCM(iB1,i2,i3,DV_P1_M1_M1) = cell.BCCM(0,i2,i3,DV_P1_M1_M1);
			
		
		}
	}
	//outlet
	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
		for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
			//FCC
			cell.FCC12(iE1,i2,i3,DV_M1_P1_ZERO) = cell.FCC12(iE1+1,i2,i3,DV_M1_P1_ZERO);
			cell.FCC12(iE1,i2,i3,DV_M1_M1_ZERO) = cell.FCC12(iE1+1,i2,i3,DV_M1_M1_ZERO);
			cell.FCC13(iE1,i2,i3,DV_M1_ZERO_P1) = cell.FCC13(iE1+1,i2,i3,DV_M1_ZERO_P1);
			cell.FCC13(iE1,i2,i3,DV_M1_ZERO_M1) = cell.FCC13(iE1+1,i2,i3,DV_M1_ZERO_M1);
			//BCC
			cell.BCCP(iE1,i2,i3,DV_M1_P1_P1) = node.BCCP(iE1+1,i2,i3,DV_M1_P1_P1);
			cell.BCCP(iE1,i2,i3,DV_M1_M1_P1) = node.BCCP(iE1+1,i2,i3,DV_M1_M1_P1);
			cell.BCCM(iE1,i2,i3,DV_M1_P1_M1) = node.BCCM(iE1+1,i2,i3,DV_M1_P1_M1);
			cell.BCCM(iE1,i2,i3,DV_M1_M1_M1) = node.BCCM(iE1+1,i2,i3,DV_M1_M1_M1);
			
		
		}
	}
	//Hack--not sure
//	for (int i3 = node.SCP.iB3; i3 < node.SCP.m3-1; i3++) {
//		for (int i2 = node.SCP.iB2; i2 < node.SCP.m2-1; i2++) {
//			copyFromGrid(0,i2,i3,node,temp);
//			copyToGrid(iB1,i2,i3,node,temp);
//		}
//	}
//	//outlet
//	for (int i3 = node.SCP.iB3; i3 < node.SCP.m3-1; i3++) {
//		for (int i2 = node.SCP.iB2; i2 < node.SCP.m2-1; i2++) {
//			copyFromGrid(iE1+1,i2,i3,cell,temp);
//			copyToGrid(iE1,i2,i3,cell,temp);	
//		}
//	}
	//For Real cells
	//Inlet
	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
	for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
	cell.SCP(iB1,i2,i3,DV_P1_ZERO_ZERO) = avg(node.SCP(0,i2,i3,DV_P1_ZERO_ZERO),node.SCP(0,i2+1,i3,DV_P1_ZERO_ZERO));		
	cell.FCC12(iB1,i2,i3,DV_P1_P1_ZERO) = avg(node.FCC12(0,i2,i3,DV_P1_P1_ZERO),node.FCC12(0,i2+1,i3,DV_P1_P1_ZERO));
	cell.FCC12(iB1,i2,i3,DV_P1_M1_ZERO) = avg(node.FCC12(0,i2,i3,DV_P1_M1_ZERO),node.FCC12(0,i2+1,i3,DV_P1_M1_ZERO));
	cell.FCC13(iB1,i2,i3,DV_P1_ZERO_P1) = avg(node.FCC13(0,i2,i3,DV_P1_ZERO_P1),node.FCC13(0,i2+1,i3,DV_P1_ZERO_P1));
	cell.FCC13(iB1,i2,i3,DV_P1_ZERO_M1) = avg(node.FCC13(0,i2,i3,DV_P1_ZERO_M1),node.FCC13(0,i2+1,i3,DV_P1_ZERO_M1));
	}
	}
	//outlet
	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
	for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
	node.SCM(iE1,i2,i3,DV_M1_ZERO_ZERO) = avg(cell.SCM(iE1+1,i2,i3,DV_M1_ZERO_ZERO),cell.SCM(iE1+1,i2-1,i3,DV_M1_ZERO_ZERO));		
	node.FCC12(iE1,i2,i3,DV_M1_P1_ZERO) = avg(cell.FCC12(iE1+1,i2,i3,DV_M1_P1_ZERO),cell.FCC12(iE1+1,i2-1,i3,DV_M1_P1_ZERO));
	node.FCC12(iE1,i2,i3,DV_M1_M1_ZERO) = avg(cell.FCC12(iE1+1,i2,i3,DV_M1_M1_ZERO),cell.FCC12(iE1+1,i2-1,i3,DV_M1_M1_ZERO));
	node.FCC13(iE1,i2,i3,DV_M1_ZERO_P1) = avg(cell.FCC13(iE1+1,i2,i3,DV_M1_ZERO_P1),cell.FCC13(iE1+1,i2-1,i3,DV_M1_ZERO_P1));
	node.FCC13(iE1,i2,i3,DV_M1_ZERO_M1) = avg(cell.FCC13(iE1+1,i2,i3,DV_M1_ZERO_M1),cell.FCC13(iE1+1,i2-1,i3,DV_M1_ZERO_M1));
	}
	}


}
			


#endif
