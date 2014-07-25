#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include<math.h>
#include<fstream>
#include "gridV.h"
#include "globals3D.h"
#include "moments3D.h"

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

void calForce(myReal del_f, int k,int flag) {
	//if(del_f < 0) {cout<<"in force calculation some dist is neg"<<k<<endl; exit(0);}
	ftempX = del_f*c_x[k];
	ftempY = del_f*c_y[k];
	ftempZ = del_f*c_z[k];
	if (flag == 0) {
		forceXN += ftempX;
		forceYN += ftempY;
		forceZN += ftempZ;
	}
	if (flag == 1) {
		forceXC += ftempX;
		forceYC += ftempY;
		forceZC += ftempZ;
	}
}

void calTorque(double x, double y, double z, int flag) {
	if (flag == 1) { //The point is a cell
		x += 0.5;
		y += 0.5;
		z += 0.5;
	}
	torqueX += (y*ftempZ - z*ftempY);
	torqueY += (z*ftempX - x*ftempZ);
	torqueZ += (x*ftempY - y*ftempX);
}


//Feq for moving
myReal getFeqMoving(int i1, int i2, int i3, lbgrid &node, int k, type flag, arrayRT &rhoA, arrayRT &thetaA, int flag1) {
	myReal rho,theta,feq[27];
	double x,y,z;
	rho = rhoA(i1,i2,i3);
	theta = thetaA(i1,i2,i3);
//	rho = 1.0;
//	theta = 1.0/5.0;
	x = (i1 - center_x);
	y = (i2 - center_y);
	z = (i3 - center_z);
	if (flag1 == 1) { //Point is a cell
		x += 0.5;
		y += 0.5;
		z += 0.5;
	}
	v_sph_x = u_sph_x + (w_sph_y*z - w_sph_z*y);
	v_sph_y = u_sph_y + (w_sph_z*x - w_sph_x*z);
	v_sph_z = u_sph_z + (w_sph_x*y - w_sph_y*x);
	if (flag == SOLID) 
		rho = 1.0;
	
	getFeq(feq,rho,v_sph_x,v_sph_y,v_sph_z,theta);
	return feq[k];
}

//MOVING SC AND FCC*********************************
void setMovingSCAndFCC(lbgrid &node, gridtype &nodeType, int i1, int i2, int i3,arrayRT &rho, arrayRT &theta, int k, int flag) {
	int k_local, k_local_opp;
	myReal del_f,feq[27];
 	if (k>=1 && k <=3) {
//			swap(node.SCM(i1,i2,i3,k-1), node.SCP(i1+c_x[k],i2+c_y[k],i3+c_z[k],k-1));
		node.SCM(i1,i2,i3,k-1) = getFeqMoving(i1,i2,i3,node,(k+3),FLUID,rho,theta,flag);
		node.SCP(i1+c_x[k],i2+c_y[k],i3+c_z[k],k-1) = getFeqMoving(i1,i2,i3,node,k,SOLID,rho,theta,flag);
		del_f = node.SCM(i1,i2,i3,k-1);
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);
	}
	//SCM:4-6
	if (k >=4 && k <=6) {
//			swap(node.SCP(i1,i2,i3,k-4), node.SCM(i1+c_x[k],i2+c_y[k],i3+c_z[k],k-4));
		node.SCP(i1,i2,i3,k-4) = getFeqMoving(i1,i2,i3,node,(k-3),FLUID,rho,theta,flag);
		node.SCM(i1+c_x[k],i2+c_y[k],i3+c_z[k],k-4) = getFeqMoving(i1,i2,i3,node,k,SOLID,rho,theta,flag);
		del_f = node.SCP(i1,i2,i3,k-4);
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);
	}
	//FCC
	//FCC 12:7-10
	if (k >= 7 && k <= 10) {
		k_local = k-7;
		k_local_opp = getFCCReverse(k_local);
//			swap(node.FCC12(i1,i2,i3,k_local_opp), node.FCC12(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local));
		node.FCC12(i1,i2,i3,k_local_opp) = getFeqMoving(i1,i2,i3,node,(k_local_opp+7),FLUID,rho,theta,flag);
		node.FCC12(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local) = getFeqMoving(i1,i2,i3,node,k,SOLID,rho,theta,flag);
		del_f = node.FCC12(i1,i2,i3,k_local_opp);
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);
	}
	//FCC 13:11-14
	if (k >= 11 && k <= 14) {
		k_local = k-11;
		k_local_opp = getFCCReverse(k_local);
//			swap(node.FCC13(i1,i2,i3,k_local_opp), node.FCC13(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local));
		node.FCC13(i1,i2,i3,k_local_opp) = getFeqMoving(i1,i2,i3,node,(k_local_opp+11),FLUID,rho,theta,flag);
		node.FCC13(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local) = getFeqMoving(i1,i2,i3,node,k,SOLID,rho,theta,flag);
		del_f = node.FCC13(i1,i2,i3,k_local_opp);
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);
	}
	//FCC 23:15-19
	if (k >= 15 && k <= 18) {
		k_local = k-15;
		k_local_opp = getFCCReverse(k_local);
//			swap(node.FCC23(i1,i2,i3,k_local_opp), node.FCC23(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local));
		node.FCC23(i1,i2,i3,k_local_opp) = getFeqMoving(i1,i2,i3,node,(k_local_opp+15),FLUID,rho,theta,flag);
		node.FCC23(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local) = getFeqMoving(i1,i2,i3,node,k,SOLID,rho,theta,flag);
		del_f = node.FCC23(i1,i2,i3,k_local_opp);
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);
	}
}

// Moving BCC***************************************************//
void setMovingBCC(lbgrid &node, lbgrid &cell, gridtype &nodeType, gridtype &cellType, int i1, int i2, int i3,arrayRT &rhoN, arrayRT &thetaN,arrayRT &rhoC, arrayRT &thetaC, int k, int flag) {
	int cx, cy, cz;
	myReal del_f;	  
	if (flag == 0) {
	if (k>= 19 && k <= 22) {
//			swap(node.BCCM(i1,i2,i3,k-19), cell.BCCP(i1+cx,i2+cy,i3+cz,k-19));
		node.BCCM(i1,i2,i3,k-19) = getFeqMoving(i1,i2,i3,node,(k+4),FLUID,rhoN,thetaN,flag);
		cell.BCCP(i1+c_x[k]-0.5,i2+c_y[k]-0.5,i3+c_z[k]-0.5,k-19) = getFeqMoving(i1,i2,i3,node,k,SOLID,rhoN,thetaN,flag);
		del_f = node.BCCM(i1,i2,i3,k-19);
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);	
	}
//BCCM:24-27
	if (k >= 23 && k <= 26) {
//			swap(node.BCCP(i1,i2,i3,k-23), cell.BCCM(i1+cx,i2+cy,i3+cz,k-23));
		node.BCCP(i1,i2,i3,k-23) = getFeqMoving(i1,i2,i3,node,(k-4),FLUID,rhoN,thetaN,flag);
		cell.BCCM(i1+c_x[k]-0.5,i2+c_y[k]-0.5,i3+c_z[k]-0.5,k-23) = getFeqMoving(i1,i2,i3,node,k,SOLID,rhoN,thetaN,flag);
		del_f = node.BCCP(i1,i2,i3,k-23);
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);
	}
	}
	else {
	if (k>= 19 && k <= 22) {
//			swap(cell.BCCM(i1,i2,i3,k-19), node.BCCP(i1+cx,i2+cy,i3+cz,k-19));
		cell.BCCM(i1,i2,i3,k-19) = getFeqMoving(i1,i2,i3,cell,(k+4),FLUID,rhoC,thetaC,flag);
		node.BCCP(i1+c_x[k]+0.5,i2+c_y[k]+0.5,i3+c_z[k]+0.5,k-19) = getFeqMoving(i1,i2,i3,cell,k,SOLID,rhoC,thetaC,flag);
		del_f = cell.BCCM(i1,i2,i3,k-19);
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);
	}
	//BCCM:24-27
	if (k >= 23 && k <= 26) {
//			swap(cell.BCCP(i1,i2,i3,k-23), node.BCCM(i1+cx,i2+cy,i3+cz,k-23));
		cell.BCCP(i1,i2,i3,k-23) = getFeqMoving(i1,i2,i3,cell,(k-4),FLUID,rhoC,thetaC,flag);
		node.BCCM(i1+c_x[k]+0.5,i2+c_y[k]+0.5,i3+c_z[k]+0.5,k-23) = getFeqMoving(i1,i2,i3,cell,k,SOLID,rhoC,thetaC,flag);
		del_f = cell.BCCP(i1,i2,i3,k-23);
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);
	}
	}
}
		

void setSCAndFCC(lbgrid &node, gridtype nodeType, int i1, int i2, int i3, int k, int flag) {
	int k_local, k_local_opp;
	myReal del_f;
	//SC swap
	//SCP:1-3
	if (k>=1 && k<=3) {
		swap(node.SCM(i1,i2,i3,k-1), node.SCP(i1+c_x[k],i2+c_y[k],i3+c_z[k],k-1));
		del_f = node.SCM(i1,i2,i3,k-1);
		#ifdef STAT_SPHERE 
			del_f = node.SCM(i1,i2,i3,k-1) - node.SCP(i1+c_x[k],i2+c_y[k],i3+c_z[k],k-1); 
		#endif
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);
	}
	//SCM:4-6
	if (k >=4 && k <=6) {
		swap(node.SCP(i1,i2,i3,k-4), node.SCM(i1+c_x[k],i2+c_y[k],i3+c_z[k],k-4));
		del_f = node.SCP(i1,i2,i3,k-4);
		#ifdef STAT_SPHERE  
			del_f = node.SCP(i1,i2,i3,k-4) - node.SCM(i1+c_x[k],i2+c_y[k],i3+c_z[k],k-4);		
		#endif
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);
	}
	//FCC
	//FCC 12:7-10
	if (k >= 7 && k <= 10) {
		k_local = k-7;
		k_local_opp = getFCCReverse(k_local);
		swap(node.FCC12(i1,i2,i3,k_local_opp), node.FCC12(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local));
		del_f = node.FCC12(i1,i2,i3,k_local_opp);
		#ifdef STAT_SPHERE  
			del_f = node.FCC12(i1,i2,i3,k_local_opp) - node.FCC12(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local);	
		#endif
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);
	}
	//FCC 13:11-14
	if (k >= 11 && k <= 14) {
		k_local = k-11;
		k_local_opp = getFCCReverse(k_local);
		swap(node.FCC13(i1,i2,i3,k_local_opp), node.FCC13(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local));
		del_f = node.FCC13(i1,i2,i3,k_local_opp);
		#ifdef STAT_SPHERE  
			del_f = node.FCC13(i1,i2,i3,k_local_opp) - node.FCC13(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local);
		#endif
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);
	}
	//FCC 23:15-18
	if (k >= 15 && k <= 18) {
		k_local = k-15;
		k_local_opp = getFCCReverse(k_local);
		swap(node.FCC23(i1,i2,i3,k_local_opp), node.FCC23(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local));
		del_f = node.FCC23(i1,i2,i3,k_local_opp);
		#ifdef STAT_SPHERE  
			del_f = node.FCC23(i1,i2,i3,k_local_opp) - node.FCC23(i1+c_x[k],i2+c_y[k],i3+c_z[k],k_local);		
		#endif
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);
	  
	}
	
}
void setBCC(lbgrid &node, lbgrid &cell, gridtype &nodeType, gridtype &cellType, int i1, int i2, int i3, int k, int flag) {
	//Btw node and cell
	//BCCP:20-23
	int cx, cy, cz;
	myReal del_f;
	if (flag == 0) {
	if (k>= 19 && k <= 22) {
		swap(node.BCCM(i1,i2,i3,k-19), cell.BCCP(i1+c_x[k]-0.5,i2+c_y[k]-0.5,i3+c_z[k]-0.5,k-19));
		del_f = node.BCCM(i1,i2,i3,k-19);
		#ifdef STAT_SPHERE 
			del_f = node.BCCM(i1,i2,i3,k-19) - cell.BCCP(i1+c_x[k]-0.5,i2+c_y[k]-0.5,i3+c_z[k]-0.5,k-19);	 	
		#endif
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);
		
	}
	//BCCM:24-27
	if (k >= 23 && k <= 26) {
		swap(node.BCCP(i1,i2,i3,k-23), cell.BCCM(i1+c_x[k]-0.5,i2+c_y[k]-0.5,i3+c_z[k]-0.5,k-23));
		del_f = node.BCCP(i1,i2,i3,k-23);
		#ifdef STAT_SPHERE 
			del_f = node.BCCP(i1,i2,i3,k-23) - cell.BCCM(i1+c_x[k]-0.5,i2+c_y[k]-0.5,i3+c_z[k]-0.5,k-23);
		#endif
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);
	
	}
	}
	else {
	if (k>= 19 && k <= 22) {
		swap(cell.BCCM(i1,i2,i3,k-19), node.BCCP(i1+c_x[k]+0.5,i2+c_y[k]+0.5,i3+c_z[k]+0.5,k-19));
		del_f = cell.BCCM(i1,i2,i3,k-19);
		#ifdef STAT_SPHERE  
			del_f = cell.BCCM(i1,i2,i3,k-19) - node.BCCP(i1+c_x[k]+0.5,i2+c_y[k]+0.5,i3+c_z[k]+0.5,k-19);
		#endif
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);
		
	}
	//BCCM:24-27
	if (k >= 23 && k <= 26) {
		swap(cell.BCCP(i1,i2,i3,k-23), node.BCCM(i1+c_x[k]+0.5,i2+c_y[k]+0.5,i3+c_z[k]+0.5,k-23));
		del_f = cell.BCCP(i1,i2,i3,k-23);
		#ifdef STAT_SPHERE  
			del_f = cell.BCCP(i1,i2,i3,k-23) - node.BCCM(i1+c_x[k]+0.5,i2+c_y[k]+0.5,i3+c_z[k]+0.5,k-23);
		#endif
		calForce(del_f,k,flag);
		calTorque((i1-center_x),(i2-center_y),(i3-center_z),flag);
	
	}
	}
}

void getRhoAndTheta(int i1, int i2, int i3, lbgrid &node,arrayRT &rhoNode, arrayRT &thetaNode) {
	myReal temp[27], rho, theta, ux, uy, uz, R, x, y, z;
	copyFromGrid(i1,i2,i3,node,temp);
	getMoments(temp,rho,ux,uy,uz,theta,R);
	rhoNode(i1,i2,i3) = rho;
	thetaNode(i1,i2,i3) = theta;
}

void setFeq(lbgrid &node, int i1, int i2, int i3) {
	myReal temp[27],rho,ux,uy,uz,R,theta;
	copyFromGrid(i1,i2,i3,node,temp);
	getMoments(temp,rho,ux,uy,uz,theta,R);
	getFeq(temp,1.0,u_inlet,uy,uz,theta);
	copyToGrid(i1,i2,i3,node,temp);
}

void setBoundaryCorrection(lbgrid &node, lbgrid &cell, gridtype nodeType, gridtype cellType) {
	int i2,i3,i1;
	i1 = node.SCP.iE1+1;
	for (i2 = node.SCP.iB1; i2 <= node.SCP.iE2; i2++) {
		for (i3 = node.SCP.iB2; i3 <= node.SCP.iE3; i3++) {
			for(int k=1; k<19; k++) {
				if (k == 1 || k == 7 || k == 9 || k == 11 || k == 12) break;
            			if(nodeType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) {
					setSCAndFCC(node,nodeType,i1,i2,i3,k,0);  
					setFeq(node,i1,i2,i3);
				}
				
				if(cellType(i1,i2,i3) == FLUID && cellType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) {
					setSCAndFCC(cell,cellType,i1,i2,i3,k,1); 
					setFeq(cell,i1,i2,i3);
				}
			}
			for(int k=19; k<27; k++) {
				if (k == 19 || k == 21 || k == 24 || k == 26) break;
				if(nodeType(i1,i2,i3) == FLUID && cellType(i1+c_x[k]-0.5,i2+c_y[k]-0.5,i3+c_z[k]-0.5) == SOLID) {
					setBCC(node,cell,nodeType,cellType,i1,i2,i3,k,0); 
					setFeq(node,i1,i2,i3);
				}
				if(cellType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k]+0.5,i2+c_y[k]+0.5,i3+c_z[k]+0.5) == SOLID) {
					setBCC(node,cell,nodeType,cellType,i1,i2,i3,k,1);
					setFeq(cell,i1,i2,i3);
				}
			}
		}
	}
	i1 = node.SCP.iB1-1;
	for (i2 = node.SCP.iB1; i2 <= node.SCP.iE2; i2++) {
		for (i3 = node.SCP.iB2; i3 <= node.SCP.iE3; i3++) {
			for(int k=1; k<19; k++) {
				if (k == 4 || k == 8 || k == 10 || k == 13 || k == 14) break;
            			if(nodeType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) {
					setSCAndFCC(node,nodeType,i1,i2,i3,k,0);  
					setFeq(node,i1,i2,i3);
				}
				if(cellType(i1,i2,i3) == FLUID && cellType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) {
					setSCAndFCC(cell,cellType,i1,i2,i3,k,1); 
					setFeq(cell,i1,i2,i3);
				}
			}
			for(int k=19; k<27; k++) {
				if (k == 20 || k == 22 || k == 23 || k == 25) break;
				if(nodeType(i1,i2,i3) == FLUID && cellType(i1+c_x[k]-0.5,i2+c_y[k]-0.5,i3+c_z[k]-0.5) == SOLID) {
					setBCC(node,cell,nodeType,cellType,i1,i2,i3,k,0);
					setFeq(node,i1,i2,i3);
				} 
				if(cellType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k]+0.5,i2+c_y[k]+0.5,i3+c_z[k]+0.5) == SOLID) {
					setBCC(node,cell,nodeType,cellType,i1,i2,i3,k,1);
					setFeq(cell,i1,i2,i3);
				}
			}
		}
	}
}
	

void setInternalBoundary(lbgrid &node, lbgrid &cell, gridtype nodeType, gridtype cellType, arrayRT &rhoNode, arrayRT &thetaNode, arrayRT &rhoCell, arrayRT &thetaCell) {
	forceXN = forceYN = forceZN = 0.0;
	forceXC = forceYC = forceZC = 0.0;
	torqueX = torqueY = torqueZ = 0.0;
	b = 0;
	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
		for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
			for (int i1 = node.SCP.iB1; i1 <= node.SCP.iE1; i1++) {
				for(int k=1; k<19; k++) {
            				if(nodeType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID)
					      setSCAndFCC(node,nodeType,i1,i2,i3,k,0);  
					if(cellType(i1,i2,i3) == FLUID && cellType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID)
					      setSCAndFCC(cell,cellType,i1,i2,i3,k,1); 
				}
				for(int k=19; k<27; k++) {
					if(nodeType(i1,i2,i3) == FLUID && cellType(i1+c_x[k]-0.5,i2+c_y[k]-0.5,i3+c_z[k]-0.5) == SOLID)
					      setBCC(node,cell,nodeType,cellType,i1,i2,i3,k,0); 
					if(cellType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k]+0.5,i2+c_y[k]+0.5,i3+c_z[k]+0.5) == SOLID)
					      setBCC(node,cell,nodeType,cellType,i1,i2,i3,k,1);
				}

                        }
		}
	}
//CHECK:FASTER--BUT NOT CORRECT
// 	int flag = 0;
// 	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
// 		for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
// 			for (int i1 = node.SCP.iB1; i1 <= node.SCP.iE1; i1++) {
// 				for(int k=1; k<19; k++) {
// 					if(nodeType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) flag = 1;
// 					if(nodeType(i1,i2,i3) == SOLID && nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == FLUID) flag = 1;
// 					if (flag == 1) break;
// 				}
// 				if (flag == 1) getRhoAndTheta(i1,i2,i3,node,rhoNode,thetaNode);
// 				for(int k=1; k<19; k++) {
// 					if(cellType(i1,i2,i3) == FLUID && cellType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID) flag = 1;
// 					if(cellType(i1,i2,i3) == SOLID && cellType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == FLUID) flag = 1;
// 					if (flag == 1) break;
// 				}
// 				if (flag == 1) getRhoAndTheta(i1,i2,i3,cell,rhoCell,thetaCell);
// 				for(int k=19; k<27; k++) {
// 					if(nodeType(i1,i2,i3) == FLUID && cellType(i1+c_x[k]-0.5,i2+c_y[k]-0.5,i3+c_z[k]-0.5) == SOLID) flag = 1;
// 					if (flag == 1) break;
// 					if(cellType(i1,i2,i3) == SOLID && nodeType(i1+c_x[k]+0.5,i2+c_y[k]+0.5,i3+c_z[k]+0.5) == FLUID) flag = 2;
// 					if (flag == 2) break;
// 				}
// 				if (flag == 1) getRhoAndTheta(i1,i2,i3,node,rhoNode,thetaNode);
// 				if (flag == 2) getRhoAndTheta(i1,i2,i3,cell,rhoCell,thetaCell);
// 				for(int k=19; k<27; k++) {
// 					if(cellType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k]+0.5,i2+c_y[k]+0.5,i3+c_z[k]+0.5) == SOLID) flag = 1;
// 					if (flag == 1) break;
// 					if(nodeType(i1,i2,i3) == SOLID && cellType(i1+c_x[k]-0.5,i2+c_y[k]-0.5,i3+c_z[k]-0.5) == FLUID) flag = 2;
// 					if (flag == 2) break;
// 				}
// 				if (flag == 1)	getRhoAndTheta(i1,i2,i3,cell,rhoCell,thetaCell);
// 				if (flag == 2)	getRhoAndTheta(i1,i2,i3,node,rhoNode,thetaNode);
// 					
// 
// 			}
// 		}
// 	}

	#ifdef MOVING_SPHERE	
	b = 1;
	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
        	for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
            		for (int i1 = node.SCP.iB1; i1 <= node.SCP.iE1; i1++) {
                		getRhoAndTheta(i1,i2,i3,node,rhoNode,thetaNode);
                		getRhoAndTheta(i1,i2,i3,cell,rhoCell,thetaCell);
            		}
        	}
    	}
	
	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
		for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
			for (int i1 = node.SCP.iB1; i1 <= node.SCP.iE1; i1++) {
				for(int k=1; k<19; k++) {
            				if(nodeType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID)
					      setMovingSCAndFCC(node,nodeType,i1,i2,i3,rhoNode,thetaNode,k,0);
					if(cellType(i1,i2,i3) == FLUID && cellType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == SOLID)
					      setMovingSCAndFCC(cell,cellType,i1,i2,i3,rhoCell,thetaCell,k,1);
				}
				for(int k=19; k<27; k++) {
					if(nodeType(i1,i2,i3) == FLUID && cellType(i1+c_x[k]-0.5,i2+c_y[k]-0.5,i3+c_z[k]-0.5) == SOLID)
					      setMovingBCC(node,cell,nodeType,cellType,i1,i2,i3,rhoNode,thetaNode,rhoCell,thetaCell,k,0);
					if(cellType(i1,i2,i3) == FLUID && nodeType(i1+c_x[k]+0.5,i2+c_y[k]+0.5,i3+c_z[k]+0.5) == SOLID)
					      setMovingBCC(node,cell,nodeType,cellType,i1,i2,i3,rhoNode,thetaNode,rhoCell,thetaCell,k,1); 
				}
			}
		}
	}
	double V = 4.187*radius*radius*radius;
//	if (u_sph_x > 0.0 || u_sph_y > 0.0 || u_sph_z > 0.0) {
//		u_sph_x = u_sph_x + (forceX+forceXOld)/(2*10*V);   //Taking sphere density as 10
//		u_sph_y = u_sph_y + (forceY+forceYOld)/(2*10*V);
//		u_sph_z = u_sph_z + (forceZ+forceZOld)/(2*10*V);
//	}
//	if (w_sph_z > 0.0 || w_sph_x > 0.0 || w_sph_y > 0.0) {
//	 	w_sph_x = w_sph_x + 2.5*(torqueX+torqueXOld)/(2*radius*radius*V*10);
//	 	w_sph_y += 2.5*(torqueY+torqueYOld)/(2*radius*radius*V*10);
//	 	w_sph_z += 2.5*(torqueZ+torqueZOld)/(2*radius*radius*V*10);
//	}
//	if (u_sph_x > 0.0 || u_sph_y > 0.0 || u_sph_z > 0.0) {
//		u_sph_x = u_sph_x + (forceX)/(10*V);   //Taking sphere density as 10
//		u_sph_y = u_sph_y + (forceY)/(10*V);
//		u_sph_z = u_sph_z + (forceZ)/(10*V);
//	}
//	if (w_sph_z > 0.0 || w_sph_x > 0.0 || w_sph_y > 0.0) {
//	 	w_sph_x = w_sph_x + 2.5*(torqueX)/(radius*radius*V*10);
//	 	w_sph_y += 2.5*(torqueY)/(radius*radius*V*10);
//	 	w_sph_z += 2.5*(torqueZ)/(radius*radius*V*10);
//	}
	#endif
	
}

void getGradRefill(int i1, int i2, int i3,myReal fgrad[27],lbgrid &node) {
	myReal temp[27],f[27],rho,ux,uy,uz,theta,R,pxx,pyy,pzz,pxy,pxz,pyz;
	copyFromGrid(i1,i2,i3,node,temp);
	getMoments(temp,rho,ux,uy,uz,theta,R);
	getP(temp,pxx,pyy,pzz,pxy,pxz,pyz);
	getGrad(f,rho,ux,uy,uz,pxx,pyy,pzz,pxy,pxz,pyz);
	for (int k = 0; k < 27; k++) 
		fgrad[k] += f[k];
}

void refill(gridtype nodeTypeOld, gridtype cellTypeOld, gridtype nodeType, gridtype cellType, lbgrid &node, lbgrid &cell) {
	myReal fgrad[27];
	for(int k = 0; k < 27; k++) fgrad[k] = 0.0;
	int count = 0;
	int cx, cy, cz;
	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
	for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
	for (int i1 = node.SCP.iB1; i1 <= node.SCP.iE1; i1++) {
		for(int k = 0; k < 27; k++) fgrad[k] = 0.0;
		count = 0;
		//If nodetype has changed, look around all fluid particles and average their grad and put it back in node
		if (nodeTypeOld(i1,i2,i3) == SOLID && nodeType(i1,i2,i3) == FLUID) {
			for (int k = 1; k < 19; k++) { //SC and FCC
				if (nodeType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == FLUID) {
					getGradRefill(i1+c_x[k],i2+c_y[k],i3+c_z[k],fgrad,node);
					count++; 
				}
			}
			for (int k = 19; k < 27; k++) { //BCC
				cx = (c_x[k]>0?0:2*c_x[k]); 
				cy = (c_y[k]>0?0:2*c_y[k]);
				cz = (c_z[k]>0?0:2*c_z[k]);
				if (cellType(i1+cx,i2+cy,i3+cz) == FLUID) {
					getGradRefill(i1+cx,i2+cy,i3+cz,fgrad,cell);
					count++;
				}
			}
			for (int k = 0; k < 27; k++) fgrad[k] /= count;
			copyToGrid(i1,i2,i3,node,fgrad);

		}

	}
	}
	}


	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
	for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
	for (int i1 = node.SCP.iB1; i1 <= node.SCP.iE1; i1++) {
	//If cell type has changed
		count = 0;
		for(int k = 0; k < 27; k++) fgrad[k] = 0.0;
		if (cellTypeOld(i1,i2,i3) == SOLID && cellType(i1,i2,i3) == FLUID) {
			for (int k = 1; k < 19; k++)  {//SC and FCC
				if (cellType(i1+c_x[k],i2+c_y[k],i3+c_z[k]) == FLUID) {
					getGradRefill(i1+c_x[k],i2+c_y[k],i3+c_z[k],fgrad,cell);
					count++; 
				}
			}
			for (int k = 19; k < 27; k++) { //BCC
				cx = (c_x[k]<0?0:2*c_x[k]); 
				cy = (c_y[k]<0?0:2*c_y[k]);
				cz = (c_z[k]<0?0:2*c_z[k]);
				if (nodeType(i1+cx,i2+cy,i3+cz) == FLUID) {
					getGradRefill(i1+cx,i2+cy,i3+cz,fgrad,node);
					count++;
				}
			}
			for (int k = 0; k < 27; k++) fgrad[k] /= count;
			copyToGrid(i1,i2,i3,cell,fgrad);
		}
	
	}
	}
	}
}

void updateCMObj(gridtype &nodeTypeOld, gridtype &nodeType, gridtype &cellTypeOld, gridtype &cellType) {
	for (int i3 = 0; i3 < nodeType.m3; i3++) {
		for (int i2 = 0; i2 < nodeType.m2; i2++) {
			for (int i1 = 0; i1 < nodeType.m1; i1++) {
				nodeTypeOld(i1,i2,i3) = nodeType(i1,i2,i3);
				cellTypeOld(i1,i2,i3) = cellType(i1,i2,i3);
			}
		}
	}
	center_x += u_sph_x*dt;
	center_y += u_sph_y*dt;
	center_z += u_sph_z*dt;
}
			

//************************************************************************************************************************************
//INLET_OUTLET

void setInletOutlet(lbgrid &node, lbgrid &cell) {
	myReal feq[27],temp[27],rho,ux,uy,uz,theta,R,pxx,pyy,pzz,pxy,pxz,pyz;
	myReal rho1,ux1,uy1,uz1,theta1,R1,pxx1,pyy1,pzz1,pxy1,pxz1,pyz1;
	//inlet-nodes
	for (int i3 = node.SCP.iB3; i3 < node.SCP.m3-1; i3++) {
		for (int i2 = node.SCP.iB2; i2 < node.SCP.m2-1; i2++) {
			copyFromGrid(node.SCP.iB1,i2,i3,node,temp);
			getMoments(temp,rho,ux,uy,uz,theta,R);
			getP(temp,pxx,pyy,pzz,pxy,pxz,pyz);
			rho = 1.0;
			//getGrad(feq,rho,u_inlet,uy,uz,pxx,pyy,pzz,pxy,pxz,pyz);
			getFeq(feq,rho,u_inlet,uy,uz,theta); 
			for (int nBi = 0 ; nBi < node.SCP.iB1; nBi++) {
				copyToGrid(nBi,i2,i3,node,feq);
				//copyToGrid(nBi,i2,i3,cell,feq);
			}
			//Cell correction
			if( i2 != node.SCP.iE2){
				copyFromGrid(node.SCP.iB1,i2+1,i3,node,temp);
				getMoments(temp,rho1,ux1,uy1,uz1,theta1,R1);
				rho1 = 1.0;
				rho = avg(rho,rho1);
				ux = avg(ux,ux1);
				uy = avg(uy,uy1);
				uz = avg(uz,uz1);
				theta = avg(theta,theta1);
				pxx = avg(pxx1,pxx);
				pyy = avg(pyy1,pyy);
				pzz = avg(pzz1,pzz);
				pxy = avg(pxy1,pxy);
				pxz = avg(pxz1,pxz);
				pyz = avg(pyz1,pyz);
			}
//			getGrad(feq,rho,u_inlet,uy,uz,pxx,pyy,pzz,pxy,pxz,pyz);
			getFeq(feq,rho,ux,uy,uz,theta); 
			for (int nBi = 0 ; nBi < node.SCP.iB1; nBi++) {
				copyToGrid(nBi,i2,i3,cell,feq);
			}
		}
	}
	//outlet-cells
	for (int i3 = node.SCP.iB3; i3 < node.SCP.m3-1; i3++) {
		for (int i2 = node.SCP.iB2; i2 < node.SCP.m2-1; i2++) {
			copyFromGrid(cell.SCP.iE1,i2,i3,cell,temp);
			getMoments(temp,rho,ux,uy,uz,theta,R);
			getP(temp,pxx,pyy,pzz,pxy,pxz,pyz);
//			getGrad(feq,rho,ux,uy,uz,pxx,pyy,pzz,pxy,pxz,pyz);
			getFeq(feq,rho,ux,uy,uz,theta); 
			for (int nBi = node.SCP.iE1+1; nBi < node.SCP.m1; nBi++) {
				copyToGrid(nBi,i2,i3,cell,feq);
			}
		
			//Node correction
			if(i2 != node.SCP.iB2) {
	 			copyFromGrid(cell.SCP.iE1,i2-1,i3,cell,temp);
				getMoments(temp,rho1,ux1,uy1,uz1,theta1,R1);
				getP(temp,pxx1,pyy1,pzz1,pxy1,pxz1,pyz1);
			
				rho = avg(rho,rho1);
				ux = avg(ux,ux1);
				uy = avg(uy,uy1);
				uz = avg(uz,uz1);
				theta = avg(theta,theta1);
				pxx = avg(pxx1,pxx);
				pyy = avg(pyy1,pyy);
				pzz = avg(pzz1,pzz);
				pxy = avg(pxy1,pxy);
				pxz = avg(pxz1,pxz);
				pyz = avg(pyz1,pyz);
			
			}
//			getGrad(feq,rho,ux,uy,uz,pxx,pyy,pzz,pxy,pxz,pyz);
			getFeq(feq,rho,ux,uy,uz,theta); 
			for (int nBi = node.SCP.iE1+1 ; nBi < node.SCP.m1; nBi++) {
				copyToGrid(nBi,i2,i3,node,feq);
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
			node.BCCP(iB1,i2,i3,DV_P1_P1_P1) = node.BCCP(0,i2,i3,DV_P1_P1_P1);
			node.BCCP(iB1,i2,i3,DV_P1_M1_P1) = node.BCCP(0,i2,i3,DV_P1_M1_P1);
			node.BCCM(iB1,i2,i3,DV_P1_P1_M1) = node.BCCM(0,i2,i3,DV_P1_P1_M1);
			node.BCCM(iB1,i2,i3,DV_P1_M1_M1) = node.BCCM(0,i2,i3,DV_P1_M1_M1);
			
			//FCC
			cell.FCC12(iB1,i2,i3,DV_P1_P1_ZERO) = cell.FCC12(0,i2,i3,DV_P1_P1_ZERO);
			cell.FCC12(iB1,i2,i3,DV_P1_M1_ZERO) = cell.FCC12(0,i2,i3,DV_P1_M1_ZERO);
			cell.FCC13(iB1,i2,i3,DV_P1_ZERO_P1) = cell.FCC13(0,i2,i3,DV_P1_ZERO_P1);
			cell.FCC13(iB1,i2,i3,DV_P1_ZERO_M1) = cell.FCC13(0,i2,i3,DV_P1_ZERO_M1);
			 
			
		
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
			cell.BCCP(iE1,i2,i3,DV_M1_P1_P1) = cell.BCCP(iE1+1,i2,i3,DV_M1_P1_P1);
			cell.BCCP(iE1,i2,i3,DV_M1_M1_P1) = cell.BCCP(iE1+1,i2,i3,DV_M1_M1_P1);
			cell.BCCM(iE1,i2,i3,DV_M1_P1_M1) = cell.BCCM(iE1+1,i2,i3,DV_M1_P1_M1);
			cell.BCCM(iE1,i2,i3,DV_M1_M1_M1) = cell.BCCM(iE1+1,i2,i3,DV_M1_M1_M1);
			//FCC
			node.FCC12(iE1,i2,i3,DV_M1_P1_ZERO) = node.FCC12(iE1+1,i2,i3,DV_M1_P1_ZERO);
			node.FCC12(iE1,i2,i3,DV_M1_M1_ZERO) = node.FCC12(iE1+1,i2,i3,DV_M1_M1_ZERO);
			node.FCC13(iE1,i2,i3,DV_M1_ZERO_P1) = node.FCC13(iE1+1,i2,i3,DV_M1_ZERO_P1);
			node.FCC13(iE1,i2,i3,DV_M1_ZERO_M1) = node.FCC13(iE1+1,i2,i3,DV_M1_ZERO_M1);	
		
		}
	}
//	//For Real cells
//	//Inlet
//	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
//	for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
//	cell.SCP(iB1,i2,i3,DV_P1_ZERO_ZERO) = avg(node.SCP(0,i2,i3,DV_P1_ZERO_ZERO),node.SCP(0,i2+1,i3,DV_P1_ZERO_ZERO));		
//	cell.FCC12(iB1,i2,i3,DV_P1_P1_ZERO) = avg(node.FCC12(0,i2,i3,DV_P1_P1_ZERO),node.FCC12(0,i2+1,i3,DV_P1_P1_ZERO));
//	cell.FCC12(iB1,i2,i3,DV_P1_M1_ZERO) = avg(node.FCC12(0,i2,i3,DV_P1_M1_ZERO),node.FCC12(0,i2+1,i3,DV_P1_M1_ZERO));
//	cell.FCC13(iB1,i2,i3,DV_P1_ZERO_P1) = avg(node.FCC13(0,i2,i3,DV_P1_ZERO_P1),node.FCC13(0,i2+1,i3,DV_P1_ZERO_P1));
//	cell.FCC13(iB1,i2,i3,DV_P1_ZERO_M1) = avg(node.FCC13(0,i2,i3,DV_P1_ZERO_M1),node.FCC13(0,i2+1,i3,DV_P1_ZERO_M1));
//	}
//	}
//	//outlet
//	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
//	for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
//	node.SCM(iE1,i2,i3,DV_M1_ZERO_ZERO) = avg(cell.SCM(iE1+1,i2,i3,DV_M1_ZERO_ZERO),cell.SCM(iE1+1,i2-1,i3,DV_M1_ZERO_ZERO));		
//	node.FCC12(iE1,i2,i3,DV_M1_P1_ZERO) = avg(cell.FCC12(iE1+1,i2,i3,DV_M1_P1_ZERO),cell.FCC12(iE1+1,i2-1,i3,DV_M1_P1_ZERO));
//	node.FCC12(iE1,i2,i3,DV_M1_M1_ZERO) = avg(cell.FCC12(iE1+1,i2,i3,DV_M1_M1_ZERO),cell.FCC12(iE1+1,i2-1,i3,DV_M1_M1_ZERO));
//	node.FCC13(iE1,i2,i3,DV_M1_ZERO_P1) = avg(cell.FCC13(iE1+1,i2,i3,DV_M1_ZERO_P1),cell.FCC13(iE1+1,i2-1,i3,DV_M1_ZERO_P1));
//	node.FCC13(iE1,i2,i3,DV_M1_ZERO_M1) = avg(cell.FCC13(iE1+1,i2,i3,DV_M1_ZERO_M1),cell.FCC13(iE1+1,i2-1,i3,DV_M1_ZERO_M1));
//	}
//	}

}
			
//void setInletOutlet(lbgrid &node, lbgrid &cell) {
//	myReal feq[27],temp[27],rho,ux,uy,uz,theta,R,pxx,pyy,pzz,pxy,pxz,pyz;
//	//inlet
//	for (int i3 = node.SCP.iB3-1; i3 < node.SCP.m3; i3++) {
//		for (int i2 = node.SCP.iB2-1; i2 < node.SCP.m2; i2++) {
//			copyFromGrid(node.SCP.iB1,i2,i3,node,temp);
//			getMoments(temp,rho,ux,uy,uz,theta,R);
//			getP(temp,pxx,pyy,pzz,pxy,pxz,pyz);
//			getGrad(feq,1.0,u_inlet,uy,uz,pxx,pyy,pzz,pxy,pxz,pyz);
////			getFeq(feq,rho,u_inlet,uy,uz,theta); 
//			for (int nBi = 0 ; nBi < node.SCP.iB1; nBi++) {
//				copyToGrid(nBi,i2,i3,node,feq);
//				copyToGrid(nBi,i2,i3,cell,feq);
//			}
//		}
//	}
//	//outlet
//	for (int i3 = node.SCP.iB3-1; i3 < node.SCP.m3; i3++) {
//		for (int i2 = node.SCP.iB2-1; i2 < node.SCP.m2; i2++) {
//			copyFromGrid(cell.SCP.iE1,i2,i3,cell,temp);
//			getMoments(temp,rho,ux,uy,uz,theta,R);
//			getP(temp,pxx,pyy,pzz,pxy,pxz,pyz);
//			getGrad(feq,rho,ux,uy,uz,pxx,pyy,pzz,pxy,pxz,pyz);
////			getFeq(feq,rho,ux,uy,uz,theta); 
//			for (int nBi = node.SCP.iE1+1 ; nBi < node.SCP.m1; nBi++) {
//				copyToGrid(nBi,i2,i3,node,feq);
//				copyToGrid(nBi,i2,i3,cell,feq);
//			}
//		}
//	}
//}

//void setInletOutletCorrection(lbgrid &node, lbgrid &cell) {
//	int iB1 = node.SCP.iB1;
//	int iE1 = node.SCP.iE1;
//	myReal temp[27];
//	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
//		for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
//			//FCC
//			node.FCC12(iB1,i2,i3,DV_P1_P1_ZERO) = node.FCC12(0,i2,i3,DV_P1_P1_ZERO);
//			node.FCC12(iB1,i2,i3,DV_P1_M1_ZERO) = node.FCC12(0,i2,i3,DV_P1_M1_ZERO);
//			node.FCC13(iB1,i2,i3,DV_P1_ZERO_P1) = node.FCC13(0,i2,i3,DV_P1_ZERO_P1);
//			node.FCC13(iB1,i2,i3,DV_P1_ZERO_M1) = node.FCC13(0,i2,i3,DV_P1_ZERO_M1);
//			//BCC
//			node.BCCP(iB1,i2,i3,DV_P1_P1_P1) = cell.BCCP(0,i2,i3,DV_P1_P1_P1);
//			node.BCCP(iB1,i2,i3,DV_P1_M1_P1) = cell.BCCP(0,i2,i3,DV_P1_M1_P1);
//			node.BCCM(iB1,i2,i3,DV_P1_P1_M1) = cell.BCCM(0,i2,i3,DV_P1_P1_M1);
//			node.BCCM(iB1,i2,i3,DV_P1_M1_M1) = cell.BCCM(0,i2,i3,DV_P1_M1_M1);
//			
//		
//		}
//	}
//	//outlet
//	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
//		for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
//			//FCC
//			cell.FCC12(iE1,i2,i3,DV_M1_P1_ZERO) = cell.FCC12(iE1+1,i2,i3,DV_M1_P1_ZERO);
//			cell.FCC12(iE1,i2,i3,DV_M1_M1_ZERO) = cell.FCC12(iE1+1,i2,i3,DV_M1_M1_ZERO);
//			cell.FCC13(iE1,i2,i3,DV_M1_ZERO_P1) = cell.FCC13(iE1+1,i2,i3,DV_M1_ZERO_P1);
//			cell.FCC13(iE1,i2,i3,DV_M1_ZERO_M1) = cell.FCC13(iE1+1,i2,i3,DV_M1_ZERO_M1);
//			//BCC
//			cell.BCCP(iE1,i2,i3,DV_M1_P1_P1) = node.BCCP(iE1+1,i2,i3,DV_M1_P1_P1);
//			cell.BCCP(iE1,i2,i3,DV_M1_M1_P1) = node.BCCP(iE1+1,i2,i3,DV_M1_M1_P1);
//			cell.BCCM(iE1,i2,i3,DV_M1_P1_M1) = node.BCCM(iE1+1,i2,i3,DV_M1_P1_M1);
//			cell.BCCM(iE1,i2,i3,DV_M1_M1_M1) = node.BCCM(iE1+1,i2,i3,DV_M1_M1_M1);
//			
//		
//		}
//	}

//	//For Real cells
//	//Inlet
//	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
//	for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
//	cell.SCP(iB1,i2,i3,DV_P1_ZERO_ZERO) = avg(node.SCP(0,i2,i3,DV_P1_ZERO_ZERO),node.SCP(0,i2+1,i3,DV_P1_ZERO_ZERO));		
//	cell.FCC12(iB1,i2,i3,DV_P1_P1_ZERO) = avg(node.FCC12(0,i2,i3,DV_P1_P1_ZERO),node.FCC12(0,i2+1,i3,DV_P1_P1_ZERO));
//	cell.FCC12(iB1,i2,i3,DV_P1_M1_ZERO) = avg(node.FCC12(0,i2,i3,DV_P1_M1_ZERO),node.FCC12(0,i2+1,i3,DV_P1_M1_ZERO));
//	cell.FCC13(iB1,i2,i3,DV_P1_ZERO_P1) = avg(node.FCC13(0,i2,i3,DV_P1_ZERO_P1),node.FCC13(0,i2+1,i3,DV_P1_ZERO_P1));
//	cell.FCC13(iB1,i2,i3,DV_P1_ZERO_M1) = avg(node.FCC13(0,i2,i3,DV_P1_ZERO_M1),node.FCC13(0,i2+1,i3,DV_P1_ZERO_M1));
//	}
//	}
//	//outlet
//	for (int i3 = node.SCP.iB3; i3 <= node.SCP.iE3; i3++) {
//	for (int i2 = node.SCP.iB2; i2 <= node.SCP.iE2; i2++) {
//	node.SCM(iE1,i2,i3,DV_M1_ZERO_ZERO) = avg(cell.SCM(iE1+1,i2,i3,DV_M1_ZERO_ZERO),cell.SCM(iE1+1,i2-1,i3,DV_M1_ZERO_ZERO));		
//	node.FCC12(iE1,i2,i3,DV_M1_P1_ZERO) = avg(cell.FCC12(iE1+1,i2,i3,DV_M1_P1_ZERO),cell.FCC12(iE1+1,i2-1,i3,DV_M1_P1_ZERO));
//	node.FCC12(iE1,i2,i3,DV_M1_M1_ZERO) = avg(cell.FCC12(iE1+1,i2,i3,DV_M1_M1_ZERO),cell.FCC12(iE1+1,i2-1,i3,DV_M1_M1_ZERO));
//	node.FCC13(iE1,i2,i3,DV_M1_ZERO_P1) = avg(cell.FCC13(iE1+1,i2,i3,DV_M1_ZERO_P1),cell.FCC13(iE1+1,i2-1,i3,DV_M1_ZERO_P1));
//	node.FCC13(iE1,i2,i3,DV_M1_ZERO_M1) = avg(cell.FCC13(iE1+1,i2,i3,DV_M1_ZERO_M1),cell.FCC13(iE1+1,i2-1,i3,DV_M1_ZERO_M1));
//	}
//	}


//}

#endif
