#ifndef _CHECK_FUNCTIONS_H
#define _CHECK_FUNCTIONS_H

#include <fstream>
#include <math.h>

#include "gridV.h"
#include "globals3D.h"
#include "moments3D.h"
using namespace std;


void checkInit(lbgrid &node, gridtype nodeType) {
	int n1, n2, n3; 
	n1 = node.SCP.m1;
	n2 = node.SCP.m2;
	n3 = node.SCP.m3;
	for (int i3 = 1; i3 < n3-1; i3++) {
		for (int i2 = 1; i2 < n2-1; i2++) {
			for (int i1 = 1; i1 < n1-1; i1++) {
				for (int k = 0; k < N; k++) {
					if (node.FCC12(i1,i2,i3,k)<=0||node.FCC13(i1,i2,i3,k)<=0|| node.FCC23(i1,i2,i3,k)<=0||node.BCCP(i1,i2,i3,k)<=0||node.BCCM(i1,i2,i3,k)<=0)
						 {//if (nodeType(i1,i2,i3) == SOLID)break;
						cout<<"type:"<<nodeType(i1,i2,i3)<<" "<<k<<"AT"<<i1<<" "<<i2<<" "<<i3<<endl; exit(0);}
				}
				for (int k = 0; k < N-1; k++) {
					if (node.SCP(i1,i2,i3,k)<=0||node.SCM(i1,i2,i3,k)<=0)
						{//if (nodeType(i1,i2,i3) == SOLID)break; 
					cout<<"type:"<<nodeType(i1,i2,i3)<<" "<<k<<"AT"<<i1<<" "<<i2<<" "<<i3<<endl; exit(0);}
				}
			}
		}
	}
}

void printGlobalRho(lbgrid &node, lbgrid &cell, gridtype nodeType, gridtype cellType, std::ofstream &file) {
	myReal temp[27], checkf[27], rho, ux, uy, uz, theta, R, global_rho, global_rho1;
	global_rho = 0;
	global_rho1 = 0;
	for (int k = 0; k < 27; k++) checkf[k] = 0.0;
	for (int i3 = node.SCM.iB3; i3 <= node.SCM.iE3; i3++) {
		for (int i2 = node.SCM.iB2; i2 <= node.SCM.iE2; i2++) {
			for (int i1 = node.SCM.iB1; i1 <= node.SCM.iE1; i1++) {
				copyFromGrid(i1,i2,i3,node,temp);
				getMoments(temp,rho,ux,uy,uz,theta,R);
				for (int k = 0; k < 27; k++)
					checkf[k] += temp[k];
//				if(nodeType(i1,i2,i3) == SOLID) rho = 0;
				global_rho += rho;
				global_rho1 += rho*(ux*ux + uy*uy + uz*uz+3.0*theta);
			}
		}
	}
	for (int i3 = node.SCM.iB3; i3 <= node.SCM.iE3; i3++) {
		for (int i2 = node.SCM.iB2; i2 <= node.SCM.iE2; i2++) {
			for (int i1 = node.SCM.iB1; i1 <= node.SCM.iE1; i1++) {
				copyFromGrid(i1,i2,i3,cell,temp);
				getMoments(temp,rho,ux,uy,uz,theta,R);
				for (int k = 0; k < 27; k++)	
					checkf[k] += temp[k];
//				if(cellType(i1,i2,i3) == SOLID) rho = 0;
				global_rho += rho;//(ux*ux + uy*uy + uz*uz)+theta;
				global_rho1 += rho*(ux*ux + uy*uy + uz*uz+3.0*theta);
			}
		}
	}
// 	for (int k = 0; k < 27; k++) 
// 		cout<<checkf[k]/(2*NX*NY*NZ)<<endl;
	//cout<<global_rho/(2*NX*NY*NZ)<<endl;// <<" ;GLOBAL THETA:"<<global_rho1/(2*NX*NY*NZ)<<endl;
	file<<global_rho/(2*NX*NY*NZ)<<endl;
	//cout<<"ITERATION OVER"<<endl;
}

void printEnergy(lbgrid &node, std::ofstream &file) {
	myReal temp[27], rho, ux, uy, uz, theta, R, E;
	E = 0;
	for (int i3 = node.SCM.iB3; i3 <= node.SCM.iE3; i3++) {
		for (int i2 = node.SCM.iB2; i2 <= node.SCM.iE2; i2++) {
			for (int i1 = node.SCM.iB1; i1 <= node.SCM.iE1; i1++) {
				copyFromGrid(i1,i2,i3,node,temp);
				getMoments(temp,rho,ux,uy,uz,theta,R);
				E += (ux*ux + uy*uy + uz*uz) * 0.5;
			}
		}
	}
	file<<E<<endl;
}

void printAllVel(lbgrid &node, gridtype nodeType) {
	myReal temp[27], rho, ux, uy, uz, theta, R;
	ofstream file;
	file.open("AllVelocities3D2.vtk");
	file<<"# vtk DataFile Version 3.0"<<endl;
	file<<"Velocity"<<endl;
	file<<"ASCII"<<endl;
	file<<"DATASET STRUCTURED_POINTS"<<endl;
	file<<"DIMENSIONS "<<NX<<" "<<NY<<" "<<NZ<<endl; 
	file<<"ORIGIN 0 0 0"<<endl;
	file<<"SPACING 1 1 1"<<endl;
	file<<"POINT_DATA "<<NX*NY*NZ<<endl;
	file<<"VECTORS velocity double"<<endl;
	for (int i3 = node.SCM.iB3; i3 <= node.SCM.iE3; i3++) {
		for (int i2 = node.SCM.iB2; i2 <= node.SCM.iE2; i2++) {
			for (int i1 = node.SCM.iB1; i1 <= node.SCM.iE1; i1++) {
				copyFromGrid(i1,i2,i3,node,temp);
				getMoments(temp,rho,ux,uy,uz,theta,R);
				if (nodeType(i1,i2,i3) == SOLID) file<<0.0<<" "<<0.0<<" "<<0.0;
				else file<<ux<<" "<<uy<<" "<<uz;
				file<<endl;
			}
		}
	}	
	file.close();
}

void printEachTime(lbgrid &node, int t, int &c) {
	if (t%100 == 0){ 
		myReal temp[27], rho, ux, uy, uz, theta, R;
		ofstream file;
		char name[20];
		sprintf(name,"CHeckVelocity%d.vtk",++c);
		file.open(name);
		for (int i3 = node.SCM.iB3; i3 <= node.SCM.iE3; i3++) {
			for (int i2 = node.SCM.iB2; i2 <= node.SCM.iE2; i2++) {
				for (int i1 = node.SCM.iB1; i1 <= node.SCM.iE1; i1++) {
					copyFromGrid(i1,i2,i3,node,temp);
					getMoments(temp,rho,ux,uy,uz,theta,R);
//					file<<ux<<" "<<uy<<" "<<uz;
//					file<<endl;
					file<<rho<<endl;
				}
			}
		}
		file.close();
	}
}

void printVelX(lbgrid &node, gridtype nodeType, int t, int &c) {
	if (t%5000 == 0){ 
		myReal temp[27], rho, ux, uy, uz, theta, R;
		ofstream file;
		char name[20];
		sprintf(name,"CHeckVelocity%d.txt",++c);
		file.open(name);
//		for (int i3 = node.SCM.iB3; i3 <= node.SCM.iE3; i3++) {
			for (int i2 = node.SCM.iB2; i2 <= node.SCM.iE2; i2++) {
//				for (int i1 = node.SCM.iB1; i1 <= node.SCM.iE1; i1++) {

					copyFromGrid(NX/2,i2,NZ/2,node,temp);
					getMoments(temp,rho,ux,uy,uz,theta,R);
					if (nodeType(NX/2,i2,NZ/2) == SOLID) ux = 0;
					file<<ux<<endl;	
				}
//			}
//		}
		file.close();
	}
}

void printCenterLocation(ofstream &f) {
	f<<center_x<<" "<<center_y<<" "<<center_z<<endl;
}

void printForceAndTorque(ofstream &f) {//(0.5*u_ref*u_ref*3.142*radius*radius)
	//f<<"TORQUE:"  << torqueX<<" "<<torqueY<<" "<<torqueZ/(0.5*w_ref*w_ref*3.14*radius*radius*radius*radius*radius)<<endl;
	f<<(forceXN+forceXC)/2.0<<" "<<(forceYN+forceYC)/2.0<<" "<<(forceZN+forceZC)/2.0<<endl;
}

void printU(ofstream &f) {
	f<<sqrt(u_sph_x*u_sph_x + u_sph_y*u_sph_y + u_sph_z*u_sph_z)/u_ref<<endl;
}

void printOmega(ofstream &f) {
	f<<sqrt(w_sph_x*w_sph_x + w_sph_y*w_sph_y + w_sph_z*w_sph_z)/w_ref<<endl;
}

void printStokesCorrection() {
	myReal factor, zeta, Fd;
	factor = 1/radius - 2.837/NY + 4.19*radius*radius/(NY*NY*NY) - 27.4*pow(radius,5)/pow(NY,6) ;
	zeta = 6*3.1428*nu/factor;
	Fd = zeta*u_ref;
	cout<<"Stokes correction:"<<Fd<<endl;
} 
#endif
