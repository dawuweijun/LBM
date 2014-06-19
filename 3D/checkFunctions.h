#ifndef _CHECK_FUNCTIONS_H
#define _CHECK_FUNCTIONS_H

#include<fstream>

#include "gridV.h"
#include "globals3D.h"
#include "moments3D.h"
using namespace std;


void checkInit(lbgrid &node) {
	int n1, n2, n3; 
	n1 = node.SCP.m1;
	n2 = node.SCP.m2;
	n3 = node.SCP.m3;
	for (int i3 = 1; i3 < n3-1; i3++) {
		for (int i2 = 1; i2 < n2-1; i2++) {
			for (int i1 = 1; i1 < n1-1; i1++) {
				for (int k = 0; k < N; k++) {
					if (node.FCC12(i1,i2,i3,k)==0||node.FCC13(i1,i2,i3,k)==0|| node.FCC23(i1,i2,i3,k)==0||node.BCCP(i1,i2,i3,k)==0||node.BCCM(i1,i2,i3,k)==0)
						cout<<"INITIALIZE IS MAKING A MISTAKE"<<endl;
				}
				for (int k = 0; k < N-1; k++) {
					if (node.SCP(i1,i2,i3,k)==0||node.SCM(i1,i2,i3,k)==0)
						cout<<"INITIALIZE IS MAKING A MISTAKE IN SC"<<endl;
				}
			}
		}
	}
}

void printGlobalRho(lbgrid &node, lbgrid &cell) {
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
				global_rho += rho;//(ux*ux + uy*uy + uz*uz)+theta;
				global_rho1 += rho*(ux*ux + uy*uy + uz*uz+3.0*theta);
			}
		}
	}
	//for (int k = 0; k < 27; k++) 
		//cout<<checkf[k]/(2*NX*NY*NZ)<<endl;
	cout<<global_rho/(2*NX*NY*NZ)<<endl;       //<<" ;GLOBAL THETA:"<<global_rho1/(2*NX*NY*NZ)<<endl;
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

void printAllVel(lbgrid &node) {
	myReal temp[27], rho, ux, uy, uz, theta, R;
	ofstream file;
	file.open("AllVelocities3D.vtk");
	for (int i3 = node.SCM.iB3; i3 <= node.SCM.iE3; i3++) {
		for (int i2 = node.SCM.iB2; i2 <= node.SCM.iE2; i2++) {
			for (int i1 = node.SCM.iB1; i1 <= node.SCM.iE1; i1++) {
				copyFromGrid(i1,i2,i3,node,temp);
				getMoments(temp,rho,ux,uy,uz,theta,R);
				file<<ux-g_x/2<<" "<<uy<<" "<<uz;
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
					file<<ux<<" "<<uy<<" "<<uz;
					file<<endl;
				}
			}
		}
		file.close();
	}
}

#endif
