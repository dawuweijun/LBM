#ifndef _MOMENTS_3D_H
#define _MOMENTS_3D_H
#include <iostream>
#include <math.h>
#include "gridV.h"
#include "globals3D.h"

void getMoments(myReal f[27], myReal& rho, myReal& u_x, myReal& u_y, myReal& u_z, myReal& theta, myReal &R) {
	rho = u_x = u_y = u_z = theta = R = 0.0;
	for (int k = 0 ; k < 27 ; k++) {
		rho += f[k];
		u_x += f[k]*c_x[k];
		u_y += f[k]*c_y[k];
		u_z += f[k]*c_z[k];
		theta += f[k]*(pow(c_x[k],2) + pow(c_y[k],2) + pow(c_z[k],2))/3.0;
		R += f[k]*(pow(c_x[k],2) + pow(c_y[k],2) + pow(c_z[k],2))*(pow(c_x[k],2) + pow(c_y[k],2) + pow(c_z[k],2));
	}
	u_x /= rho;
	u_y /= rho;
	u_z /= rho;
	u_x += g_x/2;
	theta = theta - rho*(pow(u_x,2) + pow(u_y,2) + pow(u_z,2))/3.0;	
	R /= rho;
}

void copyFromGrid(int i1, int i2, int i3, lbgrid& node, myReal f[27]) {
	int index = 0;
	f[index] = node.Zero(i1,i2,i3,0);
	index += 1;
	for (int k = 0; k < N-1; k++) 
		f[index++] = node.SCP(i1,i2,i3,k);
	for (int k = 0; k < N-1; k++) 
		f[index++] = node.SCM(i1,i2,i3,k);
	for (int k = 0; k < N; k++) 
		f[index++] = node.FCC12(i1,i2,i3,k);
	for (int k = 0; k < N; k++) 
		f[index++] = node.FCC13(i1,i2,i3,k);
	for (int k = 0; k < N; k++) 
		f[index++] = node.FCC23(i1,i2,i3,k);
	for (int k = 0; k < N; k++) 
		f[index++] = node.BCCP(i1,i2,i3,k);
	for (int k = 0; k < N; k++) 
		f[index++] = node.BCCM(i1,i2,i3,k);
}

void copyToGrid(int i1, int i2, int i3, lbgrid& node, myReal f[27]) {
	int index = 0;
	f[index] = node.Zero(i1,i2,i3,0);
	index += 1;
	for (int k = 0; k < N-1; k++) 
		node.SCP(i1,i2,i3,k) = f[index++];
	for (int k = 0; k < N-1; k++) 
		node.SCM(i1,i2,i3,k) = f[index++];
	for (int k = 0; k < N; k++) 
		node.FCC12(i1,i2,i3,k) = f[index++];
	for (int k = 0; k < N; k++) 
		node.FCC13(i1,i2,i3,k) = f[index++];
	for (int k = 0; k < N; k++) 
		node.FCC23(i1,i2,i3,k) = f[index++];
	for (int k = 0; k < N; k++) 
		node.BCCP(i1,i2,i3,k) = f[index++];
	for (int k = 0; k < N; k++) 
		node.BCCM(i1,i2,i3,k) = f[index++];
}

void getFeq(myReal feq[27], myReal rho, myReal u_x, myReal u_y, myReal u_z, myReal theta, myReal R) {
	myReal h, f_tilda,c,u,dot_prod;
	h = (15*theta*theta - R)/(9*theta*theta - R);
	u = (pow(u_x,2) + pow(u_y,2) + pow(u_z,2));
	for (int k = 0; k < 27; k++) {
		c = (pow(c_x[k],2) + pow(c_y[k],2) + pow(c_z[k],2));
		dot_prod = c_x[k]*u_x + c_y[k]*u_y + c_z[k]*u_z;
		f_tilda = w[k]*rho*(1 + (theta-theta_0)*(c*c/theta_0 - 3.0)/2.0);
		feq[k] = f_tilda*(1 + dot_prod/theta + 0.5*pow(dot_prod/theta,2) - (u*u/(2*theta))*(1-h) - h*u*u*c*c/(6*theta*theta));
	}
}


//Collision
void collide(lbgrid& node) {
	myReal temp[27], feq[27], rho, ux, uy, uz, theta, R;	
	for (int i3 = node.SCM.iB3; i3 <= node.SCM.iE1; i3++) {
		for (int i2 = node.SCM.iB2; i2 <= node.SCM.iE2; i2++) {
			for (int i1 = node.SCM.iB1; i1 <= node.SCM.iE1; i1++) {
				copyFromGrid(i1,i2,i3,node,temp);
				getMoments(temp,rho,ux,uy,uz,theta,R);
				getFeq(feq,rho,ux,uy,uz,theta,R);
				for (int k = 0; k < 27; k++) {
					 temp[k] += 2.0*beta*(feq[k] - temp[k]) + (2.0*beta*tau*w[k]*rho*c_x[k]*g_x)/theta_0;
				}
				copyToGrid(i1,i2,i3,node,temp);
			}
		}
	}
}
	

#endif
