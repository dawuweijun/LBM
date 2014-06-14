#ifndef _COLLIDE_H
#define _COLLIDE_H

#include<iostream>
#include<math.h>
#include "grid2D.h"
#include "globals.h"
#include "moments.h"

using namespace std;

void collide(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<1> &nodeZero) {
	myReal temp[9], feq[9], rho, ux, uy;
	//cout<<beta;
	for (int i = 1; i < nx+1; i++) {
		for (int j = 1; j < ny+1; j++) {
			copyFromGrid(i, j, nodeP, nodeM, nodeZero, temp);
			getMoments(temp,rho,ux,uy);
			getFeq(feq,rho,ux,uy);  
			//cout<<rho<<" ";
			for (int k = 0; k < 9; k++) {
				  temp[k] += 2.0*beta*(feq[k] - temp[k]) + (2.0*beta*tau*w[k]*rho*c_x[k]*g_x)/theta;
			}
			copyToGrid(i, j, nodeP, nodeM, nodeZero, temp);
		}
	}
	//cout<<nodeP(nx/2,ny/2,3)<<" ";
}

void collideReplica(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<1> &nodeZero, grid2D<N> &cellP, grid2D<N> &cellM, grid2D<1> &cellZero) {
	myReal temp[9], feq[9], rho, ux, uy;	
	for (int i = 1; i < nx+1; i++) {
		for (int j = 1; j < ny+1; j++) {
			copyFromGrid(i, j, nodeP, nodeM, nodeZero, temp);
			getMoments(temp,rho,ux,uy);
			getFeq(feq,rho,ux,uy);
			for (int k = 0; k < 9; k++) {
				 temp[k] += 2.0*beta*(feq[k] - temp[k]) + (2.0*beta*tau*w[k]*rho*c_x[k]*g_x)/theta;
			}
			copyToGrid(i, j, nodeP, nodeM, nodeZero, temp);
		}
	}
	for (int i = 1; i < nx+1; i++) {
		for (int j = 1; j < ny+1; j++) {
			copyFromGrid(i, j, cellP, cellM, cellZero, temp);
			getMoments(temp,rho,ux,uy);
			getFeq(feq,rho,ux,uy);
			for (int k = 0; k < 9; k++) {
				 temp[k] += 2.0*beta*(feq[k] - temp[k]) + (2.0*beta*tau*w[k]*rho*c_x[k]*g_x)/theta;
			}
			copyToGrid(i, j, cellP, cellM, cellZero, temp);
		}
	}
}
#endif
