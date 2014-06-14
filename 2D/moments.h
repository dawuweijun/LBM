#ifndef _MOMENTS_H
#define _MOMENTS_H
#include "globals.h"

void copyFromGrid(int i, int j, grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<1> &nodeZero, myReal temp[]) {
	temp[0] = nodeZero(i,j,0);
	for (int k = 1; k < 5; k++) {
		temp[k] = nodeP(i,j,k-1);
	}
	for (int k = 5; k < 9; k++) {
		temp[k] = nodeM(i,j,k-5);
	}
}

void getMoments(myReal f[9], myReal& rho, myReal& u_x, myReal& u_y) {
	rho = u_x = u_y = 0.0;
	for (int k = 0 ; k < 9 ; k++) {
		rho += f[k];
		u_x += f[k]*c_x[k];
		u_y += f[k]*c_y[k];
	}
	u_x /= rho;
	u_y /= rho;
	u_x += g_x/2;
}	

void getFeq(myReal f_eq[9], myReal rho, myReal u_x, myReal u_y) {
	myReal u_mod, dot_prod;
	u_mod = pow(u_x,2) + pow(u_y,2);
	for (int k = 0; k < 9; k++) {
		dot_prod = c_x[k]*u_x + c_y[k]*u_y;
		#ifdef GRID_STAGGERED
		f_eq[k] = w[k]*rho*(1.0 + 6.0*dot_prod + (36.0/2.0)*pow(dot_prod,2) - (6.0/2.0)*u_mod);
		#else
		f_eq[k] = w[k]*rho*(1.0 + 3.0*dot_prod + (9.0/2.0)*pow(dot_prod,2) - (3.0/2.0)*u_mod);
		#endif
	}
}

void copyToGrid(int i, int j, grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<1> &nodeZero, myReal temp[]) {
	nodeZero(i,j,0) = temp[0];
	for (int k = 1; k < 5; k++) {
		nodeP(i,j,k-1) = temp[k];
	}
	for (int k = 5; k < 9; k++) {
		nodeM(i,j,k-5) = temp[k];
	} 
}

void getGrad(myReal f_grad[], myReal temp[], myReal rho, myReal u_x, myReal u_y) {
	myReal p_xx, p_yy, p_xy, dot_prod;
	p_xx = p_yy = p_xy = 0;
	for (int k = 0 ; k < 9 ; k++) {
		p_xx += temp[k]*c_x[k]*c_x[k];
		p_yy += temp[k]*c_y[k]*c_y[k];
		p_xy += temp[k]*c_x[k]*c_y[k];
	}
	for (int k = 0 ; k < 9; k++) {
		dot_prod = c_x[k]*u_x + c_y[k]*u_y;	
     		f_grad[k] = w[k]*(rho + 3.0*rho*dot_prod + (9.0/2.0)*((p_xx-rho/3.0)*(c_x[k]*c_x[k] - 1.0/3.0) + 2.0*c_x[k]*c_y[k]*p_xy + (c_y[k]*c_y[k] - 1.0/3.0)*(p_yy-rho/3.0)));
	}

}

#endif
