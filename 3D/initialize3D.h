#ifndef _INIT_3D_H
#define _INIT_3D_H
#include <iostream>
#include<cstdlib>
#include <math.h>
#include "gridV.h"
#include "globals3D.h"
#include "moments3D.h"
using namespace std;


void initializeKida(lbgrid &node, lbgrid &cell) {
	int n1, n2, n3;
	myReal x,y,z,u1,u2,u3,rho,feq[27]; 
	n1 = node.SCP.m1;
	n2 = node.SCP.m2;
	n3 = node.SCP.m3;
	//Nodes
	for (int i3 = 1; i3 < n3-1; i3++) {
		for (int i2 = 1; i2 < n2-1; i2++) {
			for (int i1 = 1; i1 < n1-1; i1++) {
				x = (i1-0.5)*6.28/NX;
				y = (i2-0.5)*6.28/NY;
				z = (i3-0.5)*6.28/NZ;
				u1 = u_inlet*sin(x)*(cos(3.0*y)*cos(z) - cos(y)*cos(3.0*z));//;+drand48();
         			u2 = u_inlet*sin(y)*(cos(3.0*z)*cos(x) - cos(z)*cos(3.0*x));//;+drand48();
         			u3 = u_inlet*sin(z)*(cos(3.0*x)*cos(y) - cos(x)*cos(3.0*y));//;+drand48();
				getFeq(feq,1.0,u1,u2,u3,theta_0); 
				copyToGrid(i1,i2,i3,node,feq);
			}
		}
	}
	//cells
	for (int i3 = 1; i3 < n3-1; i3++) {
		for (int i2 = 1; i2 < n2-1; i2++) {
			for (int i1 = 1; i1 < n1-1; i1++) {
				x = (i1)*6.28/NX;
				y = (i2)*6.28/NY;
				z = (i3)*6.28/NZ;
				u1 = u_inlet*sin(x)*(cos(3.0*y)*cos(z) - cos(y)*cos(3.0*z));//+drand48();
         			u2 = u_inlet*sin(y)*(cos(3.0*z)*cos(x) - cos(z)*cos(3.0*x));//+drand48();
         			u3 = u_inlet*sin(z)*(cos(3.0*x)*cos(y) - cos(x)*cos(3.0*y));//+drand48();
				getFeq(feq,1.0,u1,u2,u3,theta_0); 
				copyToGrid(i1,i2,i3,cell,feq);
			}
		}
	}
}


void initializeWithZeroVel(lbgrid &node, lbgrid &cell) {
	int n1, n2, n3, x, y, z;
	myReal u1,u2,u3,rho,feq[27]; 
	n1 = node.SCP.m1;
	n2 = node.SCP.m2;
	n3 = node.SCP.m3;
	//nodes
	for (int i3 = 1; i3 < n3-1; i3++) {
		for (int i2 = 1; i2 < n2-1; i2++) {
			for (int i1 = 1; i1 < n1-1; i1++) {
				getFeq(feq,1.0,0.0,0.0,0.0,theta_0); 
				copyToGrid(i1,i2,i3,node,feq);
			}
		}
	}
	//cells
	for (int i3 = 1; i3 < n3-1; i3++) {
		for (int i2 = 1; i2 < n2-1; i2++) {
			for (int i1 = 1; i1 < n1-1; i1++) {
				getFeq(feq,1.0,0.0,0.0,0.0,theta_0); 
				copyToGrid(i1,i2,i3,cell,feq);
			}
		}
	}



//	for (int i3 = 0; i3 < n3; i3++) {
//		for (int i2 = 0; i2 < n2; i2++) {
//			for (int i1 = 0; i1 < n1; i1++) {
//				getFeq(feq,1.0,0.0,0.0,0.0,theta_0,0.0);  
//				copyToGrid(i1,i2,i3,node,feq);
//			}
//		}
//	}
//	//cells
//	for (int i3 = 0; i3 < n3; i3++) {
//		for (int i2 = 0; i2 < n2; i2++) {
//			for (int i1 = 0; i1 < n1; i1++) {
//				getFeq(feq,1.0,0.0,0.0,0.0,theta_0,0.0); 
//				copyToGrid(i1,i2,i3,cell,feq);
//			}
//		}
//	}
}
#endif
