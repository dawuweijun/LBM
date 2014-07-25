#ifndef _INIT_3D_H
#define _INIT_3D_H
#include <iostream>
#include<cstdlib>
#include <math.h>
#include "gridV.h"
#include "globals3D.h"
#include "moments3D.h"
using namespace std;


void setTypeWall(gridtype nodeType, gridtype cellType) {
//	for (int i3 = 3; i3 < nodeType.m3-3; i3++) {
//		for (int i1 = 3; i1 < nodeType.m1-3; i1++) {
//			for (int index = cellType.iE2-1;  index < cellType.m2; index++) 
//				cellType(i1,index,i3) = nodeType(i1,index,i3) = SOLID; //Solid at the top-cells
//			for (int index = nodeType.iB2+1;  index >= 0; index--) 
//				cellType(i1,index,i3) = nodeType(i1,index,i3) = SOLID; //Solid at the bottom-nodes
//		}
//	}
	for (int i3 = 5; i3 < nodeType.m3-5; i3++) {
		for (int i1 = 15; i1 < nodeType.m1-15; i1++) {
				cellType(i1,NY/2,i3) = nodeType(i1,NY/2,i3) = SOLID; 
		}
	}
	ofstream file;
	file.open("type.txt");
	for (int i2 = 0; i2 < nodeType.m2; i2++) {
		for (int i1 = 0; i1 < nodeType.m1; i1++) {
			file<<nodeType(i1,i2,NZ/2);
		}
		file<<endl;
	}
	file.close();
}
void setTypeSphere(gridtype nodeType, gridtype cellType) {
	for (int i3 = 0; i3 < nodeType.m3; i3++) {
		for (int i2 = 0; i2 < nodeType.m2; i2++) {
			for (int i1 = 0; i1 < nodeType.m1; i1++) {
				if ((i1-center_x)*(i1-center_x) + (i2-center_y)*(i2-center_y) + (i3-center_z)*(i3-center_z) <= radius*radius) 
					nodeType(i1,i2,i3) = SOLID;
				else
					nodeType(i1,i2,i3) = FLUID;
			}
		}
	}
	for (int i3 = 0; i3 < cellType.m3; i3++) {
		for (int i2 = 0; i2 < cellType.m2; i2++) {
			for (int i1 = 0; i1 < cellType.m1; i1++) {
				if ((i1+0.5-center_x)*(i1+0.5-center_x) + (i2+0.5-center_y)*(i2+0.5-center_y) + (i3+0.5-center_z)*(i3+0.5-center_z) <= radius*radius)  
					cellType(i1,i2,i3) = SOLID;
				else
					cellType(i1,i2,i3) = FLUID;
			}
		}
	}
	ofstream file;
	file.open("type.txt");
	for (int i2 = 0; i2 < nodeType.m2; i2++) {
		for (int i1 = 0; i1 < nodeType.m1; i1++) {
			file<<nodeType(i1,i2,NZ/2);
		}
		file<<endl;
	}
	file.close();
	
}

void setTypeCylinder(gridtype nodeType, gridtype cellType) {
//	for (int i3 = 0; i3 < nodeType.m3; i3++) {
//		for (int i2 = 0; i2 < nodeType.m2; i2++) {
//			for (int i1 = 0; i1 < nodeType.m1; i1++) {
//				if ((i2 - center_y)*(i2 - center_y) + (i3-center_z)*(i3-center_z) <= rad*rad) 
//					nodeType(i1,i2,i3) = FLUID;
//				else
//					nodeType(i1,i2,i3) = SOLID;
//				
//			}
//		}
//	}
//	for (int i3 = 0; i3 < nodeType.m3; i3++) {
//		for (int i2 = 0; i2 < nodeType.m2; i2++) {
//			for (int i1 = 0; i1 < nodeType.m1; i1++) {
//				if ((i2 +0.5- center_y)*(i2 +0.5- center_y) + (i3 +0.5 - center_z)*(i3 +0.5 - center_z) <= rad*rad) 
//					cellType(i1,i2,i3) = FLUID;
//				else
//					cellType(i1,i2,i3) = SOLID;
//				
//			}
//		}
//	}
	for (int i3 = 0; i3 < nodeType.m3; i3++) {
		for (int i2 = 0; i2 < nodeType.m2; i2++) {
			for (int i1 = 0; i1 < nodeType.m1; i1++) {
				if (i1 > nodeType.iB1+4){// && i1 < nodeType.iE1-4) {
				if ((i2 - center_y)*(i2 - center_y) + (i3-center_z)*(i3-center_z) <= rad*rad &&(i2 - center_y)*(i2 - center_y) + (i3-center_z)*(i3-center_z) >= (rad-3)*(rad-3) ) 
					nodeType(i1,i2,i3) = SOLID;
				else
					nodeType(i1,i2,i3) = FLUID;
				}
				
			}
		}
	}
	for (int i3 = 0; i3 < nodeType.m3; i3++) {
		for (int i2 = 0; i2 < nodeType.m2; i2++) {
			for (int i1 = 0; i1 < nodeType.m1; i1++) {
				if (i1 > nodeType.iB1+4){// && i1 < nodeType.iE1-4) {
				if ((i2 +0.5- center_y)*(i2 +0.5- center_y) + (i3 +0.5 - center_z)*(i3 +0.5 - center_z) <= rad*rad && (i2 +0.5- center_y)*(i2 +0.5- center_y) + (i3 +0.5 - center_z)*(i3 +0.5 - center_z) >= (rad-3)*(rad-3)) 
					cellType(i1,i2,i3) = SOLID;
				else
					cellType(i1,i2,i3) = FLUID;
				}
			}
		}
	}
	ofstream file;
	file.open("type.txt");
	for (int i2 = 0; i2 < nodeType.m2; i2++) {
		for (int i1 = 0; i1 < nodeType.m1; i1++) {
			file<<nodeType(i1,i2,NZ/2);
		}
		file<<endl;
	}
	file.close();
//	ofstream file;
//	file.open("type.txt");
//	for (int i2 = 0; i2 < nodeType.m2; i2++) {
//		for (int i3 = 0; i3 < nodeType.m3; i3++) {
//			file<<nodeType(nodeType.iE1-10,i2,i3);
//		}
//		file<<endl;
//	}
//	file.close();
}

void setType(gridtype nodeType, gridtype cellType, gridtype nodeTypeOld, gridtype cellTypeOld) {
	#ifdef STAT_SPHERE
	 	setTypeSphere(nodeType,cellType);
	#endif
	#ifdef MOVING_SPHERE
		setTypeSphere(nodeType,cellType);
		setTypeSphere(nodeTypeOld,cellTypeOld);
	#endif
	#ifdef CYLINDER
		setTypeCylinder(nodeType,cellType);
	#endif
	#ifdef WALL
		setTypeWall(nodeType, cellType);
	#endif
}


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
