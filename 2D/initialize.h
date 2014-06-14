#ifndef _INITIALIZE_H
#define _INITIALIZE_H

#include<iostream>
#include<math.h>
#include<fstream>
#include "grid2D.h"
#include "globals.h"
#include "collide.h"

using namespace std;

//Debug vars
int count = 0;

void initialize(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<1> &nodeZero) {
	myReal u_mod, dot_prod;
	u_mod = pow(ux_init,2) + pow(uy_init,2);
	for (int i = 0 ; i < nx+2 ; i++) {
		for (int j = 0; j < ny +2 ; j++) {
			for (int k = 0; k < 4 ; k++) {
				dot_prod = c_x[k]*ux_init + c_y[k]*uy_init;
				nodeP(i,j,k) = w[k+1]*rho_init*(1.0 + 3.0*dot_prod + (9.0/2.0)*pow(dot_prod,2) - (3.0/2.0)*u_mod);
				nodeM(i,j,k) = w[k+5]*rho_init*(1.0 + 3.0*dot_prod + (9.0/2.0)*pow(dot_prod,2) - (3.0/2.0)*u_mod);
				nodeZero(i,j,0) = w[0]*rho_init*(1.0 + 3.0*dot_prod + (9.0/2.0)*pow(dot_prod,2) - (3.0/2.0)*u_mod);
//				nodeP(i,j,k) = ++count;
//				nodeM(i,j,k) = ++count;
//				nodeZero(i,j,k) = 1;
			}
		}
	}
}

void setNodeType(GridType &nodeType) {
	for (int i = 0; i < nx+2; i++) {
		for (int j = 0; j < ny+2; j++) {
			if (pow((i - center_x),2) + pow((j - center_y),2) < pow(radius,2)) 
				nodeType(i,j) = SOLID;
			else
				nodeType(i,j) = FLUID;
		}
	}
}

void setCellType(GridType &nodeType) {
	for (int i = 0; i < nx+2; i++) {
		for (int j = 0; j < ny+2; j++) {
			if (pow((i + 0.5 - center_x),2) + pow((j + 0.5- center_y),2) < pow(radius,2)) 
				nodeType(i,j) = SOLID;
			else
				nodeType(i,j) = FLUID;
		}
	}
}

void Print(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<1> &nodeZero, int o) {
		myReal rho, u_x, u_y, temp[9];
		rho = u_x = u_y = 0;
		myReal rho_global = 0;
		ofstream file;
		file.open("Velocity7.txt");
		for (int i = nx; i > 0; i--) {
			copyFromGrid(i, o, nodeP, nodeM, nodeZero, temp);
			getMoments(temp,rho,u_x,u_y);
			//print[i-1] = u_x;
			file<<u_x - g_x/2<<endl;
		}
		file.close();
//		for (int i = nx	; i > 0; i--) {
//			for (int j = 1; j < ny+1; j++) {
//				copyFromGrid(i, o, nodeP, nodeM, nodeZero, temp);
//				getMoments(temp,rho,u_x,u_y);
//				rho_global += rho;
//			}
//		}
//		cout<<rho_global;
}

void Print2(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<1> &nodeZero, int o) {
		myReal rho, u_x, u_y, temp[9];
		rho = u_x = u_y = 0;
		ofstream file;
		file.open("Velocity7.vtk");
		for (int i = 1; i < nx+1; i++) {
			for (int j = 1; j < ny+1; j++) {
				copyFromGrid(i, j, nodeP, nodeM, nodeZero, temp);
				getMoments(temp,rho,u_x,u_y);
				file<<u_x<<"  "<<u_y<<" "<<"0.0"<<endl;
				}
		}
		file.close();

//		ofstream file;
//		file.open("VelocityMoving.txt");
//		for (int i = nx; i > 0; i--) {
//			copyFromGrid(i, ny/2, nodeP, nodeM, nodeZero, temp);
//			getMoments(temp,rho,u_x,u_y);
//			//print[i-1] = u_x;
//			file<<i<<" "<<u_x - g_x/2<<endl;
//		}
//		file.close();

//		for (int i = nx; i > 0; i--) {
//			copyFromGrid(i, o, nodeP, nodeM, nodeZero, temp);
//			getMoments(temp,rho,u_x,u_y);
//			//print[i-1] = u_x;
//			cout<<u_x - g_x/2<<endl;
//		}
}



void PrintStream(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<1> &nodeZero, GridType nodeType) {
	myReal rho, u_x, u_y, temp[9];
		rho = u_x = u_y = 0;
		ofstream file;
		file.open("Streamlines.vtk");
		for (int i = 1; i < nx+1; i++) {
			for (int j = 1; j < ny+1; j++) {
				copyFromGrid(i, j, nodeP, nodeM, nodeZero, temp);
				getMoments(temp,rho,u_x,u_y);
//				if (nodeType(i,j) == SOLID) {
//					//cout<<"SOLID";
//					u_x = u_y = 0.0;
//				}
				file<<u_x<<" "<<u_y<<" "<<"0.0"<<endl;
				}
			//file<<endl;
		}
		file.close();
//	myReal rho, u_x, u_y, temp[9];
//	double theta;
//	rho = u_x = u_y = 0;
//	myReal rho_global = 0;
//	for (int i = nx; i > 0; i--) {
//		copyFromGrid(i, o, nodeP, nodeM, nodeZero, temp);
//		getMoments(temp,rho,u_x,u_y);
//		theta = atan(u_y/u_x) * 180/3.14;
//		cout<<theta<<endl;
//	}
}


void PrintReplica(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<1> &nodeZero, grid2D<N> &cellP, grid2D<N> &cellM, grid2D<1> &cellZero, GridType nodeType, int o) {
	myReal rho, u_x, u_y, temp[9];
	rho = u_x = u_y = 0;
	int j = 1;
	ofstream out;
	out.open("VelReplica.txt");
	out<<"0.0"<<endl;
	for (int i = 1; i < nx+1; i++) {
		copyFromGrid(i, o, cellP, cellM, cellZero, temp);
		getMoments(temp,rho,u_x,u_y);
		//cout<<u_x - g_x/2<<endl;
		out<<u_x - g_x/2<<endl;
//		copyFromGrid(i, o, nodeP, nodeM, nodeZero, temp);
//		getMoments(temp,rho,u_x,u_y);
//		cout<<u_x - g_x/2<<endl;
//		j++;
//		out<<u_x - g_x/2<<endl;
//		j++;
	}
	out.close();
}

void PrintReplica2(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<1> &nodeZero, grid2D<N> &cellP, grid2D<N> &cellM, grid2D<1> &cellZero, GridType &cellType, int o) {
	myReal rho, u_x, u_y, temp[9];
	rho = u_x = u_y = 0;
	//int j = 1;
	ofstream out;
	out.open("VelReplica.vtk");
	for (int i = 1; i < nx+1; i++) {
		for (int j = 1; j < ny+1; j++) {
			copyFromGrid(i, j, cellP, cellM, cellZero, temp);
			getMoments(temp,rho,u_x,u_y);
			if (cellType(i,j) == SOLID) {
					//cout<<"SOLID";
					u_x = u_y = 0.0;
				}
			out<<u_x<<" "<<u_y<<" "<<"0.0"<<endl;
			//out<<rho<<endl;
		}
	}
	out.close();
}



//void checkPrint(grid2D<N> &nodeP, grid2D<N> &nodeM) {
//	for (int i = nx+1 ; i >=0 ; i--) {
//		for (int j = 0; j < ny +2 ; j++) {
//			cout<< nodeP(i,j,3)<<" ";
//		}
//		cout<<endl;
//	}
//	//Print its reverse
//	cout<<endl;
//	for (int i = nx+1 ; i >=0 ; i--) {
//		for (int j = 0; j < ny +2 ; j++) {
//			cout<< nodeM(i,j,3)<<" ";
//		}
//		cout<<endl;
//	}
//}

//void checkPrintReplica(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<N> &cellP, grid2D<N> &cellM) {
//	for (int i = nx+1 ; i >=0 ; i--) {
//		for (int j = 0; j < ny +2 ; j++) {
//			cout<< nodeM(i,j,3)<<" ";
//		}
//		cout<<"      ";
//		for (int j = 0; j < ny +2 ; j++) {
//			cout<< cellM(i,j,3)<<" ";
//		}
//		cout<<endl;
//	}
//	//Print its reverse
//	cout<<endl;
//	for (int i = nx+1 ; i >=0 ; i--) {
//		for (int j = 0; j < ny +2 ; j++) {
//			cout<< cellM(i,j,1)<<" ";
//		}
//		cout<<"      ";
//		for (int j = 0; j < ny +2 ; j++) {
//			cout<< cellP(i,j,1)<<" ";
//		}
//		cout<<endl;
//	}

//	cout<<endl;
//	
//}

#endif		
