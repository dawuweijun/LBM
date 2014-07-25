#ifndef _GLOBALS_3D_H
#define _GLOBALS_3D_H
#include <math.h> 

//#define STAT_SPHERE
//#define INLET
//#define MOVING_SPHERE
#define CYLINDER
//#define WALL

const int N = 4;
const int NX = 200;
const int NY = 150;
const int NZ = 150;
const int NX_G = 1;  //Number of ghost nodes in any direction
const int NY_G = 1;
const int NZ_G = 1;

const int SIMULATION_TIME = 20001;//40000;

int b = 0;

struct lbgrid {
	gridV<N> SCP;
	gridV<N> SCM; 
	gridV<N> FCC12; 
	gridV<N> FCC13;
	gridV<N> FCC23; 
	gridV<N> BCCP;
	gridV<N> BCCM; 
	gridV<1> Zero;
	//Constructor initializer list
	lbgrid(int nx, int ny, int nz, int nx_g, int ny_g, int nz_g) : SCP(nx,ny,nz,nx_g,ny_g,nz_g),SCM(nx,ny,nz,nx_g,ny_g,nz_g),FCC12(nx,ny,nz,nx_g,ny_g,nz_g),FCC13(nx,ny,nz,nx_g,ny_g,nz_g),FCC23(nx,ny,nz,nx_g,ny_g,nz_g),BCCP(nx,ny,nz,nx_g,ny_g,nz_g),BCCM(nx,ny,nz,nx_g,ny_g,nz_g),Zero(nx,ny,nz,nx_g,ny_g,nz_g) {}

};

enum type{SOLID,FLUID}; //Node type

enum SCP {         //The various distributions: P1 - 1; M1 - -1	  
     	DV_P1_ZERO_ZERO,
	DV_ZERO_P1_ZERO,
	DV_ZERO_ZERO_P1,
};
enum SCM {
	DV_M1_ZERO_ZERO,
     	DV_ZERO_M1_ZERO,
     	DV_ZERO_ZERO_M1,
};
enum FCC12 {
	DV_P1_P1_ZERO,
	DV_M1_P1_ZERO,
	DV_P1_M1_ZERO,
	DV_M1_M1_ZERO,
};
enum FCC13 {
	DV_P1_ZERO_P1,
	DV_P1_ZERO_M1,
	DV_M1_ZERO_P1,
	DV_M1_ZERO_M1,
};
enum FCC23 {
	DV_ZERO_P1_P1,
	DV_ZERO_M1_P1,
	DV_ZERO_P1_M1,
	DV_ZERO_M1_M1,
};
enum BCCP {
	DV_P1_P1_P1,
	DV_M1_P1_P1,
	DV_P1_M1_P1,
	DV_M1_M1_P1,
};
enum BCCM {
	DV_M1_M1_M1,
	DV_P1_M1_M1,
	DV_M1_P1_M1,
	DV_P1_P1_M1,
};

enum ZERO { DV_ZERO_ZERO_ZERO };

myReal c_x[27] = {0,1,0,0,-1,0,0,1,-1,1,-1,1,1,-1,-1,0,0,0,0,0.5,-0.5,0.5,-0.5,-0.5,0.5,-0.5,0.5};
myReal c_y[27] = {0,0,1,0,0,-1,0,1,1,-1,-1,0,0,0,0,1,-1,1,-1,0.5,0.5,-0.5,-0.5,-0.5,-0.5,0.5,0.5};
myReal c_z[27] = {0,0,0,1,0,0,-1,0,0,0,0,1,-1,1,-1,1,1,-1,-1,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5};
myReal w[27] = {1.0/3.0,1.0/30.0,1.0/30.0,1.0/30.0,1.0/30.0,1.0/30.0,1.0/30.0,1.0/300.0,1.0/300.0,1.0/300.0,1.0/300.0,1.0/300.0,1.0/300.0,1.0/300.0,1.0/300.0,1.0/300.0,1.0/300.0,1.0/300.0,1.0/300.0,4.0/75.0,4.0/75.0,4.0/75.0,4.0/75.0,4.0/75.0,4.0/75.0,4.0/75.0,4.0/75.0};
//myReal w[27];(nodeTypeOld, nodeType, cellTypeOld, cellType);
// 		setTypeSphe

////Simple grid
//myReal c_x[27] = {0,1,0,0,-1,0,0,1,-1,1,-1,1,1,-1,-1,0,0,0,0,1,-1,1,-1,-1,1,-1,1};
//myReal c_y[27] = {0,0,1,0,0,-1,0,1,1,-1,-1,0,0,0,0,1,-1,1,-1,1,1,-1,-1,-1,-1,1,1};
//myReal c_z[27] = {0,0,0,1,0,0,-1,0,0,0,0,1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,-1,-1};


double radius = 8.0;
int rad = min(NZ/2,NY/2) - 4;
double center_x = NX/2;
double center_y = NY/2;
double center_z = NZ/2;
double kn = 0.01;
double Re = 0.1;
double u_inlet = kn*Re/sqrt(5);


double g_x = 0.0;//(4.0*u_inlet*u_inlet)/(Re*2*rad);
double theta_0 = 1.0/5.0;

double nu = u_inlet*(2*radius)/(Re);
double tau = nu/theta_0;
double dt = 1.0; //6.28/NX;
double beta = dt/(2*tau + dt);

myReal forceXN = 0.0;
myReal forceXC = 0.0;
myReal forceXOld = 0.0;
myReal forceYN = 0.0;
myReal forceYC = 0.0;
myReal forceYOld = 0.0;
myReal forceZN = 0.0;
myReal forceZC = 0.0;
myReal forceZOld = 0.0;
myReal torqueX = 0.0;
myReal torqueXOld = 0.0;
myReal torqueY = 0.0;
myReal torqueYOld = 0.0;
myReal torqueZ = 0.0;
myReal torqueZOld = 0.0;

myReal u_sph_x = u_inlet;
myReal u_sph_y = 0.0;
myReal u_sph_z = 0.0;
myReal w_sph_x = 0.0;
myReal w_sph_y = 0.0;
myReal w_sph_z = 0.0;

myReal u_ref;
myReal w_ref;

//Final vel
myReal v_sph_x = 0.0;
myReal v_sph_y = 0.0;
myReal v_sph_z = 0.0;

//Temp vars: No physical significance
myReal ftempX = 0.0;
myReal ftempY = 0.0;
myReal ftempZ = 0.0;



void swap(int &a, int &b) {
	int temp;
	temp = a;
	a = b;
	b = temp;
}
void swap(myReal &a, myReal &b) {
	myReal temp;
	temp = a;
	a = b;
	b = temp;
}

myReal avg(myReal a, myReal b) {
	return ((a+b)/2.0);
}

myReal distance(int i1, int i2, int i3) {
	return (sqrt((i1-center_x)*(i1-center_x) + (i2-center_y)*(i2-center_y) + (i3-center_z)*(i3-center_z)));
}

int min(int a, int b) {
	if (a < b) return a;
	else return b;
}

enum DIST { 
	DVG_ZERO_ZERO_ZERO,        //The various distributions: P1 - 1; M1 - -1	  
     	DVG_P1_ZERO_ZERO,
	DVG_ZERO_P1_ZERO,
	DVG_ZERO_ZERO_P1,

	DVG_M1_ZERO_ZERO,
     	DVG_ZERO_M1_ZERO,
     	DVG_ZERO_ZERO_M1,

	DVG_P1_P1_ZERO,
	DVG_M1_P1_ZERO,
	DVG_P1_M1_ZERO,
	DVG_M1_M1_ZERO,

	DVG_P1_ZERO_P1,
	DVG_P1_ZERO_M1,
	DVG_M1_ZERO_P1,
	DVG_M1_ZERO_M1,

	DVG_ZERO_P1_P1,
	DVG_ZERO_M1_P1,
	DVG_ZERO_P1_M1,
	DVG_ZERO_M1_M1,

	DVG_P1_P1_P1,
	DVG_M1_P1_P1,
	DVG_P1_M1_P1,
	DVG_M1_M1_P1,

	DVG_M1_M1_M1,
	DVG_P1_M1_M1,
	DVG_M1_P1_M1,
	DVG_P1_P1_M1,
};
#endif
