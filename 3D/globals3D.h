#ifndef _GLOBALS_3D_H
#define _GLOBALS_3D_H
const int N = 4;
const int NX = 16;
const int NY = 16;
const int NZ = 16;
const int NX_G = 1;  //Number of ghost nodes in any direction
const int NY_G = 1;
const int NZ_G = 1;

const int SIMULATION_TIME = 15000;


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
//myReal w[27];

////Simple grid
//myReal c_x[27] = {0,1,0,0,-1,0,0,1,-1,1,-1,1,1,-1,-1,0,0,0,0,1,-1,1,-1,-1,1,-1,1};
//myReal c_y[27] = {0,0,1,0,0,-1,0,1,1,-1,-1,0,0,0,0,1,-1,1,-1,1,1,-1,-1,-1,-1,1,1};
//myReal c_z[27] = {0,0,0,1,0,0,-1,0,0,0,0,1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,-1,-1};


double radius = 4.5;
int center_x = NX/2;
int center_y = NY/2;
int center_z = NZ/2;

double u_inlet = 0.04;
double Re = 10.0;
double nu = u_inlet*NY/(Re);
double g_x = (4.0*u_inlet*u_inlet)/(Re*NY);
double theta_0 = 1.0/5.0;

double tau = nu/theta_0;
double dt = 1;//6.28/NX;
double beta = dt/(2*tau + dt);


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
//	w[0] = 1.0/3.0;
//for (int i = 1; i<6; i++)
//	w[i] = 1.0/30.0;
//for (int i = 6; i<18; i++)
//	w[i] = 1.0/300.0;
//for (int i = 18; i<26; i++)
//	w[i] = 4.0/75.0;
#endif
