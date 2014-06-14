#ifndef _GLOBALS_3D_H
#define _GLOBALS_3D_H
const int N = 4;
const int NX = 20;
const int NY = 20;
const int NZ = 20;
const int NX_G = 1;  //Number of ghost nodes 
const int NY_G = 1;
const int NZ_G = 1;


struct lbgrid {
	gridV<N-1> SCP;
	gridV<N-1> SCM; 
	gridV<N> FCC12; 
	gridV<N> FCC13;
	gridV<N> FCC23; 
	gridV<N> BCCP;
	gridV<N> BCCM; 
	gridV<1> Zero;
	//Constructor initializer list
	lbgrid(int nx, int ny, int nz, int nx_g, int ny_g, int nz_g) : SCP(nx,ny,nz,nx_g,ny_g,nz_g),SCM(nx,ny,nz,nx_g,ny_g,nz_g),FCC12(nx,ny,nz,nx_g,ny_g,nz_g),FCC13(nx,ny,nz,nx_g,ny_g,nz_g),FCC23(nx,ny,nz,nx_g,ny_g,nz_g),BCCP(nx,ny,nz,nx_g,ny_g,nz_g),BCCM(nx,ny,nz,nx_g,ny_g,nz_g),Zero(nx,ny,nz,nx_g,ny_g,nz_g) {}

};

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
	DV_P1_P1_M1,
	DV_M1_P1_M1,
	DV_P1_M1_M1,
	DV_M1_M1_M1,
};

enum ZERO { DV_ZERO_ZERO_ZERO };

myReal c_x[27] = {0,1,0,0,-1,0,0,1,-1,1,-1,1,1,-1,-1,0,0,0,0,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5};
myReal c_y[27] = {0,0,1,0,0,-1,0,1,1,-1,-1,0,0,0,0,1,-1,1,-1,0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5,-0.5};
myReal c_z[27] = {0,0,0,1,0,0,-1,0,0,0,0,1,-1,1,-1,1,1,-1,-1,0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5};
myReal w[27] = {1/3,1/30,1/30,1/30,1/30,1/30,1/30,1/300,1/300,1/300,1/300,1/300,1/300,1/300,1/300,4/75,4/75,4/75,4/75,4/75,4/75,4/75,4/75};

////Simple grid
//myReal c_x[27] = {0,1,0,0,-1,0,0,1,-1,1,-1,1,1,-1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1};
//myReal c_y[27] = {0,0,1,0,0,-1,0,1,1,-1,-1,0,0,0,0,1,-1,1,-1,1,1,-1,-1,1,1,-1,-1};
//myReal c_z[27] = {0,0,0,1,0,0,-1,0,0,0,0,1,-1,1,-1,1,1,-1,-1,1,1,1,1,-1,-1,-1,-1};

myReal g_x = 0.0;
	
myReal tau = 1;
myReal beta = 1/(2*tau + 1);
myReal theta_0 = 1/6;



void swap(int a, int b) {
	int temp;
	temp = a;
	a = b;
	b = temp;
}
#endif
