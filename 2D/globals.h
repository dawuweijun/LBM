#ifndef _GLOBALS_H
#define _GLOBALS_H

#include<math.h>

#define GRID_STAGGERED
//#define MOVING_WALL
#define BOUNCE_BACK
//#define PERIODIC
#define INLET_OUTLET
#define OBJECT_IN_FLOW


const int N = 4;
const int nx = 50;   //no of rows
const int ny = 300;

//const double Ma = 0.03;
const double Re = 60.0;

const double u_inlet = 0.04;  
const double g_x = 0.0;//(4.0*u_inlet*u_inlet)/(Re*nx);
const double u_boundary = 0.0;
const double radius = 4.5;

const int center_x = nx/2;
const int center_y = ny/2;

const int SIMULATION_TIME = 5000;

myReal rho_init = 1.0;
myReal ux_init = 0.0;
myReal uy_init = 0.0;

//myReal print[nx];

#ifdef GRID_STAGGERED
double c_x[9] = {0,1,0.5,0,-0.5,-1,-0.5,0,0.5};
double c_y[9] = {0,0,0.5,1,0.5,0,-0.5,-1,-0.5};
myReal w[9]= {4.0/9.0, 1.0/36.0,1.0/9.0,1.0/36.0,1.0/9.0,1.0/36.0,1.0/9.0,1.0/36.0,1.0/9.0};
const double theta = 1.0/6.0;

#else
int c_x[9] = {0,1,1,0,-1,-1,-1,0,1};
int c_y[9] = {0,0,1,1,1,0,-1,-1,-1};
myReal w[9]= {4.0/9.0, 1.0/9.0,1.0/36.0,1.0/9.0,1.0/36.0,1.0/9.0,1.0/36.0,1.0/9.0,1.0/36.0};
const double theta = 1.0/3.0;

#endif

double nu = u_inlet*(nx-1)/(Re);
double tau = nu/theta;
//double tau = 1.4;//g_x/(theta*sqrt(theta)) *(1/Ma);
double beta = 1/(2*tau + 1);
//double Kn = Ma/Re;
//double tau = Kn*nx/sqrt(theta) ;

////Is this required???--No:)
//myReal c_x1[5] = {0,1,0,0,-1};
//myReal c_y1[5] = {0,0,0,1,0};

#endif
