#ifndef _MOMENTS_3D_H
#define _MOMENTS_3D_H
#include <iostream>
#include<stdlib.h>
#include <math.h>
#include "gridV.h"
#include "globals3D.h"
using namespace std;

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
	node.Zero(i1,i2,i3,0) = f[index];
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

//Get moments
void getMoments(myReal f[27], myReal& rho, myReal& u_x, myReal& u_y, myReal& u_z, myReal& theta, myReal &R) {
	rho = u_x = u_y = u_z = theta = R = 0.0;
	for (int k = 0 ; k < 27 ; k++) {
		rho += f[k];
		u_x += f[k]*c_x[k];
		u_y += f[k]*c_y[k];
		u_z += f[k]*c_z[k];
		theta += f[k]*(c_x[k]*c_x[k] + c_y[k]*c_y[k] + c_z[k]*c_z[k]);
		//R += f[k]*(pow(c_x[k],2) + pow(c_y[k],2) + pow(c_z[k],2))*(pow(c_x[k],2) + pow(c_y[k],2) + pow(c_z[k],2));
	}
	u_x /= rho;
	u_y /= rho;
	u_z /= rho;
	u_x += g_x/2;
	//u_z += g_x/2;
	theta = (theta - rho*(u_x*u_x +u_y*u_y +u_z*u_z))/3.0;
//	if(theta==0) {cout<<"theta is zero";}
	theta /= rho;
	//R /= rho;
}

void getP(myReal f[27], myReal &pxx, myReal &pyy, myReal &pzz, myReal &pxy, myReal &pxz, myReal &pyz) {
	pxx = pyy = pzz = pxy = pxz = pyz = 0.0;
	for (int k = 0; k < 27; k++) {
		pxx += f[k]*c_x[k]*c_x[k];
		pyy += f[k]*c_y[k]*c_y[k];
		pzz += f[k]*c_z[k]*c_z[k];
		pxy += f[k]*c_x[k]*c_y[k];
		pxz += f[k]*c_x[k]*c_z[k];
		pyz += f[k]*c_y[k]*c_z[k];
	}
}
	

//Equilibrium
//void getFeq(myReal feq[27], myReal rho, myReal u_x, myReal u_y, myReal u_z, myReal theta) {
//	myReal h, f_tilda,c,u,dot_prod,R;
//	myReal deltheta = (theta - theta_0)/theta_0;
//	R = (3.0+(17.0*deltheta*0.25))*0.2 ;
//	h = (15.0*theta*theta - R)/(9.0*theta*theta - R);
//	u = u_x*u_x +u_y*u_y +u_z*u_z;
//	for (int k = 0; k < 27; k++) {
//		c = c_x[k]*c_x[k] + c_y[k]*c_y[k] + c_z[k]*c_z[k];
//		dot_prod = c_x[k]*u_x + c_y[k]*u_y + c_z[k]*u_z;
//		f_tilda = w[k]*rho*(1 + ((theta-theta_0)/(2.0*theta_0))*((c/theta_0) - 3.0));
//		feq[k] = f_tilda*(1 + dot_prod/theta + 0.5*(dot_prod/theta)*(dot_prod/theta) - (u/(2.0*theta))*(1.0-h) - h*u*c/(6.0*theta*theta));
//	}
//}
 inline void getFeq(myReal Feq[27], double rho, double ux, double uy, double uz, double theta){

  double theta0 = 1.0/5.0;
  double deltheta = (theta - theta0)/theta0;
  double theta2 = theta*theta;
  double theta0Inv = 1.0/theta0; 
  double threeByfour = 3.0/4.0;
  
  double fzero = 1.0 /3.0* rho *(1.0 + deltheta*0.5*(-3.0));
  
  double fSC =  1.0 /30.0* rho *(1.0 + deltheta*0.5*(theta0Inv -3.0));

  double fFCC = 1.0 /300.0* rho *(1.0 + deltheta*0.5*(2.0*theta0Inv -3.0));
  
  double fBCC = 4.0 /75.0* rho *(1.0 + deltheta*0.5*(theta0Inv*threeByfour -3.0));

  double uSq = ux*ux + uy*uy + uz*uz;

  double R = (3.0+(17.0*deltheta*0.25))*0.2 ;

  double h4 = (R-(15.0*theta2))/(R-(9.0*theta2)) ; 

  
  double factor2 = -uSq*0.5/theta * (1.0 - h4 ) ;
  
  Feq[DVG_ZERO_ZERO_ZERO]=  fzero*(1.0 + factor2);

  ux = ux/(theta), uy = uy/(theta), uz = uz/(theta);
  
  factor2 = -uSq*0.5/theta * (1.0 - h4 + h4/(3.0*theta) ) ;
  
  Feq[DVG_P1_ZERO_ZERO]	= fSC * (1.0 + ux + 0.5*ux*ux + factor2) ;
  Feq[DVG_M1_ZERO_ZERO]	= fSC * (1.0 - ux + 0.5*ux*ux + factor2); 

  Feq[DVG_ZERO_P1_ZERO]	= fSC * (1.0 + uy + 0.5*uy*uy + factor2); 
  Feq[DVG_ZERO_M1_ZERO]	= fSC * (1.0 - uy + 0.5*uy*uy + factor2); 

  Feq[DVG_ZERO_ZERO_P1]	= fSC * (1.0 + uz + 0.5*uz*uz + factor2); 
  Feq[DVG_ZERO_ZERO_M1]	= fSC * (1.0 - uz + 0.5*uz*uz + factor2); 
  
  //FCC 
  factor2 = -uSq*0.5/theta * (1.0 - h4 + (2.0*h4)/(3.0*theta) ) ;

    Feq[DVG_P1_P1_ZERO]= fFCC* (1.0 + ux+ uy  + 0.5*(ux +uy)*(ux+uy) + factor2); 
    Feq[DVG_M1_M1_ZERO]= fFCC* (1.0 - ux- uy + 0.5*(-ux -uy)*(-ux-uy) + factor2); 

    Feq[DVG_P1_M1_ZERO]=  fFCC* (1.0 + ux- uy + 0.5*(ux -uy)*(ux-uy) + factor2); 
    Feq[DVG_M1_P1_ZERO]=  fFCC* (1.0 - ux+ uy + 0.5*(-ux +uy)*(-ux+uy) + factor2); 

    Feq[DVG_P1_ZERO_P1]=  fFCC* (1.0 + ux+ uz + 0.5*(ux +uz)*(ux+uz) + factor2); 
    Feq[DVG_M1_ZERO_M1]=  fFCC* (1.0 - ux- uz + 0.5*(-ux -uz)*(-ux-uz) + factor2); 

    Feq[DVG_P1_ZERO_M1]= fFCC* (1.0 + ux- uz + 0.5*(ux -uz)*(ux-uz) + factor2); 
    Feq[DVG_M1_ZERO_P1]= fFCC* (1.0 - ux+ uz + 0.5*(-ux +uz)*(-ux+uz) + factor2); 

    Feq[DVG_ZERO_P1_P1]=fFCC* (1.0 + uy+ uz + 0.5*(uz +uy)*(uz+uy) + factor2); 
    Feq[DVG_ZERO_M1_M1]=fFCC* (1.0 - uy- uz + 0.5*(-uz -uy)*(-uz-uy) + factor2); 

    Feq[DVG_ZERO_P1_M1]=  fFCC* (1.0 + uy- uz + 0.5*(-uz +uy)*(-uz+uy) + factor2); 
    Feq[DVG_ZERO_M1_P1]= fFCC* (1.0 - uy+ uz + 0.5*(uz -uy)*(uz-uy) + factor2); 
    
    //BCC
    ux =ux/(2.0), uy = uy/(2.0), uz = uz/(2.0);
  
    factor2 = -uSq*0.5/theta * (1.0 - h4 + h4/(4.0*theta) ) ;

    Feq[DVG_P1_P1_P1]= fBCC* (1.0 +ux+uy+uz + 0.5*(ux +uy+uz)*(ux+uy+uz) + factor2) ;
    Feq[DVG_M1_M1_M1]= fBCC* (1.0 -ux-uy-uz + 0.5*(ux +uy+uz)*(ux+uy+uz) + factor2) ;

    Feq[DVG_P1_P1_M1]= fBCC* (1.0 +ux+uy-uz + 0.5*(ux +uy-uz)*(ux+uy-uz) + factor2) ;
    Feq[DVG_P1_M1_P1]=fBCC* (1.0 +ux-uy+uz + 0.5*(ux -uy+uz)*(ux-uy+uz) + factor2)  ;
    Feq[DVG_M1_P1_P1]=fBCC* (1.0 - ux+uy+uz + 0.5*(-ux +uy+uz)*(-ux+uy+uz) + factor2); 
    
    Feq[DVG_M1_M1_P1]= fBCC* (1.0 -ux-uy+uz + 0.5*(-ux -uy+uz)*(-ux-uy+uz) + factor2) ;
    Feq[DVG_M1_P1_M1]= fBCC* (1.0 -ux+uy-uz + 0.5*(-ux +uy-uz)*(-ux+uy-uz) + factor2) ;
    Feq[DVG_P1_M1_M1]= fBCC* (1.0 +ux-uy-uz + 0.5*(ux -uy-uz)*(ux-uy-uz) + factor2) ;
    
} 

void getFeqIndividual(myReal &feq, myReal rho, myReal u_x, myReal u_y, myReal u_z, myReal theta, int k) {
	myReal h, f_tilda,c,u,dot_prod,R;
	myReal deltheta = (theta - theta_0)/theta_0;
	R = (3.0+(17.0*deltheta*0.25))*0.2 ;
	h = (15.0*theta*theta - R)/(9.0*theta*theta - R);
	u = u_x*u_x +u_y*u_y +u_z*u_z;
	c = c_x[k]*c_x[k] + c_y[k]*c_y[k] + c_z[k]*c_z[k];
	dot_prod = c_x[k]*u_x + c_y[k]*u_y + c_z[k]*u_z;
	f_tilda = w[k]*rho*(1 + ((theta-theta_0)/(2.0*theta_0))*((c/theta_0) - 3.0));
	feq = f_tilda*(1 + dot_prod/theta + 0.5*(dot_prod/theta)*(dot_prod/theta) - (u/(2.0*theta))*(1.0-h) - h*u*c/(6.0*theta*theta));
}


//Siddhart's Grad function
void getGrad(myReal grad[27], myReal rho, myReal uX, myReal uY, myReal uZ, myReal pXX, myReal pYY, myReal pZZ, myReal pXY, myReal pXZ, myReal pYZ)   
{ 
  double dot;
  myReal w[27], cX[27], cY[27], cZ[27];
  
w[DVG_ZERO_ZERO_ZERO] = 1.0/3.0 ;
w[DVG_P1_ZERO_ZERO] = 1.0/30.0 ; w[DVG_M1_ZERO_ZERO] = 1.0/30.0; w[DVG_ZERO_P1_ZERO] = 1.0/30.0 ;
w[DVG_ZERO_M1_ZERO] = 1.0/30.0 ; w[DVG_ZERO_ZERO_P1] = 1.0/30.0; w[DVG_ZERO_ZERO_M1] = 1.0/30.0;
 
w[DVG_P1_P1_ZERO] = 1.0/300.0; w[DVG_M1_M1_ZERO] =  1.0/300.0; w[DVG_P1_M1_ZERO] =  1.0/300.0;
w[DVG_M1_P1_ZERO] = 1.0/300.0; w[DVG_P1_ZERO_P1] =  1.0/300.0; w[DVG_M1_ZERO_M1] =  1.0/300.0;
w[DVG_P1_ZERO_M1] = 1.0/300.0; w[DVG_M1_ZERO_P1] =  1.0/300.0; w[DVG_ZERO_P1_P1] =  1.0/300.0;
w[DVG_ZERO_M1_M1] =  1.0/300.0; w[DVG_ZERO_P1_M1] = 1.0/300.0; w[DVG_ZERO_M1_P1] =  1.0/300.0;


w[DVG_P1_P1_P1] = 4.0/75.0; w[DVG_M1_M1_M1] = 4.0/75.0; w[DVG_P1_P1_M1] = 4.0/75.0;
w[DVG_M1_M1_P1] = 4.0/75.0; w[DVG_P1_M1_P1] = 4.0/75.0; w[DVG_M1_P1_M1] = 4.0/75.0;
w[DVG_P1_M1_M1] = 4.0/75.0; w[DVG_M1_P1_P1] = 4.0/75.0; 

  
cX[DVG_ZERO_ZERO_ZERO]	= 0.0;
cX[DVG_P1_ZERO_ZERO]	= 1.0; cX[DVG_M1_ZERO_ZERO] = -1.0;	cX[DVG_ZERO_P1_ZERO] = 0.0;
cX[DVG_ZERO_M1_ZERO]	= 0.0 ; cX[DVG_ZERO_ZERO_P1] = 0.0; cX[DVG_ZERO_ZERO_M1] = 0.0; 

cX[DVG_P1_P1_ZERO]	= 1.0; cX[DVG_M1_M1_ZERO] = -1.0; cX[DVG_P1_M1_ZERO]	= 1.0;
cX[DVG_M1_P1_ZERO] = -1.0; cX[DVG_P1_ZERO_P1] = 1.0 ; cX[DVG_M1_ZERO_M1] = -1.0;
cX[DVG_P1_ZERO_M1] = 1.0; cX[DVG_M1_ZERO_P1] = -1.0; cX[DVG_ZERO_P1_P1]	= 0.0;
cX[DVG_ZERO_M1_M1] = 0.0 ; cX[DVG_ZERO_P1_M1] = 0.0; cX[DVG_ZERO_M1_P1] = 0.0 ;

cX[DVG_P1_P1_P1] =	0.5; cX[DVG_M1_M1_M1] =	-0.5; cX[DVG_P1_P1_M1] =  0.5;
cX[DVG_M1_M1_P1] =	-0.5 ; cX[DVG_P1_M1_P1] =	0.5 ; cX[DVG_M1_P1_M1] = -0.5 ;
cX[DVG_P1_M1_M1] =	0.5; cX[DVG_M1_P1_P1] =	-0.5 ;

cY[DVG_ZERO_ZERO_ZERO] =	0.0 ;
cY[DVG_P1_ZERO_ZERO] = 0.0; cY[DVG_M1_ZERO_ZERO] = 0.0; cY[DVG_ZERO_P1_ZERO] = 1.0 ;
cY[DVG_ZERO_M1_ZERO] = -1.0 ; cY[DVG_ZERO_ZERO_P1] =	0.0 ; cY[DVG_ZERO_ZERO_M1] =0.0; 
cY[DVG_P1_P1_ZERO] = 1.0; cY[DVG_M1_M1_ZERO] =	-1.0 ;	cY[DVG_P1_M1_ZERO] =-1.0 ;
cY[DVG_M1_P1_ZERO] = 1.0; cY[DVG_P1_ZERO_P1] =	0.0; cY[DVG_M1_ZERO_M1] =0.0 ;
cY[DVG_P1_ZERO_M1] =	0.0 ; cY[DVG_M1_ZERO_P1] = 0.0; cY[DVG_ZERO_P1_P1] =  1.0;
cY[DVG_ZERO_M1_M1] =	-1.0 ; cY[DVG_ZERO_P1_M1] = 1.0; cY[DVG_ZERO_M1_P1] = -1.0;

cY[DVG_P1_P1_P1] =	0.5; cY[DVG_M1_M1_M1] =	-0.5; cY[DVG_P1_P1_M1] =  0.5;
cY[DVG_M1_M1_P1] =	-0.5 ; cY[DVG_P1_M1_P1] =	-0.5 ; cY[DVG_M1_P1_M1] = 0.5 ;
cY[DVG_P1_M1_M1] =	-0.5; cY[DVG_M1_P1_P1] =	0.5; 

cZ[DVG_ZERO_ZERO_ZERO] = 0.0;
cZ[DVG_P1_ZERO_ZERO]	= 0.0;  cZ[DVG_M1_ZERO_ZERO] = 0.0 ;	cZ[DVG_ZERO_P1_ZERO] = 0.0 ;
cZ[DVG_ZERO_M1_ZERO] = 0.0; cZ[DVG_ZERO_ZERO_P1] = 1.0;	cZ[DVG_ZERO_ZERO_M1] = -1.0; 

cZ[DVG_P1_P1_ZERO] = 0.0; cZ[DVG_M1_M1_ZERO] = 0.0; cZ[DVG_P1_M1_ZERO] = 0.0;
cZ[DVG_M1_P1_ZERO] = 0.0; cZ[DVG_P1_ZERO_P1] = 1.0;	cZ[DVG_M1_ZERO_M1] = -1.0;
cZ[DVG_P1_ZERO_M1] = -1.0; cZ[DVG_M1_ZERO_P1] = 1.0; cZ[DVG_ZERO_P1_P1] = 1.0;
cZ[DVG_ZERO_M1_M1] = -1.0; cZ[DVG_ZERO_P1_M1] = -1.0; cZ[DVG_ZERO_M1_P1] = 1.0;

cZ[DVG_P1_P1_P1] = 0.5 ; cZ[DVG_M1_M1_M1] = -0.5; cZ[DVG_P1_P1_M1] = -0.5;
cZ[DVG_M1_M1_P1] = 0.5; cZ[DVG_P1_M1_P1] = 0.5 ; cZ[DVG_M1_P1_M1] = -0.5;
cZ[DVG_P1_M1_M1] = -0.5; cZ[DVG_M1_P1_P1] = 0.5;


for(int dv=0;dv<27;dv++)
{
dot = cX[dv] *uX + cY[dv] * uY + cZ[dv]*uZ;
grad[dv] = w[dv] *(rho + 5.0*rho*dot + (25.0/2.0)*((pXX -rho/5.0)*(cX[dv]*cX[dv] - 1.0/5.0) + 2.0*cX[dv]*cY[dv]*pXY
+ 2.0*cX[dv]*cZ[dv]*pXZ + 2.0*cY[dv]*cZ[dv]*pYZ  + (cY[dv]*cY[dv] - 1.0/5.0)*(pYY-rho/5.0) + (pZZ -rho/5.0)*(cZ[dv]*cZ[dv] - 1.0/5.0)));
}

}

void collide(lbgrid& node) {
	myReal temp[27], feq[27], rho, ux, uy, uz, theta, R;	
	for (int i3 = node.SCM.iB3; i3 <= node.SCM.iE3; i3++) {
		for (int i2 = node.SCM.iB2; i2 <= node.SCM.iE2; i2++) {
			for (int i1 = node.SCM.iB1; i1 <= node.SCM.iE1; i1++) {
				copyFromGrid(i1,i2,i3,node,temp);
				getMoments(temp,rho,ux,uy,uz,theta,R);
				//if (rho>10) {cout<<"Rho is large"<<i1<<" "<<i2<<" "<<i3<<endl; exit(0);}
				getFeq(feq,rho,ux,uy,uz,theta);
				for (int k = 0; k < 27; k++) {
				//if (temp[k] < 0) { cout<<k<<"AT"<<i1<<" "<<i2<<" "<<i3<<endl; exit(0);}
					 temp[k] += 2.0*beta*(feq[k] - temp[k]) + (2.0*beta*tau*w[k]*rho*c_x[k]*g_x)/theta_0;
				}
				copyToGrid(i1,i2,i3,node,temp);
			}
		}
	}
}



#endif
