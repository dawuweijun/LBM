#include <iostream>
#include <math.h>
#include <unistd.h>
using namespace std;

#include "gridV.h"
#include "gridV_BC.h"
#include "globals3D.h"
#include "stream3D.h"
#include "boundary3D.h"
#include "initialize3D.h"
#include "moments3D.h"
#include "checkFunctions.h"

int main() {

	ofstream file,file2;
	file.open("Energy.txt");
	file2.open("Force.txt");
	int c = 0;


	lbgrid node(NX,NY,NZ,NX_G,NY_G,NZ_G);
	lbgrid cell(NX,NY,NZ,NX_G,NY_G,NZ_G);
	gridtype nodeType(NX,NY,NZ,NX_G,NY_G,NZ_G);
	gridtype cellType(NX,NY,NZ,NX_G,NY_G,NZ_G);

	initializeWithZeroVel(node,cell);
	//setTypeSphere(nodeType,cellType);
	setTypeWall(nodeType, cellType);
	//cout<<beta<<endl;
	for (int t = 0; t < SIMULATION_TIME; t++) {
//		cout<<"After periodic-----";
//		printGlobalRho(node,cell);
		collide(node);
		collide(cell);

		setInletOutlet(node,cell);
		periodic23(node);
		periodic23(cell);
		
		stream(node,cell);
		setInletOutletCorrection(node,cell);
		setInternalBoundary(node,cell,nodeType,cellType);
		if (t%100 == 0) {
			printGlobalRho(node,cell,nodeType,cellType,file);
			file2<<forceX/(u_inlet*u_inlet*radius*radius*3.14*0.5)<<endl;		
		}
		//printEachTime(node,t,c);
		//printVelX(node,nodeType,t,c);
		if(t%100==0) printAllVel(node,nodeType);// printEnergy(node,file);
		//cout<<t<<endl;
	}
	//printAllVel(node,nodeType);
	file.close();
	file2.close();

}
		
//*******To check if feq is working**********	

//	//CHECK FEQ
//	myReal f[27];	
//	myReal rho,ux,uy,uz,t,r;
//	rho=0.0871;ux=0.032543;uy=0.02112;uz=0.01326;t=1.0/12.0;r=0.0;
//	cout<<rho<<" "<<ux<<" "<<uy<<" "<<uz<<" "<<t<<" "<<endl;
//	getFeq(f,rho,ux,uy,uz,t);
//	//cout<<f[3]<<endl;
//	getMoments(f,rho,ux,uy,uz,t,r);
//	cout<<rho<<" "<<ux<<" "<<uy<<" "<<uz<<" "<<t<<" "<<endl;
//	//END CHECK FEQ	

		
	
