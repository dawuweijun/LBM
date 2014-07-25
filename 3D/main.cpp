//MAIN
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
	
	//CHECK BETA OR KNUDSEN
	cout<<"BETA:"<<beta<<endl;
	//END CHECK

	ofstream fEnergy, fCenter, fForce, fU, fOmega;
	fEnergy.open("Energy.txt");
	fCenter.open("Center_loc.txt");
	fForce.open("Force.txt");
	fU.open("Trans_vel.txt");
	fOmega.open("Rot_vel.txt");
	int c = 0;

	//REFERENCE VELOCITIES
	u_ref = u_inlet;
	w_ref = w_sph_z;
	//END REF VEL

	lbgrid node(NX,NY,NZ,NX_G,NY_G,NZ_G);
	lbgrid cell(NX,NY,NZ,NX_G,NY_G,NZ_G);
	gridtype nodeType(NX,NY,NZ,NX_G,NY_G,NZ_G);
	gridtype cellType(NX,NY,NZ,NX_G,NY_G,NZ_G);

	gridtype nodeTypeOld(NX,NY,NZ,NX_G,NY_G,NZ_G);
	gridtype cellTypeOld(NX,NY,NZ,NX_G,NY_G,NZ_G);

	arrayRT rhoNode(NX,NY,NZ,NX_G,NY_G,NZ_G);
	arrayRT thetaNode(NX,NY,NZ,NX_G,NY_G,NZ_G);
	arrayRT rhoCell(NX,NY,NZ,NX_G,NY_G,NZ_G);
	arrayRT thetaCell(NX,NY,NZ,NX_G,NY_G,NZ_G);

	//Initialize
	initializeWithZeroVel(node,cell);
	setType(nodeType,cellType,nodeTypeOld,cellTypeOld);
	//End initialize
		
	//First print-check
	printStokesCorrection();
	printU(fU);
	printOmega(fOmega);
	//End check

	//START SIMULATION
	for (int t = 0; t < SIMULATION_TIME; t++) {
		//refill(nodeTypeOld, cellTypeOld, nodeType, cellType, node, cell);
		checkInit(node,nodeType);

		collide(node);
		collide(cell);
		#ifdef INLET
			setInletOutlet(node,cell);
			periodic23(node);
			periodic23(cell);	
		#else
			periodic123(node);
			periodic123(cell);
		#endif

		stream(node,cell);

//		forceXOld = forceX;
//		torqueZOld = torqueZ;
//		forceYOld = forceY;
//		torqueYOld = torqueY;
//		forceZOld = forceZ;
//		torqueXOld = torqueX;

		setInternalBoundary(node, cell, nodeType, cellType, rhoNode, thetaNode, rhoCell, thetaCell);
//		setBoundaryCorrection(node,cell,nodeType,cellType);
		#ifdef INLET 
			setInletOutletCorrection(node,cell);
		#endif

		#ifdef MOVING_SPHERE
			updateCMObj(nodeTypeOld, nodeType, cellTypeOld, cellType);
			setTypeSphere(nodeType, cellType);
		#endif
		
		//CHECK CENTER
		if (center_x > NX - 2*radius || center_y > NY - 2*radius || center_z >NZ - 2*radius || center_x < 2*radius || center_y < 2*radius || center_z < 2*radius) break;
		//END CENTER CHECK


		//******PRINT ROUTINES*********
		if (t%10 == 0) {
//			printCenterLocation(fCenter);
			printForceAndTorque(fForce);			
			printGlobalRho(node,cell,nodeType,cellType,fEnergy);
//			printAllVel(node,nodeType);
//			printU(fU);
//			printOmega(fOmega);			
//			printVelX(node,nodeType,t,c);
		}
	
		//*******END PRINT ***********


	}
	//END SIMULATION
	//close files
	fCenter.close();
	fEnergy.close();
	fForce.close();
	fU.close();
	fOmega.close();
}

//END MAIN
		
	
