#include <iostream>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

#include "grid2D.h"
#include "collide.h"
#include "stream.h"
#include "boundary.h"
#include "initialize.h"
#include "globals.h"
#include "moments.h"
//#include "moments.h"

using namespace std;
int c = 0;
int main() {
	#ifdef GRID_STAGGERED
	grid2D<4> lbgridP(nx,ny);
	grid2D<4> lbgridM(nx,ny);
	grid2D<1> lbgridZero(nx,ny);
	grid2D<4> lbgridReplicaP(nx,ny);
	grid2D<4> lbgridReplicaM(nx,ny);
	grid2D<1> lbgridReplicaZero(nx,ny);
	GridType lbgridType(nx,ny);
	GridType lbgridReplicaType(nx,ny);
	initialize(lbgridP,lbgridM,lbgridZero);
	initialize(lbgridReplicaP, lbgridReplicaM, lbgridReplicaZero);
	setNodeType(lbgridType);
	setCellType(lbgridReplicaType);
	#else
	grid2D<4> lbgridP(nx,ny);
	grid2D<4> lbgridM(nx,ny);
	grid2D<1> lbgridZero(nx,ny);
	GridType lbgridType(nx,ny);
	initialize(lbgridP,lbgridM,lbgridZero);
	setNodeType(lbgridType);
	#endif
//	ofstream f;
//	f.open("type.txt");
//	for (int i = 1; i < nx+1; i++) {
//		for (int j = 1; j < ny+1; j++) {
//			f<<lbgridType(i,j);
//		}
//		f<<endl;
//	}
//	for (int i = 1; i < nx+1; i++) {
//		for (int j = 1; j < ny+1; j++) {
//			f<<lbgridReplicaType(i,j);
//		}
//		f<<endl;
//	}
//	f.close();

	for (int t = 0; t < SIMULATION_TIME; t++) {

		#ifdef GRID_STAGGERED
			prepareBC(lbgridP);
			prepareBC(lbgridM);
			prepareBC(lbgridReplicaP);
			prepareBC(lbgridReplicaM);
			#ifdef PERIODIC
				periodicY(lbgridP);
				periodicY(lbgridM);
				periodicY(lbgridReplicaP);
				periodicY(lbgridReplicaM);	
			#endif
			#ifdef INLET_OUTLET
				setInletOutletReplica(lbgridP, lbgridM, lbgridZero, lbgridReplicaP, lbgridReplicaM, lbgridReplicaZero);
			#endif

			collideReplica(lbgridP,lbgridM,lbgridZero,lbgridReplicaP, lbgridReplicaM, lbgridReplicaZero);
			streamReplica(lbgridP, lbgridM, lbgridReplicaP, lbgridReplicaM);

			#ifdef BOUNCE_BACK
				applyBounceBackReplica(lbgridP, lbgridM, lbgridReplicaP, lbgridReplicaM);
			#endif
			#ifdef INLET_OUTLET
			setInletCorrectionReplica(lbgridP, lbgridM, lbgridReplicaP, lbgridReplicaM);
			#endif
			#ifdef MOVING_WALL
				applyMovingBoundaryReplica(lbgridP, lbgridM, lbgridZero, lbgridReplicaP, lbgridReplicaM, lbgridReplicaZero);
			#endif

			#ifdef OBJECT_IN_FLOW
				setObjectBCReplica(lbgridP, lbgridM, lbgridZero, lbgridReplicaP, lbgridReplicaM, lbgridReplicaZero, lbgridType, lbgridReplicaType);
			#endif


		#else
			prepareBC(lbgridP);
			prepareBC(lbgridM);
			#ifdef PERIODIC
				periodicY(lbgridP);
				periodicY(lbgridM);
				periodicY(lbgridZero);
			#endif
			#ifdef INLET_OUTLET
				setInletOutlet(lbgridP, lbgridM, lbgridZero);
			#endif

			collide(lbgridP,lbgridM, lbgridZero);
			stream(lbgridP, lbgridM);

			#ifdef BOUNCE_BACK
				applyBounceBack(lbgridP, lbgridM);
			#endif
			#ifdef MOVING_WALL
				applyMovingBoundary(lbgridP, lbgridM, lbgridZero);
			#endif
			
			#ifdef INLET_OUTLET
			setInletCorrection(lbgridP, lbgridM);
			#endif 

			#ifdef OBJECT_IN_FLOW
				setObjectBC(lbgridP, lbgridM, lbgridZero, lbgridType);
			#endif
		#endif

//		if (t%100 == 0) {
//			Print(lbgridP, lbgridM, lbgridZero,ny/2); 
//			cout<<endl;
//			//sleep(2);
//		}
		//checkPrint(lbgridP, lbgridM);
//		if (t%100 == 0){ 
//			myReal rho, u_x, u_y, temp[9];
//			rho = u_x = u_y = 0;
//			ofstream file;
//			char name[20];
//			sprintf(name,"CHeckVelocity%d.txt",++c);
//			file.open(name);
//			for (int i = nx; i > 0; i--) {
//				for (int j = 1; j < ny+1; j++) {
//					copyFromGrid(i, j, lbgridP, lbgridM, lbgridZero, temp);
//					getMoments(temp,rho,u_x,u_y);
//					file<<u_x<<"  ";
//				}
//				file<<endl;
//			}
//			file.close();
//		}

	}

		
	//Print(lbgridP, lbgridM, lbgridZero, ny/2);
	//Print2(lbgridP, lbgridM, lbgridZero,ny/2);
	//PrintStream(lbgridP, lbgridM, lbgridZero, lbgridType);
	PrintReplica2(lbgridP, lbgridM, lbgridZero,lbgridReplicaP, lbgridReplicaM, lbgridReplicaZero, lbgridReplicaType, ny/2);
		
	
	//cout<<tau;
}
	
