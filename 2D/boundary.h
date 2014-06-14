#ifndef _BOUNDARY_H
#define _BOUNDARY_H

#include<iostream>
#include<math.h>
#include "grid2D.h"
#include "globals.h"
#include "moments.h"

using namespace std;

void applyBounceBack(grid2D<N> &nodeP, grid2D<N> &nodeM) {
	//top and bottom
	for (int j = 1; j < ny+1; j++) {
		for (int k = 1; k < N; k++) {
			nodeM(nx,j,k) = nodeP(nx+1,j,k);
			nodeP(1,j,k) = nodeM(0,j,k);
		}
	}
	//sides
//	for (int i = 1; i < nx+1; i++) {
//			nodeP(i,ny,3) = nodeM(i,ny+1,3);
//			nodeM(i,ny,0) = nodeP(i,ny+1,0);
//			nodeM(i,ny,1) = nodeP(i,ny+1,1);
//			nodeP(i,1,0) = nodeM(i,0,0);
//			nodeP(i,1,1) = nodeM(i,0,1);
//			nodeM(i,1,3) = nodeP(i,0,3);
//	}
	
}

void setInletOutlet(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<1> &nodeZero) {
	myReal f_grad[9], temp[9], feq[9], rho, ux, uy;
	//Inlet
	for (int i = 0; i < nx+2; i++) {
		copyFromGrid(i, 1, nodeP, nodeM, nodeZero, temp);
		getMoments(temp,rho,ux,uy);
		getFeq(feq,rho,u_inlet,0.0);
		//getGrad(feq,temp,1.3,u_inlet,0.0);
		copyToGrid(i, 0, nodeP, nodeM, nodeZero, feq);
	}
	//Outlet
	for (int i = 0; i < nx+2; i++) {
		copyFromGrid(i, ny, nodeP, nodeM, nodeZero, temp);
		getMoments(temp,rho,ux,uy);
		getFeq(feq,rho,ux,uy);
		//getGrad(feq,temp,rho,ux,uy);
		copyToGrid(i, ny+1, nodeP, nodeM, nodeZero, feq);
	}
}

void setInletOutletReplica(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<1> &nodeZero, grid2D<N> &cellP, grid2D<N> &cellM, grid2D<1> &cellZero) {
	myReal f_grad[9], temp[9],  feq[9], rho, ux, uy, avg;
	//Inlet
	for (int i = 1; i < nx+1; i++) {
		copyFromGrid(i, 1, nodeP, nodeM, nodeZero, temp);
		getMoments(temp,rho,ux,uy);
		getFeq(feq,rho,u_inlet,0.0);
		//getFeq(feq,1.003,ux,uy);
		//getGrad(feq,temp,1.024,u_inlet,0.0);
		copyToGrid(i, 0, nodeP, nodeM, nodeZero, feq);
		copyToGrid(i, 0, cellP, cellM, cellZero, feq);
	}
	for (int i = 1; i < nx+1; i++) {
		avg = (nodeP(i,1,0) + nodeP(i+1,1,0))/2;
		cellP(i,0,0) = avg; 
	}
	//Outlet
	for (int i = 1; i < nx+1; i++) {
		copyFromGrid(i, ny, cellP, cellM, cellZero, temp);
		getMoments(temp,rho,ux,uy);
		getFeq(feq,rho,ux,uy);
		//getGrad(feq,temp,rho,ux,uy);
		copyToGrid(i, ny+1, cellP, cellM, cellZero, feq);
		copyToGrid(i, ny+1, nodeP, nodeM, nodeZero, feq);
	}
	for (int i = 1; i < nx+1; i++) {
		avg = (cellM(i,ny,0) + cellM(i-1,ny,0))/2;
		nodeM(i,ny+1,0) = avg;
	}
}


void setInletCorrection(grid2D<N> &nodeP, grid2D<N> &nodeM) {
	for (int i = 1; i < nx+1; i++) {
		nodeP(i,1,1) = nodeP(i,0,1);
		nodeM(i,1,3) = nodeM(i,0,3);
		nodeP(i,ny,3) = nodeP(i,ny+1,3);
		nodeM(i,ny,1) = nodeM(i,ny+1,1);
	}
}

void setInletCorrectionReplica(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<N> &cellP, grid2D<N> &cellM) {
	for (int i = 1; i < nx+1; i++) {
		nodeP(i,1,1) = cellP(i-1,0,1);
		nodeM(i,1,3) = cellM(i,0,3);
		cellP(i,ny,3) = nodeP(i,ny+1,3);
		cellM(i,ny,1) = nodeM(i+1,ny+1,1);
	}
}


void applyMovingBoundary(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<1> &nodeZero) {
	myReal u, rho, feq[9];
	u = 0.04;
	for (int j = 1; j < ny+1; j++) {
		rho = 0.0;
		for (int k = 0; k<N; k++) {
			rho = rho + nodeP(nx,j,k) + nodeM(nx,j,k);
		}
		rho += nodeZero(nx,j,0);
		getFeq(feq, rho, u, 0);
		for (int k = 1; k<N; k++) {
			nodeM(nx,j,k) = feq[k+5];
		}
	}

//	for (int j = 1; j < ny+1; j++) {
//		rho = 0.0;
//		for (int k = 0; k<N; k++) {
//			rho = rho + nodeP(1,j,k) + nodeM(1,j,k);
//		}
//		rho += nodeZero(1,j,0);
//		getFeq(feq, rho, -u, 0);
//		for (int k = 1; k<N; k++) {
//			nodeP(1,j,k) = feq[k+1];
//		}
//	}
}

void applyMovingBoundaryReplica(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<1> &nodeZero, grid2D<N> &cellP, grid2D<N> &cellM, grid2D<1> &cellZero) {
	myReal u, rho, rho1, feq[9];
	u = 0.04;
	for (int j = 1; j < ny+1; j++) {
		rho = 0.0;
		rho1 = 0.0;
		for (int k = 0; k<N; k++) {
			rho = rho + cellP(nx,j,k) + cellM(nx,j,k);
			rho1 = rho1 + nodeP(nx,j,k) + nodeM(nx,j,k);
		}
		rho += cellZero(nx,j,0);
		rho1 += nodeZero(nx,j,0);
		getFeq(feq, rho, u, 0);
		for (int k = 1; k<N; k++) {
			cellM(nx,j,k) = feq[k+5];
		}
		getFeq(feq, rho1, u, 0);
		nodeM(nx,j,2) = feq[7];
	}

	for (int j = 1; j < ny+1; j++) {
		rho = 0.0;
		rho1 = 0.0;
		for (int k = 0; k<N; k++) {
			rho = rho + nodeP(1,j,k) + nodeM(1,j,k);
			rho1 = rho1 + cellP(1,j,k) + cellM(1,j,k);
		}
		rho += nodeZero(1,j,0);
		rho1 += cellZero(1,j,0);
		getFeq(feq, rho, -u, 0);
		for (int k = 1; k<N; k++) {
			nodeP(1,j,k) = feq[k+1];
		}
		getFeq(feq, rho1, -u, 0);
		cellP(1,j,2) = feq[3];
	}
}


void setObjectBC(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<1> &nodeZero, GridType &nodeType) {
	myReal temp[9], temp1[9], shift;
	for (int i = 1; i < nx+1; i++) {
		for (int j = 1; j < ny+1; j++) {
			for (int k = 1; k < 5; k++) {
				if (nodeType(i,j) == FLUID && nodeType(i+c_y[k],j+c_x[k]) == SOLID) {
						//cout<<"Cylinder reached at"<<i<<","<<j<<endl;
						copyFromGrid(i, j, nodeP, nodeM, nodeZero, temp);
						copyFromGrid(i+c_y[k],j+c_x[k], nodeP, nodeM, nodeZero, temp1);
						shift = temp[k+4];
						temp[k+4] = temp1[k];
						temp1[k] = shift;
						copyToGrid(i, j, nodeP, nodeM, nodeZero, temp);
						copyToGrid(i+c_y[k],j+c_x[k], nodeP, nodeM, nodeZero, temp1);
				}
				if (nodeType(i,j) == SOLID && nodeType(i+c_y[k],j+c_x[k]) == FLUID) {
						//cout<<"Cylinder exited"<<endl;
						copyFromGrid(i, j, nodeP, nodeM, nodeZero, temp);
						copyFromGrid(i+c_y[k],j+c_x[k], nodeP, nodeM, nodeZero, temp1);
						shift = temp[k+4];
						temp[k+4] = temp1[k];
						temp1[k] = shift;
						copyToGrid(i, j, nodeP, nodeM, nodeZero, temp);
						copyToGrid(i+c_y[k],j+c_x[k], nodeP, nodeM, nodeZero, temp1);
				}
			}
		}
	}
}


void setObjectBCReplica(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<1> &nodeZero, grid2D<N> &cellP, grid2D<N> &cellM, grid2D<1> &cellZero, GridType &nodeType, GridType &cellType) {
	myReal temp[9], temp1[9], shift;
	int k;
	//cout<<"Set cylinder bdry"<<endl;
	for (int i = 1; i < nx+1; i++) {
		for (int j = 1; j < ny+1; j++) {
			if (nodeType(i,j) == FLUID && cellType(i,j) == SOLID) {
				shift = nodeM(i,j,1);
				nodeM(i,j,1) = cellP(i,j,1);
				cellP(i,j,1) = shift;
			}
			if (nodeType(i,j) == FLUID && cellType(i,j-1) == SOLID) {
				shift = nodeM(i,j,3);
				nodeM(i,j,3) = cellP(i,j-1,3);
				cellP(i,j-1,3) = shift;
			}
			if (nodeType(i,j) == FLUID && nodeType(i,j+1) == SOLID) {
				shift = nodeM(i,j,0);
				nodeM(i,j,0) = nodeP(i,j+1,0);
				nodeP(i,j+1,0) = shift;
			}
			if (cellType(i,j) == FLUID && cellType(i,j+1) == SOLID) {
				shift = cellM(i,j,0);
				cellM(i,j,0) = cellP(i,j+1,0);
				cellP(i,j+1,0) = shift;
			}
			if (cellType(i,j) == FLUID && nodeType(i+1,j+1) == SOLID) {
				shift = cellM(i,j,1);
				cellM(i,j,1) = nodeP(i+1,j+1,1);
				nodeP(i+1,j+1,1) = shift;
			}
			if (cellType(i,j) == FLUID && nodeType(i+1,j) == SOLID) {
				shift = cellM(i,j,3);
				cellM(i,j,3) = nodeP(i+1,j,3);
				nodeP(i+1,j,3) = shift;
			}
			//Other way
			if (nodeType(i,j) == SOLID && cellType(i,j) == FLUID) {
				shift = nodeM(i,j,1);
				nodeM(i,j,1) = cellP(i,j,1);
				cellP(i,j,1) = shift;
			}
			if (nodeType(i,j) == SOLID && cellType(i,j-1) == FLUID) {
				shift = nodeM(i,j,3);
				nodeM(i,j,3) = cellP(i,j-1,3);
				cellP(i,j-1,3) = shift;
			}
			if (nodeType(i,j) == SOLID && nodeType(i,j+1) == FLUID) {
				shift = nodeM(i,j,0);
				nodeM(i,j,0) = nodeP(i,j+1,0);
				nodeP(i,j+1,0) = shift;
			}
			if (cellType(i,j) == SOLID && cellType(i,j+1) == FLUID) {
				shift = cellM(i,j,0);
				cellM(i,j,0) = cellP(i,j+1,0);
				cellP(i,j+1,0) = shift;
			}
			if (cellType(i,j) == SOLID && nodeType(i+1,j+1) == FLUID) {
				shift = cellM(i,j,1);
				cellM(i,j,1) = nodeP(i+1,j+1,1);
				nodeP(i+1,j+1,1) = shift;
			}
			if (cellType(i,j) == SOLID && nodeType(i+1,j) == FLUID) {
				shift = cellM(i,j,3);
				cellM(i,j,3) = nodeP(i+1,j,3);
				nodeP(i+1,j,3) = shift;
			}
			
			
			
		}
	}
	
}

void applyBounceBackReplica(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<N> &cellP, grid2D<N> &cellM) {
	//top and bottom
	for (int j = 1; j < ny+1; j++) {
		nodeM(nx,j,2) = nodeP(nx+1,j,2);
		cellM(nx,j,1) = cellP(nx+1,j,1);
		cellM(nx,j,2) = cellP(nx+1,j,2);
		cellM(nx,j,3) = cellP(nx,j,3);

		cellP(1,j,2) = cellM(0,j,2);
		nodeP(1,j,1) = nodeM(0,j,1);
		nodeP(1,j,2) = nodeM(0,j,2);
		nodeP(1,j,3) = nodeM(0,j,1);
	}
	//sides
//	for (int i = 1; i < nx+1; i++) {
//		nodeM(i,ny,0) = nodeP(i,ny+1,0);
//		cellM(i,ny,1) = cellP(i,ny+1,1);
//		cellM(i,ny,2) = cellP(i,ny+1,2);
//		cellP(i,ny,3) = cellP(i,ny+1,3);
//		
//		cellP(i,1,0) = cellM(i,0,0);
//		nodeP(i,1,1) = nodeM(i,0,1);
//		nodeP(i,1,2) = nodeM(i,0,2);
//		nodeM(i,1,3) = nodeP(i,0,3);
//	}	
}




//void setCorners(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<1> &nodeZero){
////	myReal f_grad[9], temp[9], feq[9], rho, ux, uy;
////	copyFromGrid(1, 1, nodeP, nodeM, nodeZero, temp);
////	getMoments(temp,rho,ux,uy);
////	getFeq(feq,rho,u_inlet,0);
////	nodeM(0,0,3) = feq[8];
////	copyFromGrid(nx, 1, nodeP, nodeM, nodeZero, temp);
////	getMoments(temp,rho,ux,uy);
////	getFeq(feq,rho,u_inlet,0);
////	nodeP(nx,1,1) = feq[2];
////	copyFromGrid(nx, ny, nodeP, nodeM, nodeZero, temp);
////	getMoments(temp,rho,ux,uy);
////	getFeq(feq,rho,u_inlet,0);
////	nodeP(0,0,3) = feq[4];
////	copyFromGrid(1, ny, nodeP, nodeM, nodeZero, temp);
////	getMoments(temp,rho,ux,uy);
////	getFeq(feq,rho,u_inlet,0);
////	nodeM(0,0,1) = feq[6];
//	nodeM(1,1,3) = nodeM(0,0,3);
//	nodeP(nx,1,1) = nodeP(nx+1,0,1);
//	nodeP(nx,ny,3) = nodeP(nx+1,ny+1,3);
//	nodeM(1,ny,1) = nodeM(0,ny+1,1);
//	
//}

//k = 2;
//			while (k < 5) {
//				if (nodeType(i,j) == FLUID && cellType(i+c_y1[k],j+c_x1[k]) == SOLID) {
//						//cout<<"Works";
//						//cout<<"Cylinder reached at"<<i<<","<<j<<endl;
//						copyFromGrid(i, j, nodeP, nodeM, nodeZero, temp);
//						copyFromGrid(i+c_y1[k],j+c_x1[k], cellP, cellM, cellZero, temp1);
//						shift = temp[k+4];
//						temp[k+4] = temp1[k];
//						temp1[k] = shift;
//						copyToGrid(i, j, nodeP, nodeM, nodeZero, temp);
//						copyToGrid(i+c_y1[k],j+c_x1[k], cellP, cellM, cellZero, temp1);
//				}
//				if (cellType(i,j) == FLUID && cellType(i+c_y1[k]+1,j+c_x1[k]+1) == SOLID) {
//						//cout<<"Cylinder reached at"<<i<<","<<j<<endl;
//						copyFromGrid(i, j, cellP, cellM, cellZero, temp);
//						copyFromGrid(i+c_y1[k],j+c_x1[k], nodeP, nodeM, nodeZero, temp1);
//						shift = temp[k+4];
//						temp[k+4] = temp1[k];
//						temp1[k] = shift;
//						copyToGrid(i, j, cellP, cellM, cellZero, temp);
//						copyToGrid(i+c_y1[k],j+c_x1[k], nodeP, nodeM, nodeZero, temp1);
//				}
//				if (nodeType(i,j) == SOLID && cellType(i+c_y1[k],j+c_x1[k]) == FLUID) {
//						//cout<<"Cylinder reached at"<<i<<","<<j<<endl;
//						copyFromGrid(i, j, nodeP, nodeM, nodeZero, temp);
//						copyFromGrid(i+c_y1[k],j+c_x1[k], cellP, cellM, cellZero, temp1);
//						shift = temp[k+4];
//						temp[k+4] = temp1[k];
//						temp1[k] = shift;
//						copyToGrid(i, j, nodeP, nodeM, nodeZero, temp);
//						copyToGrid(i+c_y1[k],j+c_x1[k], cellP, cellM, cellZero, temp1);
//				}
//				if (cellType(i,j) == SOLID && nodeType(i+c_y1[k]+1,j+c_x1[k]+1) == FLUID) {
//						//cout<<"Cylinder reached at"<<i<<","<<j<<endl;
//						copyFromGrid(i, j, cellP, cellM, cellZero, temp);
//						copyFromGrid(i+c_y1[k],j+c_x1[k], nodeP, nodeM, nodeZero, temp1);
//						shift = temp[k+4];
//						temp[k+4] = temp1[k];
//						temp1[k] = shift;
//						copyToGrid(i, j, cellP, cellM, cellZero, temp);
//						copyToGrid(i+c_y1[k],j+c_x1[k], nodeP, nodeM, nodeZero, temp1);
//				}
//				k = k+2;
//				
//			}
//			k = 1;
//			while (k < 4) {
//				if (nodeType(i,j) == FLUID && nodeType(i+c_y1[k],j+c_x1[k]) == SOLID) {
//						//cout<<"Cylinder reached at"<<i<<","<<j<<endl;
//						copyFromGrid(i, j, nodeP, nodeM, nodeZero, temp);
//						copyFromGrid(i+c_y1[k],j+c_x1[k], nodeP, nodeM, nodeZero, temp1);
//						shift = temp[k+4];
//						temp[k+4] = temp1[k];
//						temp1[k] = shift;
//						copyToGrid(i, j, nodeP, nodeM, nodeZero, temp);
//						copyToGrid(i+c_y1[k],j+c_x1[k],  nodeP, nodeM, nodeZero, temp1);
//				}
//				if (cellType(i,j) == FLUID && cellType(i+c_y1[k],j+c_x1[k]) == SOLID) {
//						//cout<<"Cylinder reached at"<<i<<","<<j<<endl;
//						copyFromGrid(i, j, cellP, cellM, cellZero, temp);
//						copyFromGrid(i+c_y1[k],j+c_x1[k], cellP, cellM, cellZero, temp1);
//						shift = temp[k+4];
//						temp[k+4] = temp1[k];
//						temp1[k] = shift;
//						copyToGrid(i, j, cellP, cellM, cellZero, temp);
//						copyToGrid(i+c_y1[k],j+c_x1[k], cellP, cellM, cellZero, temp1);
//				}
//				if (nodeType(i,j) == SOLID && nodeType(i+c_y1[k],j+c_x1[k]) == FLUID) {
//						//cout<<"Cylinder reached at"<<i<<","<<j<<endl;
//						copyFromGrid(i, j, nodeP, nodeM, nodeZero, temp);
//						copyFromGrid(i+c_y1[k],j+c_x1[k], nodeP, nodeM, nodeZero, temp1);
//						shift = temp[k+4];
//						temp[k+4] = temp1[k];
//						temp1[k] = shift;
//						copyToGrid(i, j, nodeP, nodeM, nodeZero, temp);
//						copyToGrid(i+c_y1[k],j+c_x1[k],  nodeP, nodeM, nodeZero, temp1);
//				}
//				if (cellType(i,j) == SOLID && cellType(i+c_y1[k],j+c_x1[k]) == FLUID) {
//						//cout<<"Cylinder reached at"<<i<<","<<j<<endl;
//						copyFromGrid(i, j, cellP, cellM, cellZero, temp);
//						copyFromGrid(i+c_y1[k],j+c_x1[k], cellP, cellM, cellZero, temp1);
//						shift = temp[k+4];
//						temp[k+4] = temp1[k];
//						temp1[k] = shift;
//						copyToGrid(i, j, cellP, cellM, cellZero, temp);
//						copyToGrid(i+c_y1[k],j+c_x1[k], cellP, cellM, cellZero, temp1);
//				}
//				k = k+2;
//			}


//			if (nodeType(i,j) == FLUID && cellType(i,j) == SOLID) {
//				shift = nodeM(i,j,1);
//				nodeM(i,j,1) = cellP(i,j,1);
//				cellP(i,j,1) = shift;
//			}
//			if (nodeType(i,j) == FLUID && cellType(i,j-1) == SOLID) {
//				shift = nodeM(i,j,3);
//				nodeM(i,j,3) = cellP(i,j,3);
//				cellP(i,j,3) = shift;
//			}
//			if (nodeType(i,j) == FLUID && nodeType(i,j+1) == SOLID) {
//				shift = nodeM(i,j,0);
//				nodeM(i,j,0) = nodeP(i,j,0);
//				nodeP(i,j,0) = shift;
//			}
//			if (cellType(i,j) == FLUID && cellType(i,j+1) == SOLID) {
//				shift = cellM(i,j,0);
//				cellM(i,j,0) = cellP(i,j,0);
//				cellP(i,j,0) = shift;
//			}
//			if (cellType(i,j) == FLUID && nodeType(i+1,j+1) == SOLID) {
//				shift = cellM(i,j,1);
//				cellM(i,j,1) = nodeP(i,j,1);
//				nodeP(i,j,1) = shift;
//			}
//			if (cellType(i,j) == FLUID && nodeType(i+1,j) == SOLID) {
//				shift = cellM(i,j,3);
//				cellM(i,j,3) = nodeP(i,j,3);
//				nodeP(i,j,3) = shift;
//			}
#endif
