#ifndef _STREAM_H
#define _STREAM_H

#include<iostream>
#include<math.h>
#include "grid2D.h"
#include "globals.h"



void stream(grid2D<N> &nodeP, grid2D<N> &nodeM) {
	for (int i = 1; i < nx+1; i++) {
		for (int j = 1; j < ny+1; j++) {
			for (int k = 0; k < N; k++) {
				nodeM(i,j,k) = nodeM(i-c_y[k+5],j-c_x[k+5],k);
			}
		}
	}
	for (int i = nx	; i > 0; i--) {
		for (int j = ny; j > 0; j--) {
			for (int k = 0; k < N; k++) {
				nodeP(i,j,k) = nodeP(i-c_y[k+1],j-c_x[k+1],k);
			}
		}
	}
}

void streamReplica(grid2D<N> &nodeP, grid2D<N> &nodeM, grid2D<N> &cellP, grid2D<N> &cellM) {
	for (int i = nx; i > 0; i--) {
		for (int j = ny; j > 0; j--) {
			cellP(i,j,1) = nodeP(i,j,1);
			cellP(i,j,3) = nodeP(i,j+1,3);
			cellP(i,j,0) = cellP(i,j-1,0);
			cellP(i,j,2) = cellP(i-1,j,2);	
		}
		for (int j = ny; j > 0; j--) {
			nodeP(i,j,0) = nodeP(i,j-1,0);
			nodeP(i,j,2) = nodeP(i-1,j,2);
			nodeP(i,j,1) = cellP(i-1,j-1,1);
			nodeP(i,j,3) = cellP(i-1,j,3);
		}
	}
	for (int i = 1; i < nx+1; i++) {
		for (int j = 1; j < ny+1; j++) {
			nodeM(i,j,0) = nodeM(i,j+1,0);
			nodeM(i,j,2) = nodeM(i+1,j,2);
			nodeM(i,j,1) = cellM(i,j,1);
			nodeM(i,j,3) = cellM(i,j-1,3);
		}
		for (int j = 1; j < ny+1; j++) {
			cellM(i,j,1) = nodeM(i+1,j+1,1);
			cellM(i,j,3) = nodeM(i+1,j,3);
			cellM(i,j,0) = cellM(i,j+1,0);
			cellM(i,j,2) = cellM(i+1,j,2);
		}
	}

}

#endif
