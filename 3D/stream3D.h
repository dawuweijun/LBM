//STREAM FUNCTION: ADVECTION
#ifndef _STREAM_3D_H
#define _STREAM_3D_H
#include <iostream>
#include "gridV.h"
#include "globals3D.h"

//using namespace std; 

void stream(lbgrid &node, lbgrid &cell) {
	//SCP SHIFT
	//Positive (1,0,0)...
	for (int i3 = node.SCP.iE3; i3 >= node.SCP.iB3; i3--) {
		for (int i2 = node.SCP.iE2; i2 >= node.SCP.iB2; i2--) {
			for (int i1 = node.SCP.iE1; i1 >= node.SCP.iB1; i1--) {
				 //if(node.SCP(0,i2,i3,DV_ZERO_P1_ZERO)==0)cout<<"Stream"<<i3<<endl;
				 node.SCP(i1,i2,i3, DV_P1_ZERO_ZERO) = node.SCP(i1-node.SCP.nB1,i2,i3,DV_P1_ZERO_ZERO);       
				 node.SCP(i1,i2,i3, DV_ZERO_P1_ZERO) = node.SCP(i1,i2-node.SCP.nB2,i3,DV_ZERO_P1_ZERO);
		 		 node.SCP(i1,i2,i3, DV_ZERO_ZERO_P1) = node.SCP(i1,i2,i3-node.SCP.nB3,DV_ZERO_ZERO_P1);
			}
		}
	}
	for (int i3 = node.SCP.iE3; i3 >= node.SCP.iB3; i3--) {
		for (int i2 = node.SCP.iE2; i2 >= node.SCP.iB2; i2--) {
			for (int i1 = node.SCP.iE1; i1 >= node.SCP.iB1; i1--) {
				 cell.SCP(i1,i2,i3, DV_P1_ZERO_ZERO) = cell.SCP(i1-cell.SCP.nB1,i2,i3,DV_P1_ZERO_ZERO);       
				 cell.SCP(i1,i2,i3, DV_ZERO_P1_ZERO) = cell.SCP(i1,i2-cell.SCP.nB2,i3,DV_ZERO_P1_ZERO);
		 		 cell.SCP(i1,i2,i3, DV_ZERO_ZERO_P1) = cell.SCP(i1,i2,i3-cell.SCP.nB3,DV_ZERO_ZERO_P1);
			}
		}
	}
	//Negative (-1,0,0)...

	for (int i3 = node.SCM.iB3; i3 <= node.SCM.iE3; i3++) {
		for (int i2 = node.SCM.iB2; i2 <= node.SCM.iE2; i2++) {
			for (int i1 = node.SCM.iB1; i1 <= node.SCM.iE1; i1++) {
				 node.SCM(i1,i2,i3,DV_M1_ZERO_ZERO) = node.SCM(i1+node.SCM.nB1,i2,i3,DV_M1_ZERO_ZERO);
	      			 node.SCM(i1,i2,i3,DV_ZERO_M1_ZERO) = node.SCM(i1,i2+node.SCM.nB2,i3,DV_ZERO_M1_ZERO);
	       			 node.SCM(i1,i2,i3,DV_ZERO_ZERO_M1) = node.SCM(i1,i2,i3+node.SCM.nB3,DV_ZERO_ZERO_M1);
				 
			}
		}
	}
	
	for (int i3 = node.SCM.iB3; i3 <= node.SCM.iE3; i3++) {
		for (int i2 = node.SCM.iB2; i2 <= node.SCM.iE2; i2++) {
			for (int i1 = node.SCM.iB1; i1 <= node.SCM.iE1; i1++) {
//				 //count++;
				 cell.SCM(i1,i2,i3,DV_M1_ZERO_ZERO) = cell.SCM(i1+cell.SCM.nB1,i2,i3,DV_M1_ZERO_ZERO);
	      			 cell.SCM(i1,i2,i3,DV_ZERO_M1_ZERO) = cell.SCM(i1,i2+cell.SCM.nB2,i3,DV_ZERO_M1_ZERO);
	       			 cell.SCM(i1,i2,i3,DV_ZERO_ZERO_M1) = cell.SCM(i1,i2,i3+cell.SCM.nB3,DV_ZERO_ZERO_M1);
			}
		}
	}
				
	//FCC SHIFT
	//12
	for (int i3 = node.FCC12.iB3; i3 <= node.FCC12.iE3; i3++) {
		for (int i2 = node.FCC12.iB2; i2 <= node.FCC12.iE2; i2++) {
			for (int i1 = node.FCC12.iB1; i1 <= node.FCC12.iE1; i1++) {
				node.FCC12(i1,i2,i3,DV_P1_M1_ZERO) = node.FCC12(i1-node.FCC12.nB1,i2+node.FCC12.nB2,i3,DV_P1_M1_ZERO);
				node.FCC12(i1,i2,i3,DV_M1_M1_ZERO) = node.FCC12(i1+node.FCC12.nB1,i2+node.FCC12.nB2,i3,DV_M1_M1_ZERO);
 		}
		}
		for (int i2 = node.FCC12.iE2; i2 >= node.FCC12.iB2; i2--) {
			for (int i1 = node.FCC12.iE1; i1 >= node.FCC12.iB1; i1--) {
				node.FCC12(i1,i2,i3,DV_P1_P1_ZERO) = node.FCC12(i1-node.FCC12.nB1,i2-node.FCC12.nB2,i3,DV_P1_P1_ZERO);
				node.FCC12(i1,i2,i3,DV_M1_P1_ZERO) = node.FCC12(i1+node.FCC12.nB1,i2-node.FCC12.nB2,i3,DV_M1_P1_ZERO);
 		}
		}
	}
	
	
	
	for (int i3 = node.FCC12.iB3; i3 <= node.FCC12.iE3; i3++) {
		for (int i2 = node.FCC12.iB2; i2 <= node.FCC12.iE2; i2++) {
			for (int i1 = node.FCC12.iB1; i1 <= node.FCC12.iE1; i1++) {
		 
				cell.FCC12(i1,i2,i3,DV_P1_M1_ZERO) = cell.FCC12(i1-cell.FCC12.nB1,i2+cell.FCC12.nB2,i3,DV_P1_M1_ZERO);
				cell.FCC12(i1,i2,i3,DV_M1_M1_ZERO) = cell.FCC12(i1+cell.FCC12.nB1,i2+cell.FCC12.nB2,i3,DV_M1_M1_ZERO);
			}
		}
		for (int i2 = node.FCC12.iE2; i2 >= node.FCC12.iB2; i2--) {
			for (int i1 = node.FCC12.iE1; i1 >= node.FCC12.iB1; i1--) {
			 	cell.FCC12(i1,i2,i3,DV_P1_P1_ZERO) = cell.FCC12(i1-cell.FCC12.nB1,i2-cell.FCC12.nB2,i3,DV_P1_P1_ZERO);
				cell.FCC12(i1,i2,i3,DV_M1_P1_ZERO) = cell.FCC12(i1+cell.FCC12.nB1,i2-cell.FCC12.nB2,i3,DV_M1_P1_ZERO);
			}
		}
	}
	
	
	
	
	
//////////	//13
	for (int i2 = node.FCC13.iB2; i2 <= node.FCC13.iE2; i2++) {
		for (int i3 = node.FCC13.iB3; i3 <= node.FCC13.iE3; i3++) {
			for (int i1 = node.FCC13.iB1; i1 <= node.FCC13.iE1; i1++) {
				node.FCC13(i1,i2,i3,DV_P1_ZERO_M1) = node.FCC13(i1-node.FCC13.nB1,i2,i3+node.FCC13.nB3,DV_P1_ZERO_M1);
				node.FCC13(i1,i2,i3,DV_M1_ZERO_M1) = node.FCC13(i1+node.FCC13.nB1,i2,i3+node.FCC13.nB3,DV_M1_ZERO_M1);

				cell.FCC13(i1,i2,i3,DV_P1_ZERO_M1) = cell.FCC13(i1-cell.FCC13.nB1,i2,i3+cell.FCC13.nB3,DV_P1_ZERO_M1);
				cell.FCC13(i1,i2,i3,DV_M1_ZERO_M1) = cell.FCC13(i1+cell.FCC13.nB1,i2,i3+cell.FCC13.nB3,DV_M1_ZERO_M1);
			}
		}
		for (int i3 = node.FCC13.iE3; i3 >= node.FCC13.iB3; i3--) {
			for (int i1 = node.FCC13.iE1; i1 >= node.FCC13.iB1; i1--) {
				node.FCC13(i1,i2,i3,DV_P1_ZERO_P1) = node.FCC13(i1-node.FCC13.nB1,i2,i3-node.FCC13.nB3,DV_P1_ZERO_P1);
				node.FCC13(i1,i2,i3,DV_M1_ZERO_P1) = node.FCC13(i1+node.FCC13.nB1,i2,i3-node.FCC13.nB3,DV_M1_ZERO_P1);

				cell.FCC13(i1,i2,i3,DV_P1_ZERO_P1) = cell.FCC13(i1-cell.FCC13.nB1,i2,i3-cell.FCC13.nB3,DV_P1_ZERO_P1);
				cell.FCC13(i1,i2,i3,DV_M1_ZERO_P1) = cell.FCC13(i1+cell.FCC13.nB1,i2,i3-cell.FCC13.nB3,DV_M1_ZERO_P1);
			}
		}
	}
////////	//23
	for (int i1 = node.FCC23.iB1; i1 <= node.FCC23.iE1; i1++) {
		for (int i3 = node.FCC23.iB3; i3 <= node.FCC23.iE3; i3++) {
			for (int i2 = node.FCC23.iB2; i2 <= node.FCC23.iE2; i2++) {
				node.FCC23(i1,i2,i3,DV_ZERO_P1_M1) = node.FCC23(i1,i2-node.FCC23.nB2,i3+node.FCC23.nB3,DV_ZERO_P1_M1);
				node.FCC23(i1,i2,i3,DV_ZERO_M1_M1) = node.FCC23(i1,i2+node.FCC23.nB2,i3+node.FCC23.nB3,DV_ZERO_M1_M1);

				cell.FCC23(i1,i2,i3,DV_ZERO_P1_M1) = cell.FCC23(i1,i2-cell.FCC23.nB2,i3+cell.FCC23.nB3,DV_ZERO_P1_M1);
				cell.FCC23(i1,i2,i3,DV_ZERO_M1_M1) = cell.FCC23(i1,i2+cell.FCC23.nB2,i3+cell.FCC23.nB3,DV_ZERO_M1_M1);
			}
		}
		for (int i3 = node.FCC23.iE3; i3 >= node.FCC23.iB3; i3--) {
			for (int i2= node.FCC23.iE2; i2 >= node.FCC23.iB2; i2--) {
				node.FCC23(i1,i2,i3,DV_ZERO_P1_P1) = node.FCC23(i1,i2-node.FCC23.nB2,i3-node.FCC23.nB3,DV_ZERO_P1_P1);
				node.FCC23(i1,i2,i3,DV_ZERO_M1_P1) = node.FCC23(i1,i2+node.FCC23.nB2,i3-node.FCC23.nB3,DV_ZERO_M1_P1);

				cell.FCC23(i1,i2,i3,DV_ZERO_P1_P1) = cell.FCC23(i1,i2-cell.FCC23.nB2,i3-cell.FCC23.nB3,DV_ZERO_P1_P1);
				cell.FCC23(i1,i2,i3,DV_ZERO_M1_P1) = cell.FCC23(i1,i2+cell.FCC23.nB2,i3-cell.FCC23.nB3,DV_ZERO_M1_P1);
			}
		}
	}

////	//BCC SHIFT---->node. from cell. and cell. from node.
////	//Positive 	
	for (int i3 = node.BCCP.iE3; i3 >= node.BCCP.iB3; i3--) {
		for (int i2 = node.BCCP.iE2; i2 >= node.BCCP.iB2; i2--) {
			for (int i1 = node.BCCP.iE1; i1 >= node.BCCP.iB1; i1--) {				
				cell.BCCP(i1,i2,i3,DV_P1_P1_P1) = node.BCCP(i1,i2,i3,DV_P1_P1_P1);
				cell.BCCP(i1,i2,i3,DV_M1_P1_P1) = node.BCCP(i1+node.BCCP.nB1,i2,i3,DV_M1_P1_P1);
				cell.BCCP(i1,i2,i3,DV_M1_M1_P1) = node.BCCP(i1+node.BCCP.nB1,i2+node.BCCP.nB2,i3,DV_M1_M1_P1);
				cell.BCCP(i1,i2,i3,DV_P1_M1_P1) = node.BCCP(i1,i2+node.BCCP.nB2,i3,DV_P1_M1_P1);
			}
		}
		for (int i2 = node.BCCP.iE2; i2 >= node.BCCP.iB2; i2--) {
			for (int i1 = node.BCCP.iE1; i1 >= node.BCCP.iB1; i1--) {
				node.BCCP(i1,i2,i3,DV_P1_P1_P1) = cell.BCCP(i1-cell.BCCP.nB1,i2-cell.BCCP.nB2,i3-cell.BCCP.nB3,DV_P1_P1_P1);
				node.BCCP(i1,i2,i3,DV_M1_P1_P1) = cell.BCCP(i1,i2-cell.BCCP.nB2,i3-cell.BCCP.nB3,DV_M1_P1_P1);
				node.BCCP(i1,i2,i3,DV_M1_M1_P1) = cell.BCCP(i1,i2,i3-cell.BCCP.nB3,DV_M1_M1_P1);
				node.BCCP(i1,i2,i3,DV_P1_M1_P1) = cell.BCCP(i1-cell.BCCP.nB1,i2,i3-cell.BCCP.nB3,DV_P1_M1_P1);
			}
		}
	}
	//Negative
	for (int i3 = node.BCCM.iB3; i3 <= node.BCCM.iE3; i3++) {
		for (int i2 = node.BCCM.iB2; i2 <= node.BCCM.iE2; i2++) {
			for (int i1 = node.BCCM.iB1; i1 <= node.BCCM.iE1; i1++) {
				node.BCCM(i1,i2,i3,DV_M1_M1_M1) = cell.BCCM(i1,i2,i3,DV_M1_M1_M1);
				node.BCCM(i1,i2,i3,DV_M1_P1_M1) = cell.BCCM(i1,i2-cell.BCCM.nB2,i3,DV_M1_P1_M1); 
				node.BCCM(i1,i2,i3,DV_P1_M1_M1) = cell.BCCM(i1-cell.BCCM.nB1,i2,i3,DV_P1_M1_M1);
				node.BCCM(i1,i2,i3,DV_P1_P1_M1) = cell.BCCM(i1-cell.BCCM.nB1,i2-cell.BCCM.nB2,i3,DV_P1_P1_M1);
			}
		}
		for (int i2 = node.BCCM.iB2; i2 <= node.BCCM.iE2; i2++) {
			for (int i1 = node.BCCM.iB1; i1 <= node.BCCM.iE1; i1++) {
				cell.BCCM(i1,i2,i3,DV_M1_M1_M1) = node.BCCM(i1+node.BCCM.nB1,i2+node.BCCM.nB2,i3+node.BCCM.nB3,DV_M1_M1_M1);
				cell.BCCM(i1,i2,i3,DV_M1_P1_M1) = node.BCCM(i1+node.BCCM.nB1,i2,i3+node.BCCM.nB3,DV_M1_P1_M1);
				cell.BCCM(i1,i2,i3,DV_P1_M1_M1) = node.BCCM(i1,i2+node.BCCM.nB2,i3+node.BCCM.nB3,DV_P1_M1_M1);
				cell.BCCM(i1,i2,i3,DV_P1_P1_M1) = node.BCCM(i1,i2,i3+node.BCCM.nB3,DV_P1_P1_M1);
			}
		}
	}

//Stream finished
}
		
#endif	
			
				
	
