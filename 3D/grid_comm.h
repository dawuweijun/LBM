#ifndef _COMMUNICATE_H_
#define  _COMMUNICATE_H_
#include<iostream>
#include<stdlib.h>

#include "gridV.h"
#include "gridV_Probe.h"
#include "globals3D.h"
#include "gridV_Pack.h"

#include"mpi.h"
#include<procInfo3D.h> 

void recvCommunicate3DFace1(int myRank,  int *myneighbours, MPI::Request *request, myReal *tmpDataReceive1M, myReal *tmpDataReceive1P, int m1, int m2, int m3){ 
      int tag = 0;
	  int dataSizePlane23 = 18*m2*m3;
      request[0]=MPI::COMM_WORLD.Irecv(tmpDataReceive1P, dataSizePlane23, MPI::DOUBLE, myneighbours[NEB_M1_ZERO_ZERO], tag); 
      tag = 1;
      request[1]=MPI::COMM_WORLD.Irecv(tmpDataReceive1M, dataSizePlane23, MPI::DOUBLE, myneighbours[NEB_P1_ZERO_ZERO], tag); 
     } 
     
void recvCommunicate3DFace2(int myRank,  int *myneighbours, MPI::Request *request, myReal *tmpDataReceive2M, myReal *tmpDataReceive2P, int m1, int m2, int m3){
      int tag = 2;
	  int dataSizePlane13 = 18*m1*m3;
      request[0]=MPI::COMM_WORLD.Irecv(tmpDataReceive2P, dataSizePlane13, MPI::DOUBLE, myneighbours[NEB_ZERO_M1_ZERO], tag);  
      tag = 3;
      request[1]=MPI::COMM_WORLD.Irecv(tmpDataReceive2M, dataSizePlane13, MPI::DOUBLE, myneighbours[NEB_ZERO_P1_ZERO], tag); 
    }  
      
void recvCommunicate3DFace3(int myRank,  int *myneighbours, MPI::Request *request, myReal *tmpDataReceive3M, myReal *tmpDataReceive3P, int m1, int m2, int m3){
      int tag = 4;
	  int dataSizePlane12 = 18*m1*m2;
      request[0]=MPI::COMM_WORLD.Irecv(tmpDataReceive3P, dataSizePlane12, MPI::DOUBLE, myneighbours[NEB_ZERO_ZERO_M1], tag);
      tag = 5;
      request[1]=MPI::COMM_WORLD.Irecv(tmpDataReceive3M, dataSizePlane12, MPI::DOUBLE, myneighbours[NEB_ZERO_ZERO_P1], tag);
    }


void sendCommunicateFace1(int myRank,  int *myneighbours, MPI::Request *request, myReal *tmpDataSend1P, myReal *tmpDataSend1M, int m1, int m2, int m3){
    int tag = 0;
	int dataSizePlane23 = 18*m2*m3;
    request[0]= MPI::COMM_WORLD.Isend(tmpDataSend1P,dataSizePlane23 , MPI::DOUBLE, myneighbours[NEB_P1_ZERO_ZERO], tag); 
    tag  = 1;
    request[1]= MPI::COMM_WORLD.Isend(tmpDataSend1M,dataSizePlane23 , MPI::DOUBLE, myneighbours[NEB_M1_ZERO_ZERO], tag); 
    }
    
void sendCommunicateFace2(int myRank,  int *myneighbours, MPI::Request *request, myReal *tmpDataSend2P, myReal *tmpDataSend2M, int m1, int m2, int m3){
     int tag = 2;
	 int dataSizePlane13 = 18*m1*m3;
     request[0]= MPI::COMM_WORLD.Isend(tmpDataSend2P,dataSizePlane13 , MPI::DOUBLE, myneighbours[NEB_ZERO_P1_ZERO], tag);
     tag  = 3;
     request[1]=MPI::COMM_WORLD.Isend(tmpDataSend2M,dataSizePlane13 , MPI::DOUBLE, myneighbours[NEB_ZERO_M1_ZERO], tag); 
    }
    
void sendCommunicateFace3(int myRank,  int *myneighbours, MPI::Request *request, myReal *tmpDataSend3P, myReal *tmpDataSend3M, int m1, int m2, int m3){
      int tag = 4;
	  int dataSizePlane12 = 18*m1*m2;
      request[0]=MPI::COMM_WORLD.Isend(tmpDataSend3P,dataSizePlane12 , MPI::DOUBLE, myneighbours[NEB_ZERO_ZERO_P1], tag);
      tag  = 5;
      request[1]=MPI::COMM_WORLD.Isend(tmpDataSend3M,dataSizePlane12 , MPI::DOUBLE, myneighbours[NEB_ZERO_ZERO_M1], tag);
    }
#endif	