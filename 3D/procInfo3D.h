#ifndef _PROCINFO3D_H_
#define _PROCINFO3D_H_
#include"mpi.h"
#include<iostream>
#define TOPOLOGY_DIM 3
 enum direction{Y1, Y2, Y3};

#define ERROR_MSG()  std::cout<<"I am out o here"<<std::endl;
 enum nebName {
      NEB_ZERO_ZERO_ZERO,
    // Energy shell =1
      NEB_P1_ZERO_ZERO ,
      NEB_M1_ZERO_ZERO ,
      NEB_ZERO_P1_ZERO ,
      NEB_ZERO_M1_ZERO ,
      NEB_ZERO_ZERO_P1 ,
      NEB_ZERO_ZERO_M1,
    //Energy Shell =2
      NEB_P1_P1_ZERO,
      NEB_M1_M1_ZERO,
      NEB_P1_M1_ZERO,
      NEB_M1_P1_ZERO,
      NEB_P1_ZERO_P1,
      NEB_M1_ZERO_M1,
      NEB_P1_ZERO_M1,
      NEB_M1_ZERO_P1,
      NEB_ZERO_P1_P1,
      NEB_ZERO_M1_M1,
      NEB_ZERO_P1_M1,
      NEB_ZERO_M1_P1,
   //Energy Shell =3
      NEB_P1_P1_P1,
      NEB_M1_M1_M1,
      NEB_P1_M1_P1,
      NEB_M1_P1_M1,
      NEB_P1_M1_M1,
      NEB_M1_P1_P1,
      NEB_M1_M1_P1,
      NEB_P1_P1_M1		
  };  
  
struct procInfo3D{
       procInfo3D(int n1, int n2, int n3, bool wrap1=true, bool wrap2=true, bool wrap3=true, bool reorderGiven =true){
	  // Number of Processor in each direction
	  dimSizes [Y1] = n1;
	  dimSizes [Y2] = n2;
	  dimSizes [Y3] = n3; 

	  wrapAround[Y1]= wrap1;
	  wrapAround[Y2]= wrap2;
	  wrapAround[Y3]= wrap3; 
	  reorder = reorderGiven;
	  gridComm  =  MPI::COMM_WORLD.Create_cart (TOPOLOGY_DIM, dimSizes, wrapAround, reorder);
	  int myGridRank =gridComm.Get_rank() ;
	  gridComm.Get_coords( myGridRank , TOPOLOGY_DIM ,coordinates ) ;
	  for(int i=0;i<27;i++) nList[i]=-300;
      }
      bool reorder ; // 1 means MPI can decide how to arrange neighbours
      bool wrapAround[TOPOLOGY_DIM] ;// Periodicity info
      int dimSizes[TOPOLOGY_DIM];// Size in each dir 
      int coordinates[TOPOLOGY_DIM] ; // Virtual coordinate of the processor
      MPI::Cartcomm gridComm;
      void getnList();
      int nList[27];
      int ncord[27];
      int coord[3];
};

#include<procInfo3D.C>

#endif
