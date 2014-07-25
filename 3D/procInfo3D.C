void procInfo3D::getnList(){ 
    
	  int   gridRank ;

      gridRank =gridComm.Get_cart_rank(coordinates) ;
      nList[NEB_ZERO_ZERO_ZERO] = gridRank;

      int disp = 1;

      /************Shift create neighbour fist output =Rank -disp and second is rank + disp  in direction i  ***********/

      gridComm.Shift(Y1, disp, nList[NEB_M1_ZERO_ZERO],  nList[NEB_P1_ZERO_ZERO]);
      gridComm.Shift(Y2, disp, nList[NEB_ZERO_M1_ZERO],  nList[NEB_ZERO_P1_ZERO]);
      gridComm.Shift(Y3, disp, nList[NEB_ZERO_ZERO_M1],  nList[NEB_ZERO_ZERO_P1]);  

int coordDesired[TOPOLOGY_DIM], coordX1N[TOPOLOGY_DIM], coordX2N[TOPOLOGY_DIM], coordX3N[TOPOLOGY_DIM];

gridComm.Get_coords(nList[NEB_P1_ZERO_ZERO], TOPOLOGY_DIM, coordX1N);
gridComm.Get_coords(nList[NEB_ZERO_P1_ZERO], TOPOLOGY_DIM, coordX2N);
gridComm.Get_coords(nList[NEB_ZERO_ZERO_P1], TOPOLOGY_DIM, coordX3N);

coordDesired[Y1]= coordX1N[Y1];
coordDesired[Y2]= coordX2N[Y2];
coordDesired[Y3]= coordX3N[Y3];
/*****************BODY CENTERED*********************/
nList[NEB_P1_P1_P1]=gridComm.Get_cart_rank(coordDesired);
       /*********************/
gridComm.Get_coords(nList[NEB_M1_ZERO_ZERO], TOPOLOGY_DIM, coordX1N);
coordDesired[Y1]= coordX1N[Y1];
nList[NEB_M1_P1_P1]=gridComm.Get_cart_rank(coordDesired);
       /*********************/
gridComm.Get_coords(nList[NEB_ZERO_M1_ZERO], TOPOLOGY_DIM, coordX2N);
coordDesired[Y2]= coordX2N[Y2];
nList[NEB_M1_M1_P1]=gridComm.Get_cart_rank(coordDesired);
        /*********************/       
 gridComm.Get_coords(nList[NEB_ZERO_ZERO_M1], TOPOLOGY_DIM, coordX3N);
coordDesired[Y3]= coordX3N[Y3];
nList[NEB_M1_M1_M1]=gridComm.Get_cart_rank(coordDesired);
        /*********************/
gridComm.Get_coords(nList[NEB_P1_ZERO_ZERO], TOPOLOGY_DIM, coordX1N);
gridComm.Get_coords(nList[NEB_ZERO_P1_ZERO], TOPOLOGY_DIM, coordX2N);
coordDesired[Y1]= coordX1N[Y1];
coordDesired[Y2]= coordX2N[Y2];
nList[NEB_P1_P1_M1]=gridComm.Get_cart_rank(coordDesired);
        /*********************/
gridComm.Get_coords(nList[NEB_ZERO_ZERO_P1], TOPOLOGY_DIM, coordX3N);
gridComm.Get_coords(nList[NEB_ZERO_M1_ZERO], TOPOLOGY_DIM, coordX2N);
coordDesired[Y2]= coordX2N[Y2];
coordDesired[Y3]= coordX3N[Y3];
nList[NEB_P1_M1_P1]=gridComm.Get_cart_rank(coordDesired);
         /*********************/
gridComm.Get_coords(nList[NEB_ZERO_ZERO_M1], TOPOLOGY_DIM, coordX3N);
coordDesired[Y3]= coordX3N[Y3];
nList[NEB_P1_M1_M1]=gridComm.Get_cart_rank(coordDesired);
       /*********************/
gridComm.Get_coords(nList[NEB_M1_ZERO_ZERO], TOPOLOGY_DIM, coordX1N);
gridComm.Get_coords(nList[NEB_ZERO_P1_ZERO], TOPOLOGY_DIM, coordX2N);
coordDesired[Y1]= coordX1N[Y1];
coordDesired[Y2]= coordX2N[Y2];
nList[NEB_M1_P1_M1]=gridComm.Get_cart_rank(coordDesired);
           /*********************/
/*****************FACE CENTERED*********************/
gridComm.Get_coords(nList[NEB_ZERO_ZERO_ZERO], TOPOLOGY_DIM, coordX3N);
coordDesired[Y3]= coordX3N[Y3];
nList[NEB_M1_P1_ZERO]=gridComm.Get_cart_rank(coordDesired);
        /*********************/
gridComm.Get_coords(nList[NEB_P1_ZERO_ZERO], TOPOLOGY_DIM, coordX1N);
coordDesired[Y1]= coordX1N[Y1];
nList[NEB_P1_P1_ZERO]=gridComm.Get_cart_rank(coordDesired);
    /*********************/
gridComm.Get_coords(nList[NEB_ZERO_M1_ZERO], TOPOLOGY_DIM, coordX2N);
coordDesired[Y2]= coordX2N[Y2];
nList[NEB_P1_M1_ZERO]=gridComm.Get_cart_rank(coordDesired);
  /*********************/
gridComm.Get_coords(nList[NEB_M1_ZERO_ZERO], TOPOLOGY_DIM, coordX1N);
coordDesired[Y1]= coordX1N[Y1];
nList[NEB_M1_M1_ZERO]=gridComm.Get_cart_rank(coordDesired);
  /*********************/ 
gridComm.Get_coords(nList[NEB_ZERO_ZERO_P1], TOPOLOGY_DIM, coordX3N);
coordDesired[Y1]= coordX3N[Y1];
coordDesired[Y3]= coordX3N[Y3];
nList[NEB_ZERO_M1_P1]=gridComm.Get_cart_rank(coordDesired);
  /*********************/ 
gridComm.Get_coords(nList[NEB_ZERO_ZERO_M1], TOPOLOGY_DIM, coordX3N);
coordDesired[Y3]= coordX3N[Y3];
nList[NEB_ZERO_M1_M1]=gridComm.Get_cart_rank(coordDesired);
  /*********************/ 
gridComm.Get_coords(nList[NEB_ZERO_P1_ZERO], TOPOLOGY_DIM, coordX2N);
coordDesired[Y2]= coordX2N[Y2];
nList[NEB_ZERO_P1_M1]=gridComm.Get_cart_rank(coordDesired);
  /*********************/ 

gridComm.Get_coords(nList[NEB_ZERO_ZERO_P1], TOPOLOGY_DIM, coordX3N);
coordDesired[Y3]= coordX3N[Y3];
nList[NEB_ZERO_P1_P1]=gridComm.Get_cart_rank(coordDesired);
  /*********************/  
gridComm.Get_coords(nList[NEB_P1_ZERO_ZERO], TOPOLOGY_DIM, coordX1N);
coordDesired[Y1]= coordX1N[Y1];
coordDesired[Y2]= coordX1N[Y2];
nList[NEB_P1_ZERO_P1]=gridComm.Get_cart_rank(coordDesired);
  /*********************/ 
gridComm.Get_coords(nList[NEB_ZERO_ZERO_M1], TOPOLOGY_DIM, coordX3N);
coordDesired[Y3]= coordX3N[Y3];
nList[NEB_P1_ZERO_M1]=gridComm.Get_cart_rank(coordDesired);
  /*********************/ 

gridComm.Get_coords(nList[NEB_M1_ZERO_ZERO], TOPOLOGY_DIM, coordX1N);
coordDesired[Y1]= coordX1N[Y1];
nList[NEB_M1_ZERO_M1]=gridComm.Get_cart_rank(coordDesired);
  /*********************/ 
gridComm.Get_coords(nList[NEB_ZERO_ZERO_P1], TOPOLOGY_DIM, coordX3N);
coordDesired[Y3]= coordX3N[Y3];
nList[NEB_M1_ZERO_P1]=gridComm.Get_cart_rank(coordDesired);
  /*********************/ 
int myGridRank;   
}


