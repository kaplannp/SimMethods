#include "sim.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
namespace sim{
  Simulation::Simulation(double end_t, double dt, double vel) : t{0}, 
    END_T{end_t}, DT{dt}{
    cout.precision(2);
    //initialize the grids to zero
    for(GridRowT& row : iGrid)
      for(double& val : row)
        val = 0.0;
    for(GridRowT& row : jGrid)
      for(double& val : row)
        val = 0.0;
    for(GridRowT& row : pGrid)
      for(double& val : row)
        val = 0.0;
    //update the flow of first column (constant)
    for(int i = 0; i < nCols; i++)
      jGrid[i][0] = vel;
  }

  void Simulation::simulate(){
    cout << "simulating" << endl;
    while(t < END_T){
      update(DT);
    }
  }

  void Simulation::update(double delT){
    GridT* iDiffI = calcDiffGI(iGrid);
    GridT* jDiffI = calcDiffGI(jGrid);
    GridT* iDiffJ = calcDiffGJ(iGrid);
    GridT* jDiffJ = calcDiffGJ(jGrid);
    GridT* iDiff2I = calcDiff2GI(iGrid);
    GridT* jDiff2I = calcDiff2GI(jGrid);
    GridT* iDiff2J = calcDiff2GJ(iGrid);
    GridT* jDiff2J = calcDiff2GJ(jGrid);
    //TODO flow updates
    GridT* iFlowUpdate = calcFlowUpdate(iGrid, jGrid, *iDiffI, *iDiffJ, *iDiff2I,
        *iDiff2J, delT);
    GridT* jFlowUpdate = calcFlowUpdate(jGrid, iGrid, *jDiffI, *jDiffJ, *jDiff2I,
        *jDiff2J, delT);
    GridT* iStar = sumGrids(iGrid, *iFlowUpdate);
    GridT* jStar = sumGrids(jGrid, *jFlowUpdate);
    iGrid = *iStar;
    jGrid = *jStar;
    cout << "----------------------------------------------" << endl;
    printGrid(iGrid, jGrid);
    delete iDiffI;
    delete jDiffI;
    delete iDiffJ;
    delete jDiffJ;
    delete iDiff2I;
    delete jDiff2I;
    delete iDiff2J;
    delete jDiff2J;
    delete iFlowUpdate;
    delete jFlowUpdate;
    delete iStar;
    delete jStar;
    t += delT;
  }

  GridT* Simulation::sumGrids(GridT& first, GridT& second){
    GridT* sum = new GridT();
    for(int i = 0; i < nRows; i++){
      for(int j = 0; j < nCols; j++){
        (*sum)[i][j] = first[i][j] + second[i][j];
      }
    }
    return sum;
  }

  void Simulation::calcPressure(GridT& iStar, GridT& pressure){
    //Define initial error and tolerance for convergence (error > tol necessary initially)
    double error = 1.0;
    double tol = 1e-3;

    //Evaluate derivative of starred velocities
    GridT* iStarDiffI = calcDiffGI(iStar);
    GridT* iStarDiffJ = calcDiffGJ(iStar);

    GridT pOld;

    //Continue iterative solution until error becomes smaller than tolerance
    int i = 0;
    while(error>tol){
      i++;
        
      //Save current pressure as p_old
      pOld = pressure;
      
      //Evaluate second derivative of pressure from p_old
      //p2_xy=(p_old[2:,1:cols+1]+p_old[0:rows,1:cols+1])/dy**2+(p_old[1:rows+1,2:]+p_old[1:rows+1,0:cols])/dx**2
    //  
    //  #Calculate new pressure 
    //  p[1:rows+1,1:cols+1]=(p2_xy)*factor-(rho*factor/dt)*(ustar1_x+vstar1_y)
    //  
    //  #Find maximum error between old and new pressure matrices
    //  error=np.amax(abs(p-p_old))
    //  
    //  #Apply pressure boundary conditions
    //  SetPBoundary(space,left,right,top,bottom)
    //  
    //  #Escape condition in case solution does not converge after 500 iterations
      if(i>500){
        tol*=10;
      }
    }
    delete iStarDiffI;
    delete iStarDiffJ;
  }


  //Math helpers
  GridT* Simulation::calcDiffGJ(GridT& grid){
    GridT* diffGrid = new GridT();
    //skip the first row b/c its
    //Note you start by calculating the flow out of last, but you'll bounce it
    //back later
    //This calculates j component
    //Assumption: the fluid never bounces back against
    for(int i = 0; i < nRows; i++){
      //go until the last col. skip first because you're not updating it. Ever.
      for(int j = 1; j < nCols - 1; j++){
        (*diffGrid)[i][j] = (grid[i][j+1] - grid[i][j-1]) 
          / (2*BOX_SPACING);
      }
    }
    ////Now Reflect TODO you may find this useful later in life
    //for(int j = 0; j < nCols; j++){
    //  grid[1][j][0] = abs(grid[1][j][0]); //positive
    //  grid[nRows-2][j][0] = -abs(grid[1][j][0]); //neg
    //}
    return diffGrid;
  }

  GridT* Simulation::calcDiffGI(GridT& grid){
    GridT* diffGrid = new GridT();
    //now calculate the i component
    //you skip i=1 and because these values are assumed to be zero
    //TODO dunno if you should update the last one
    for(int j = 0; j < nCols; j++){
      for(int i = 1; i < nRows; i++){
        (*diffGrid)[i][j] = (grid[i][j+1] - grid[i-1][j]) 
          / (2*BOX_SPACING);
      }
    }
    return diffGrid;
  }

  GridT* Simulation::calcDiff2GI(GridT& grid){
    GridT* diff2Grid = new GridT();
    for(int j = 0; j < nCols; j++){
      for(int i = 1; i < nRows-1; i++){
        (*diff2Grid)[i][j] = (grid[i+1][j] - 2*grid[i][j] + 
            grid[i-1][j]) / pow(BOX_SPACING,2);
      }
    }
    return diff2Grid;
  }

  GridT* Simulation::calcDiff2GJ(GridT& grid){
    GridT* diff2Grid = new GridT();
    for(int i = 0; i < nRows; i++){
      for(int j = 1; j < nCols-1; j++){
        (*diff2Grid)[i][j] = (grid[i][j+1] - 2*grid[i][j] + 
            grid[i][j-1]) / pow(BOX_SPACING,2);
      }
    }
    return diff2Grid;
  }

  GridT* Simulation::calcFlowUpdate(GridT& grid1, GridT& grid2, 
      GridT& diffI, GridT& diffJ, GridT& diff2I, GridT& diff2J, double dt){
    GridT* update = new GridT();
    for(int i = 0; i < nRows; i++){
      //set the 0th col to 0
      (*update)[i][0] = 0;
      (*update)[i][0] = 0;
      for(int j = 1; j < nCols; j++){
        for(int k = 0; k < 2; k++){
          //compute update 1
          (*update)[i][j] = dt * (
              -grid1[i][j] * diffJ[i][j]
              -grid2[i][j] * diffI[i][j]
              +VISCOSITY * (diff2I[i][j] + diff2J[i][j])
              );
        }
      }
    }
    return update;
  }

  double Simulation::vecNorm(const array<double, 2>& vec){
    double sum = 0;
    for(double x : vec)
      sum += pow(x,2);
    return sqrt(sum);
  }

  double Simulation::round2Zero(double val){
    return (abs(val) < 1e-9 ) ? 0 : val;
  }

  void Simulation::printGrid(GridT& iGrid, GridT& jGrid){
    for(int i = 0; i < nRows; i++){
      for(int j = 0; j < nCols; j++){
        cout << setw(4) << "(" << round2Zero(iGrid[i][j]) << ", " << 
              round2Zero(jGrid[i][j]) << ")" << setw(4);
      }
      cout << endl;
    }

  }
}
