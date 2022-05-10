#include "sim.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
namespace sim{
  Simulation::Simulation(double end_t, double dt, double vel) : t{0}, 
    END_T{end_t}, DT{dt}{
    cout.precision(2);
    for(GridRowT& row : grid)
      for(GridElT& el : row)
        for(double& val : el)
          val = 0.0;
    //initialize pressure grid
    for(PGridRowT& row : pGrid)
      for(double& val : row)
        val = 0.0;
    //update the flow of first column (constant)
    for(int i = 0; i < nCols; i++)
      grid[i][0][1] = vel;
  }

  void Simulation::simulate(){
    cout << "simulating" << endl;
    while(t < END_T){
      update(DT);
    }
  }

  void Simulation::update(double delT){
    GridT* diffI = calcDiffGI(grid);
    GridT* diffJ = calcDiffGJ(grid);
    GridT* diff2I = calcDiff2GI(grid);
    GridT* diff2J = calcDiff2GJ(grid);
    GridT* flowUpdate = calcFlowUpdate(*diffI, *diffJ, *diff2I, *diff2J, delT);
    GridT* uStar = sumGrids(grid, *flowUpdate);
    grid = *uStar;
    cout << "----------------------------------------------" << endl;
    printGrid(grid);
    delete diffI;
    delete diffJ;
    delete diff2I;
    delete diff2J;
    delete flowUpdate;
    delete uStar;
    t += delT;
  }

  GridT* Simulation::sumGrids(GridT& first, GridT& second){
    GridT* sum = new GridT();
    for(int i = 0; i < nRows; i++){
      for(int j = 0; j < nCols; j++){
        for(int k = 0; k < 2; k++){
          (*sum)[i][j][k] = first[i][j][k] + second[i][j][k];
        }
      }
    }
    return sum;
  }

  void Simulation::assignGrid(GridT& target, GridT& value){
    for(int i = 0; i < nRows; i++){
      for(int j = 0; j < nCols; j++){
        for(int k = 0; k < 2; k++){
          target[i][j][k] = value[i][j][k];
        }
      }
    }
  }

  void Simulation::calcPressure(GridT& uStar, PGridT& pressure){
    //Define initial error and tolerance for convergence (error > tol necessary initially)
    double error = 1.0;
    double tol = 1e-3;

    //Evaluate derivative of starred velocities
    GridT* uStarDiffI = calcDiffGI(uStar);
    GridT* uStarDiffJ = calcDiffGJ(uStar);

    //Continue iterative solution until error becomes smaller than tolerance
    //i=0
    //while(error>tol):
    //    i+=1
    //    
    //    #Save current pressure as p_old
    //    p_old=p.astype(float,copy=True)
    //    
    //    #Evaluate second derivative of pressure from p_old
    //    p2_xy=(p_old[2:,1:cols+1]+p_old[0:rows,1:cols+1])/dy**2+(p_old[1:rows+1,2:]+p_old[1:rows+1,0:cols])/dx**2
    //    
    //    #Calculate new pressure 
    //    p[1:rows+1,1:cols+1]=(p2_xy)*factor-(rho*factor/dt)*(ustar1_x+vstar1_y)
    //    
    //    #Find maximum error between old and new pressure matrices
    //    error=np.amax(abs(p-p_old))
    //    
    //    #Apply pressure boundary conditions
    //    SetPBoundary(space,left,right,top,bottom)
    //    
    //    #Escape condition in case solution does not converge after 500 iterations
    //    if(i>500):
    //        tol*=10
    delete uStarDiffI;
    delete uStarDiffJ;
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
        (*diffGrid)[i][j][1] = (grid[i][j+1][1] - grid[i][j-1][1]) 
          / (2*BOX_SPACING);
        (*diffGrid)[i][j][0] = (grid[i][j+1][0] - grid[i][j-1][0]) 
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
        (*diffGrid)[i][j][1] = (grid[i][j+1][1] - grid[i-1][j][1]) 
          / (2*BOX_SPACING);
        (*diffGrid)[i][j][0] = (grid[i][j+1][0] - grid[i-1][j][0]) 
          / (2*BOX_SPACING);
      }
    }
    return diffGrid;
  }

  GridT* Simulation::calcDiff2GI(GridT& grid){
    GridT* diff2Grid = new GridT();
    for(int j = 0; j < nCols; j++){
      for(int i = 1; i < nRows-1; i++){
        for(int k = 0; k < 2; k++){
          (*diff2Grid)[i][j][k] = (grid[i+1][j][k] - 2*grid[i][j][k] + 
              grid[i-1][j][k]) / pow(BOX_SPACING,2);
        }
      }
    }
    return diff2Grid;
  }

  GridT* Simulation::calcDiff2GJ(GridT& grid){
    GridT* diff2Grid = new GridT();
    for(int i = 0; i < nRows; i++){
      for(int j = 1; j < nCols-1; j++){
        for(int k = 0; k < 2; k++){
          (*diff2Grid)[i][j][k] = (grid[i][j+1][k] - 2*grid[i][j][k] + 
              grid[i][j-1][k]) / pow(BOX_SPACING,2);
        }
      }
    }
    return diff2Grid;
  }

  GridT* Simulation::calcFlowUpdate(GridT& diffI, GridT& diffJ, GridT& diff2I, 
          GridT& diff2J, double dt){
    GridT* update = new GridT();
    for(int i = 0; i < nRows; i++){
      //set the 0th col to 0
      (*update)[i][0][0] = 0;
      (*update)[i][0][1] = 0;
      for(int j = 1; j < nCols; j++){
        for(int k = 0; k < 2; k++){
          (*update)[i][j][k] = dt * (
              -grid[i][j][k] * diffJ[i][j][k]
              -grid[i][j][(k+1)%2] * diffI[i][j][k]
              +VISCOSITY * (diff2I[i][j][k] + diff2J[i][j][k])
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

  void Simulation::printGrid(GridT& grid){
    for(int i = 0; i < nRows; i++){
      for(int j = 0; j < nCols; j++){
        cout << setw(4) << "(" << round2Zero(grid[i][j][0]) << ", " << 
              round2Zero(grid[i][j][1]) << ")" << setw(4);
      }
      cout << endl;
    }

  }
}
