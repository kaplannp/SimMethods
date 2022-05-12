#include "sim.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <functional>
#include <fstream>

using namespace std;
namespace sim{
  Simulation::Simulation(double end_t, double dt, double vel, string fileName) :
    t{0}, END_T{end_t}, DT{dt}, fileStream(fileName), vel{vel}{
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
    for(int i = 0; i < nRows; i++)
      jGrid[i][0] = vel;
  }

  void Simulation::simulate(){
    cout << "simulating" << endl;
    while(t < END_T){
      cout << "t is " << t << endl;
      update(DT);
    }
  }

  void Simulation::update(double delT){
    cout << "----------------------------------------------" << endl;
    printGrid(iGrid, jGrid);
    GridT* iDiffI = calcDiffGI(iGrid);
    GridT* jDiffI = calcDiffGI(jGrid);
    GridT* iDiffJ = calcDiffGJ(iGrid);
    GridT* jDiffJ = calcDiffGJ(jGrid);
    GridT* iDiff2I = calcDiff2GI(iGrid);
    GridT* jDiff2I = calcDiff2GI(jGrid);
    GridT* iDiff2J = calcDiff2GJ(iGrid);
    GridT* jDiff2J = calcDiff2GJ(jGrid);
    GridT* iFlowUpdate = calcFlowUpdate(
        iGrid, jGrid, *iDiffI, *iDiffJ, *iDiff2I, *iDiff2J, delT);
    GridT* jFlowUpdate = calcFlowUpdate(
        jGrid, iGrid, *jDiffI, *jDiffJ, *jDiff2I, *jDiff2J, delT);
    GridT* iStar = sumGrids(iGrid, *iFlowUpdate);
    //applyIBoundaryConditions(*iStar);
    GridT* jStar = sumGrids(jGrid, *jFlowUpdate);
    GridT* updatePressure = calcPressure(*iStar, *jStar, pGrid, delT);
    //as of now spatial boundaries are junk, but pressure bounds are good stuff
    assignGrid(pGrid, *updatePressure);
    GridT* pDiffI = calcDiffGI(pGrid);
    GridT* pDiffJ = calcDiffGJ(pGrid);
    GridT* pDiffIMult = applyUnaryOp(
        *pDiffI, [=](double x){ return x * - delT / RHO;});
    GridT* pDiffJMult = applyUnaryOp(
        *pDiffJ, [=](double x){ return x * - delT / RHO;});
    GridT* updatedIGrid = sumGrids(*iStar, *pDiffIMult);
    GridT* updatedJGrid = sumGrids(*jStar, *pDiffJMult);
    //now you must apply boundary conditions to i and j grid
    applyIBoundaryConditions(*updatedIGrid);
    applyJBoundaryConditions(*updatedJGrid);
    assignGrid(iGrid, *updatedIGrid);
    assignGrid(jGrid, *updatedJGrid);
    delete iDiffI; iDiffI = nullptr;
    delete jDiffI; jDiffI = nullptr;
    delete iDiffJ; iDiffJ = nullptr;
    delete jDiffJ; jDiffJ = nullptr;
    delete iDiff2I; iDiff2I = nullptr;
    delete jDiff2I; jDiff2I = nullptr;
    delete iDiff2J; iDiff2J = nullptr;
    delete jDiff2J; jDiff2J = nullptr;
    delete iFlowUpdate; iFlowUpdate = nullptr;
    delete jFlowUpdate; jFlowUpdate = nullptr;
    delete iStar; iStar = nullptr;
    delete jStar; jStar = nullptr;
    delete pDiffI; pDiffI = nullptr;
    delete pDiffJ; pDiffJ = nullptr;
    delete updatedIGrid; updatedIGrid = nullptr;
    delete updatedJGrid; updatedJGrid = nullptr;
    delete pDiffIMult; pDiffIMult = nullptr;
    delete pDiffJMult; pDiffJMult = nullptr;
    delete updatePressure; updatePressure = nullptr;
    t += delT;
    //Log the current state
    fileStream << "---" << endl;
    printGrid(iGrid, fileStream);
    fileStream << "---" << endl;
    printGrid(jGrid, fileStream);
    fileStream << "---" << endl;
    printGrid(pGrid, fileStream);
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

  GridT* Simulation::subGrids(GridT& first, GridT& second){
    GridT* sum = new GridT();
    for(int i = 0; i < nRows; i++){
      for(int j = 0; j < nCols; j++){
        (*sum)[i][j] = first[i][j] - second[i][j];
      }
    }
    return sum;
  }

  GridT* Simulation::applyUnaryOp(GridT grid, function<double(double)> op){
    GridT* outGrid = new GridT();
    for(int i = 0; i < nCols; i++)
      for(int j = 0; j < nRows; j++)
        (*outGrid)[i][j] = op(grid[i][j]);
    return outGrid;
  }

  GridT* Simulation::calcPressure(GridT& iStar, GridT& jStar, GridT& pressure,
      double dt){
    GridT* retGrid = new GridT();
    assignGrid(*retGrid, pressure);
    //Define initial error and tolerance for convergence 
    //(error > tol necessary initially)
    double error = 1.0;
    double tol = 1e-3;

    //Evaluate derivative of starred velocities
    GridT* iStarDiffI = calcDiffGI(iStar);
    GridT* jStarDiffJ = calcDiffGJ(jStar);

    GridT pOld;

    //Continue iterative solution until error becomes smaller than tolerance
    int i = 0;
    while(error>tol){
      i++;
        
      //Save current pressure as p_old
      assignGrid(pOld,*retGrid);
      
      //Evaluate second derivative of pressure from p_old by doing this
      //p2_xy=(p_old[2:,1:cols+1]+p_old[0:rows,1:cols+1])/dy**2
      //     +(p_old[1:rows+1,2:]+p_old[1:rows+1,0:cols])/dx**2
      GridT* pDiff2I = new GridT();
      GridT* pDiff2J = new GridT();
      fillGrid(*pDiff2I, 0.0);
      fillGrid(*pDiff2J, 0.0);
      for(int i = 1; i < nRows - 1; i++){
        for(int j = 0; j < nCols; j++){
          (*pDiff2I)[i][j] = (pOld[i-1][j] + pOld[i+1][j]) / pow(BOX_SPACING,2);
        }
      }
      //changes in here
      for(int i = 0; i < nRows; i++){
        for(int j = 1; j < nCols - 1; j++){
          (*pDiff2J)[i][j] = (pOld[i][j-1] + pOld[i][j+1]) / pow(BOX_SPACING,2);
        }
      }
      //TODO Perhaps set pressure to 0 along boundary?

      GridT* pDiff2IJ = sumGrids(*pDiff2I, *pDiff2J);
      
      //Calculate new pressure doing this
      //p[1:rows+1,1:cols+1]=(p2_xy)*factor-(rho*factor/dt)*(ustar1_x+vstar1_y)
      double factor = pow(BOX_SPACING,2) / 4;
      for(int i = 0; i < nRows; i++){
        for(int j = 0; j < nCols; j++){
          (*retGrid)[i][j] = ((*pDiff2IJ)[i][j]) * factor 
            -((RHO*factor/dt) 
              *((*iStarDiffI)[i][j] + (*jStarDiffJ)[i][j])
             );
        }
      }

      //Find maximum error between old and new pressure matrices
      GridT* pDifference = subGrids(*retGrid, pOld);
      error = 0; // will be overwritten to be max of pDifference
      for(int i = 0; i < nRows; i++){
        for(int j = 0; j < nCols; j++){
          error = max(error, pow((*pDifference)[i][j],2));
        }
      }
      //Dunno what this is about 
      //Apply pressure boundary conditions
      applyPBoundaryConditions(*retGrid);
      //Escape condition in case solution does not converge after 500 iterations
      if(i>500){
        tol*=10;
      }
      delete pDiff2I; pDiff2I = nullptr;
      delete pDiff2J; pDiff2J = nullptr;
      delete pDiff2IJ; pDiff2IJ = nullptr;
      delete pDifference; pDifference = nullptr;
    }
    delete iStarDiffI; iStarDiffI = nullptr;
    delete jStarDiffJ; jStarDiffJ = nullptr;
    //printGrid(*retGrid);
    return retGrid;
  }

  void Simulation::fillGrid(GridT& grid, double val){
    for(GridRowT& row : grid)
      for(double& item : row)
        item = val;
  }

  void Simulation::assignGrid(GridT& target, GridT& value){
    for(int i = 0; i < nRows; i++){
      for(int j = 0; j < nCols; j++){
        target[i][j] = value[i][j];
      }
    }
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
      (*diffGrid)[i][0] = 0;
      (*diffGrid)[i][nCols-1] = 0;
      for(int j = 1; j < nCols - 1; j++){
        (*diffGrid)[i][j] = (grid[i][j+1] - grid[i][j-1]) 
          / (2*BOX_SPACING);
      }
    }
    return diffGrid;
  }

  GridT* Simulation::calcDiffGI(GridT& grid){
    GridT* diffGrid = new GridT();
    //now calculate the i component
    //you skip i=1 and because these values are assumed to be zero
    //TODO dunno if you should update the last one
    for(int j = 0; j < nCols; j++){
      (*diffGrid)[0][j] = 0;
      (*diffGrid)[nRows-1][j] = 0;
      for(int i = 1; i < nRows - 1; i++){
        (*diffGrid)[i][j] = (grid[i][j+1] - grid[i-1][j]) 
          / (2*BOX_SPACING);
      }
    }
    return diffGrid;
  }

  GridT* Simulation::calcDiff2GI(GridT& grid){
    GridT* diff2Grid = new GridT();
    for(int j = 0; j < nCols; j++){
      (*diff2Grid)[0][j] = 0;
      (*diff2Grid)[nRows-1][j] = 0;
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
      (*diff2Grid)[i][0] = 0;
      (*diff2Grid)[i][nCols-1] = 0;
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
    //these values are dependent on derivatives, so are undefined at bounds
    for(int j = 0; j < nCols; j++){
      (*update)[0][j] = 0;
      (*update)[nRows-1][j] = 0;
    }
    for(int i = 1; i < nRows - 1; i++){
      //set the 0th col to 0
      (*update)[i][0] = 0;
      for(int j = 1; j < nCols - 1; j++){
        //compute update 1
        (*update)[i][j] = dt * (
            -grid1[i][j] * diffJ[i][j]
            -grid2[i][j] * diffI[i][j]
            +VISCOSITY * (diff2I[i][j] + diff2J[i][j])
            );
      }
    }
    return update;
  }

  void Simulation::applyIBoundaryConditions(GridT& iGrid){
    for(int j = 0; j < nCols; j++){
      iGrid[0][j] = 0; //no ups and downs for top
      iGrid[nRows-1][j] = 0; //no ups and downs for bottom
    }
    for(int i = 0; i < nRows; i++){
      iGrid[i][0] = 0; //no ups and downs for left
      iGrid[i][nCols-1] = iGrid[i][nCols-2]; //Var deriv is 0 for right bound
    }
  }

  void Simulation::applyPBoundaryConditions(GridT& pGrid){
    for(int j = 0; j < nCols; j++){
      pGrid[0][j] = pGrid[1][j]; //top pressure deriv is 0
      pGrid[nRows-1][j] = pGrid[nRows-2][j]; //bottom pressure deriv is 0
    }
    for(int i = 0; i < nRows; i++){
      pGrid[i][0] = pGrid[i][1]; //left bound has const gradient
      pGrid[i][nCols-1] = 0; //right bound has no pressure
    }
  }

  void Simulation::applyJBoundaryConditions(GridT& iGrid){
    for(int j = 0; j < nCols; j++){
      iGrid[0][j] = 0; //no side to sides for top
      iGrid[nRows-1][j] = 0; //no side to sides for bottom
    }
    for(int i = 0; i < nRows; i++){
      iGrid[i][0] = vel; //const 5 coming from the left
      iGrid[i][nCols-1] = iGrid[i][nCols-2]; //Var deriv is 0 for right bound
    }
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
    
  void Simulation::printGrid(GridT& grid, ostream& out){
    for(int i = 0; i < nRows; i++){
      for(int j = 0; j < nCols; j++){
        out << round2Zero(grid[i][j]);
        //print a sep if you're not at end
        if(j != nCols - 1){
          out << "  ";
        }
      }
      out << endl;
    }
  }

  Simulation::~Simulation(){
    fileStream.close();
  }
}
