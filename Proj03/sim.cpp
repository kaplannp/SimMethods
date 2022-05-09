#include "sim.h"
#include <iostream>

using namespace std;
namespace sim{
  Simulation::Simulation(double end_t, double dt, double vel) : t{0}, END_T{end_t}, DT{dt}{
    cout.precision(2);
    for(GridRowT col : grid)
      for(GridElT el : col)
        el.fill(0);
    //update the flow of first column (constant)
    int j = 0;
    for(int i = 0; i < nCols; i++)
      grid[i][j][1] = vel;
  }

  void Simulation::simulate(){
    cout << "simulating" << endl;
    while(t < END_T){
      update(DT);
    }
  }

  void Simulation::update(double delT){
    t += delT;
  }

  //Math helpers
  GridT* Simulation::calcDiffYG(){
    //TODO sorta. Not quite there
    GridT* diffGrid = new GridT();
    //skip the first row b/c its
    //Note you start by calculating the flow out of last, but you'll bounce it
    //back later
    for(int i = 1; i < nRows; i++){
      GridRowT row = grid[i];
      //go until the last col. skip first because you're not updating it. Ever.
      for(int j = 1; j < nCols; j++){
        row[j][1] = (row[j][1] - row[j-1][1]) / BOX_SPACING;
      }
    }
    return diffGrid;
  }

  void Simulation::printGrid(){
    for(int i = 0; i < nRows; i++){
      cout << endl;
      for(int j = 0; j < nCols; j++){
        cout << "(" << grid[i][j][0] << ", " << grid[i][j][1] << ")" << '\t';
      }
    }

  }
}
