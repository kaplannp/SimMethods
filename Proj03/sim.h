#define J_GRIDSIZE 10
#define I_GRIDSIZE 10
#include <array>
#include <functional>
#include <iostream>
#include <fstream>

using namespace std;

namespace sim{
  typedef array<double, J_GRIDSIZE> GridRowT;
  typedef array<GridRowT, I_GRIDSIZE> GridT;
  class Simulation{
    private:
      GridT iGrid;
      GridT jGrid;
      GridT pGrid;
      double t;
      const double END_T;
      const double DT;
      const double vel;
      const int nRows = I_GRIDSIZE;
      const int nCols = J_GRIDSIZE;
      const double BOX_SPACING = 1;
      const double VISCOSITY = 1e2;
      const double RHO = 1e-9;
      ofstream fileStream;

      //Bunch of space for computation
      //updates
      GridT* iDiffI;
      GridT* jDiffI;
      GridT* iDiffJ;
      GridT* jDiffJ;
      GridT* iDiff2I;
      GridT* jDiff2I;
      GridT* iDiff2J;
      GridT* jDiff2J;
      GridT* iFlowUpdate;
      GridT* jFlowUpdate;
      GridT* iStar;
      GridT* jStar;
      GridT* updatePressure;
      GridT* pDiffI;
      GridT* pDiffJ;
      GridT* pDiffIMult;
      GridT* pDiffJMult;
      GridT* updatedIGrid;
      GridT* updatedJGrid;
      GridT* iStarDiffI;
      GridT* jStarDiffJ;
      GridT* pDiff2I;
      GridT* pDiff2J;
      GridT* pDiff2IJ;
      GridT* pDifference;
      /*
       * workhorse;
       * @param delT the time to update
       */
      void update(double delT);


      //GOLDEN RULE of calcDiff is that you fill the undefined with 0
      /*
       * These methods are used to calculate changes in the flow per
       * units of space. DiffGI gives the differences along the y axis
       */
      void calcDiffGI(GridT& target, GridT& grid);
      void calcDiffGJ(GridT& target, GridT& grid);

      /*
       * These methods are used to calculate the double derivative
       * of flow per unit of space. Diff2GI gives the difference along the 
       * y axis
       */
      void calcDiff2GI(GridT& target, GridT& grid);
      void calcDiff2GJ(GridT& target, GridT& grid);

      /*
       * Calculate the increment of flow contributed by movement of flow.
       * Let the direction of the flows kept in grid1 be the x direction,
       * then this returns the flow update in the x direction
       * @param grid1 the flows in x direction
       * @param grid2 the flows in the y direction
       * @param diffI the 1st derivate of x in the i direction
       * @param diffJ the 1st derivative of x in the j direction
       * @param diff2I the 2nd derivative of x in the i direction
       * @param diff2J the 2nd derivative of x in the j direction
       * @param dt the change in time
       * 
       * @
       */
      void calcFlowUpdate(GridT& target, GridT& grid1, GridT& grid2, GridT& diffI, 
          GridT& diffJ, GridT& diff2I, GridT& diff2J, double dt);

      /*
       * used to sum two grids together
       */
      void sumGrids(GridT& target, GridT& first, GridT& second);
      /*
       * used to sub two grids together. It's about time we invented
       * a binary and unary operation applyer function or found one in the std,
       * but not about that at the moment
       */
      void subGrids(GridT& target, GridT& first, GridT& second);

      /*
       * I got fed up with writing new functions. Try this
       * @param GridT the grid
       * @param op a function taking in a double and returning a double
       *   to be applied to the grid
       */
      void applyUnaryOp(GridT& target, GridT& grid, 
          function<double(double)> op);

      /*
       * reflect from upper and lower any i current.
       * Set the i current for j = 0 to 0
       */
      void applyIBoundaryConditions(GridT& iGrid); 

      /*
       * set the boundary to zero everywhere inplace
       */
      void applyJBoundaryConditions(GridT& pGrid); 

      /*
       * set the upper and lower bound pressures equal to the pressure 
       * immediately below.
       *
       * do not modify pressures to right or left
       */
      void applyPBoundaryConditions(GridT& pGrid); 


      /*
       * Iteratively solve the Pressure Poisson distribution
       * and store the result to the pressure reference passed
       * @param uStar the starting value of uStar (just so we don't need to 
       *   recompute)
       * @param pressure so you have a good starting point
       * Big thanks to Gaurav. This function is a translation of his python
       * version
       *   https://towardsdatascience.com/computational-fluid-dynamics-using-\
       *   python-modeling-laminar-flow-272dad1ebec
       */
      void calcPressure(GridT& target, GridT& iStar, GridT& jStar, 
          GridT& pressure, double dt);

      /*
       * fills a grid with the value
       */
      void fillGrid(GridT& grid, double val);

      /* 
       * used to copy the values from value to target
       */
      void assignGrid(GridT& target, GridT& value);

      /*
       * calculates the norm of a 2d vector
       */
      double vecNorm(const array<double, 2>& vec);
      /*
       * useful for pretty printing
       */
      double round2Zero(double val);


    public:
      /*
       * @param endT the time at which to end the simulation
       * @param dt the time passing between timesteps
       */
      Simulation(double endT, double dt, double vel, string fileName);
      /*
       * run the simulation
       */
      void simulate();

      /*
       * prints out the grid
       */
      void printGrid(GridT& iGrid, GridT& jGrid);
      void printGrid(GridT& grid, ostream& out=cout);

      /*
       * closes the file
       */
      ~Simulation();
  };
}
