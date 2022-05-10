#define J_GRIDSIZE 10
#define I_GRIDSIZE 10
#include <array>

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
      const int nRows = I_GRIDSIZE;
      const int nCols = J_GRIDSIZE;
      const double BOX_SPACING = 1;
      const double VISCOSITY = 1;
      const double RHO = 1;
      /*
       * workhorse of the simulation method. updates for a timestep
       * @param delT the time to update
       */
      void update(double delT);

      /*
       * These methods are used to calculate changes in the flow per
       * units of space. DiffGI gives the differences along the y axis
       */
      GridT* calcDiffGI(GridT& grid);
      GridT* calcDiffGJ(GridT& grid);

      /*
       * These methods are used to calculate the double derivative
       * of flow per unit of space. Diff2GI gives the difference along the 
       * y axis
       */
      GridT* calcDiff2GI(GridT& grid);
      GridT* calcDiff2GJ(GridT& grid);

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
      GridT* calcFlowUpdate(GridT& grid1, GridT& grid2, GridT& diffI, GridT& diffJ, GridT& diff2I, 
          GridT& diff2J, double dt);

      /*
       * used to sum two grids together
       */
      GridT* sumGrids(GridT& first, GridT& second);

      /*
       * used to assign the values of one grid to another grid
       */
      void assignGrid(GridT& target, GridT& value);

      /*
       * Iteratively solve the Pressure Poisson distribution
       * @param uStar the starting value of uStar (just so we don't need to 
       *   recompute)
       * @param pressure so you have a good starting point
       * Big thanks to Gaurav. This function is a translation of his python
       * version
       *   https://towardsdatascience.com/computational-fluid-dynamics-using-\
       *   python-modeling-laminar-flow-272dad1ebec
       */
      void calcPressure(GridT& iStar, GridT& pressure);

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
      Simulation(double endT, double dt, double vel);
      /*
       * run the simulation
       */
      void simulate();

      /*
       * prints out the grid
       */
      void printGrid(GridT& iGrid, GridT& jGrid);
  };
}
