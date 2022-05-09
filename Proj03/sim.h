#define J_GRIDSIZE 10
#define I_GRIDSIZE 10
#include <array>

using namespace std;

namespace sim{
  typedef array<double, 2> GridElT;
  typedef array<GridElT, J_GRIDSIZE> GridRowT;
  typedef array<GridRowT, I_GRIDSIZE> GridT;
  class Simulation{
    private:
      GridT grid;
      double t;
      const double END_T;
      const double DT;
      const double BOX_SPACING = 1;
      const int nRows = I_GRIDSIZE;
      const int nCols = J_GRIDSIZE;
      /*
       * workhorse of the simulation method. updates for a timestep
       * @param delT the time to update
       */
      void update(double delT);

      GridT* calcDiffYG();

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
      void printGrid();
  };
}
