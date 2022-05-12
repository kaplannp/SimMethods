#include "sim.h"
#include <iostream>

using namespace std;
using namespace sim;

int main(){
  const double endT{1};
  const double dt{1e-3};
  const double vel{5};
  Simulation sim{endT, dt, vel, "simData.out"};
  sim.simulate();
  return 0;
}

