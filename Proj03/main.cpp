#include "sim.h"
#include <iostream>

using namespace std;
using namespace sim;

int main(){
  const double endT{.3};
  const double dt{0.001};
  const double vel{5};
  Simulation sim{endT, dt, vel, "simData.out"};
  sim.simulate();
  return 0;
}

