#include "sim.h"
#include <iostream>

using namespace std;
using namespace sim;

int main(){
  const double endT{10};
  const double dt{0.1};
  const double vel{5};
  Simulation sim{endT, dt, vel};
  sim.simulate();
  sim.printGrid();
}

