#include<vector>
#include<iostream>
#include <boost/numeric/odeint.hpp>
//TODO I don't really want all of boost, but this breaks and the above works
//#include <boost/numeric/odeint/integrate/XYZ.hpp>

using namespace std;

//The twiddlies
const double startT = 0.0;
const double endT = 0.6;
const double firstDT = 0.05;
// px, py, vx, vy
vector<double> state = {0.0,1.0,1.0,2.0};

/*
 * params:
 *   double t: the timestep
 *   double y: the state
 */
//TODO what the hell is this /* t */?
vector<double> derivative(const vector<double> &x, vector<double> &dxdt, const double /* t */){
//vector<double> derivative(const vector<double> &x, vector<double> &dxdt){
  //velocity = 
  return {5.9};
}

int main(){

  boost::numeric::odeint::integrate(derivative, state, startT, endT, firstDT);
  for (double s: state){
    std::cout << s << std::endl;
  }
  return 0;
}
