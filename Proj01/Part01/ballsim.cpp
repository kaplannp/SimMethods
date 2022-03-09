#include<vector>
#include<iostream>
#include<cmath>
#include <boost/numeric/odeint.hpp> // odeint function definitions

using namespace std;
using namespace boost::numeric::odeint;

typedef std::vector<double> stateType;

//The twiddlies
double startT = 0.0;
double endT = 1.5;
double firstDT = 0.05; // no longer used for timmestep version
double n_steps = 20000000.0;
//the constants
double dragC = 0.1;
double gravity = 9.8;
// px, py, vx, vy
stateType state = {0.0,10.0,1.0,2.0};

/*
 * drag coefficient is computed according to dragC
 * params:
 *   v: the velocity of the object
 * returns: a change in speed (magnitude)
 */
double windResist(const double v){
  return dragC*v*v;
}

/*
 * get the magnitude of a vector
 * params:
 *   v: the vector to sum
 * returns: double the l2 magnitude (not sqrt)
 */
double vecNorm(const vector<double>& v){
  double sum = 0;
  for(double d: v){
    sum += d*d;
  }
  return sum;
}

/*
 * params:
 *   double t: the timestep
 *   double y: the state
 */
// /*t*/ in func dec means take a nameless argument
void derivative(const stateType &x, stateType &dxdt, 
    const double /* t */){
  double velNorm2 = vecNorm({x[2], x[3]});
  double velNorm = sqrt(velNorm2);
  //Wind resistance on velocity
  dxdt[2] = -windResist(velNorm);
  dxdt[3] = dxdt[2];
  dxdt[2] *= -x[2] / velNorm;
  dxdt[3] *= -x[3] / velNorm;
  //gravity
  dxdt[3] -= gravity;

  //positions
  dxdt[0] = x[2];
  dxdt[1] = x[3];

  //gravity on pos
}

vector<stateType> storedData = vector<stateType>(n_steps);
/*
 * print out data
 */
void callbackData(const stateType &x , double /* t */){
  storedData.push_back(x);
  //for(double d: x)
  //  cout << d << ", ";
  //cout << endl;
}

int main(){
  double dt = (endT - startT) / n_steps;
  boost::numeric::odeint::integrate_n_steps(
      runge_kutta_cash_karp54<stateType>(),
      derivative, state, startT, dt, n_steps,
      callbackData);
  for(double d: storedData.back())
    cout << d << ", ";
  cout << endl;
  return 0;
}
