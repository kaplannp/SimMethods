
#include <iostream>

#include <boost/numeric/odeint.hpp> // odeint function definitions

using namespace std;
using namespace boost::numeric::odeint;

// Defining a shorthand for the type of the mathematical state
typedef std::vector< double > state_type;

// System to be solved: dx/dt = -2 x
void my_system( const state_type &x , state_type &dxdt , const double t )
{
    dxdt[0] =  0*x[0]  +   1*x[1];
    dxdt[1] = -1*x[0]  - 2.2*x[1];
}

// Observer, prints time and state when called (during integration)
void my_observer( const state_type &x, const double t )
{
    std::cout  << t << "   " << x[0] << "   " << x[1] << std::endl;
}

// ------  Main
int main()
{
    state_type x0(2); // Initial condition, vector of 2 elements (position and velocity)
    x0[0] = 0.0;
    x0[1] = 1.0;

    // Integration parameters
    double t0 = 0.0;
    double t1 = 20.0;
    double dt = 1.0;

    // ----  Steppers definition  ----
    // Basic stepper:
    // follows given timestep size "dt"
    typedef runge_kutta4<state_type> rk4;

    // Error stepper, used to create the controlled stepper
    typedef runge_kutta_cash_karp54< state_type > rkck54;

    // Controlled stepper: 
    // it's built on an error stepper and allows us to have the output at each 
    // internally defined (refined) timestep, via integrate_adaptive call 
    typedef controlled_runge_kutta< rkck54 > ctrl_rkck54;

    // ----  Run integrations with the different steppers  ----

    // Run integrator with rk4 stepper
    std::cout << "==========  rk4 - basic stepper  ====================" << std::endl;
    integrate_adaptive( rk4(), my_system, x0, t0, t1, dt, my_observer );

    // Run integrator with controlled stepper
    std::cout << "==========  ctrl_rkck54 - Controlled Stepper  =======" << std::endl;
    x0[0] = 0.0; // Reset initial conditions
    x0[1] = 1.0;
    integrate_adaptive( ctrl_rkck54(), my_system, x0 , t0 , t1 , dt, my_observer );
}
