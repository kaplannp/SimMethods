#include <iostream>
#include <vector>
#include <random>
#include <cmath>

#define N_SPECIES 1


using namespace std;

typedef double stateT[N_SPECIES];

/*
 * Global constants used in the simulation
 */
stateT X0 = {4};
double T0 = 0.0;
double ST = 0.1; //Sample time interval
double ENDT = 10;
double C[N_SPECIES] = {.1};


/*
 * Global variables used in the simulation
 */
stateT X;
double t;
stateT* history = new stateT[(int)(ENDT/ST)];

/*
 * performs weighted random sampling
 * @param prop the propensities for each species
 * @param rand a random number between 0 and 1
 */
int wrs(double prop[N_SPECIES], double rand){
  double sum = 0;
  sum = accumulate(prop, prop+N_SPECIES, sum, plus<double>());
  double sample = sum*rand; //this is the random val between 0 and sum(prop)
  //cout << "sample is " << sample << endl;
  //cout << "prop is " << prop[0] << endl;

  //Now search for the index closest to the input
  double sumSoFar = 0;
  int mu = 0;
  while(sumSoFar < sample){
    sumSoFar += prop[mu];
    mu++;
  }
  return --mu;
}


/*
 * Records the current state
 * @param: X the state
 */
void recordState(stateT X){
  static int n_calls = 0;
  for(int i = 0; i < N_SPECIES; i++){
    history[n_calls][i] = X[i];
  }
  n_calls++;
}

/*
 * Used to calculate a propensity array given the state
 * @param X the current state
 * @returns: double[N_SPECIES] the propensities at each state
 */
double* calcProp(stateT X){
  double* props = new double[N_SPECIES];
  for(int i = 0; i < N_SPECIES; i++){
    props[i] = X[i] * C[i];
  }
  return props;
}

/*
 * Runs the gillespie algorithm
 */
void gillespie(stateT X0, double tFinal, double t0=0){
  //Fire up the random number generator
  default_random_engine generator; //NOTE: We are seeded!
  uniform_real_distribution<double> distribution(0.0,1.0);
  //Initialization
  copy(X0, X0+N_SPECIES, X);
  t = t0;
  double ns = t0; // next sample time

  //Algorithm
  while(t < tFinal){
    if(t > ns){
      recordState(X);
      ns += ST;
    }
    double rand1 = distribution(generator);
    double* prop = calcProp(X);
    int mu = wrs(prop, rand1);
    cout << "mu this time was " <<  mu << endl;
    double rand2 = distribution(generator);
    double sum = 0;
    double dt = (1/accumulate(prop, prop+N_SPECIES, sum)) * log(1/rand2);
    //x = applyReaction[mu](x)
    //t += tau
    t += .05;
    delete[] prop;
  }
}
int main(){
  gillespie(X0, ENDT, T0);
}

