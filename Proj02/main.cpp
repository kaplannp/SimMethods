#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <functional>
#include <fstream>
#include <chrono>

#define N_SPECIES 3


using namespace std;

typedef double stateT[N_SPECIES];

/*
 * Global constants used in the simulation
 */
stateT X0 = {400.0, 0, 0};
double T0 = 0.0;
double ST = 0.001; //Sample time interval
double ENDT = .1;
double C[N_SPECIES] = {100, 50, 20};
int N_RUNS = 1000;


/*
 * Global variables used in the simulation
 */
stateT X;
double t;
int histSize = (int)(ENDT/ST);
stateT* history = new stateT[histSize];

/* Functions for reactions */
auto rxn0 = [](stateT& X){ 
  double newX0 = X[0] / 1.01f;
  X[1] += X[0] - newX0;
  X[0] = newX0;
};
auto rxn1 = [](stateT& X){ 
  double newX1 = X[1] / 1.01f;
  X[2] += X[1] - newX1;
  X[1] = newX1;
};
auto rxn2 = [](stateT& X){ 
  double newX2 = X[2] / 1.01f;
  X[2] = newX2;
};

vector<function<void(stateT&)>> applyReaction=vector<function<void(stateT&)>>();

typedef std::chrono::high_resolution_clock Clock;
Clock::time_point beginning = Clock::now();

/*
 * performs weighted random sampling
 * @param prop the propensities for each species
 * @param rand a random number between 0 and 1
 */
int wrs(double prop[N_SPECIES], double rand){
  double sum = 0;
  sum = accumulate(prop, prop+N_SPECIES, sum, plus<double>());
  double sample = sum*rand; //this is the random val between 0 and sum(prop)

  //Now search for the index closest to the input
  double sumSoFar = 0;
  int mu = 0;
  while(sumSoFar < sample){
    sumSoFar += prop[mu];
    mu++;
  }
  return --mu;
}


int nCalls = 0; //global for the number of calls to record State.
//This is used to index into history
/*
 * Records the current state
 * @param: X the state
 */
void recordState(const stateT& X){
  for(int i = 0; i < N_SPECIES; i++){
    history[nCalls][i] = X[i];
  }
  nCalls++;
}

/*
 * Used to calculate a propensity array given the state
 * @param X the current state
 * @returns: double[N_SPECIES] the propensities at each state
 */
double* calcProp(const stateT& X){
  double* props = new double[N_SPECIES];
  for(int i = 0; i < N_SPECIES; i++){
    props[i] = X[i] * C[i];
  }
  return props;
}



/**
 * This is unacceptable! C++ is going to make me write my own csv writing
 * function! I have issues! This is such a common and simple task that I would
 * like an stl function for parsing my csv, or at the very least a nice boost
 * reading writing class, but no!
 *
 * Preamble expounded, this function writes history in tsv format to filename
 * @param filename the name of the file to be written.
 *
 */
void writeHis(string filename){
  ofstream file;
  file.open(filename);
  for(int i = 0; i < histSize; i++){
    file << history[i][0];
    for(int j = 1; j < N_SPECIES; j++){
      file << "\t" <<  history[i][j];
    }
    file << endl;
  }
  file.close();
}

/*
 * Runs the gillespie algorithm recording state into global variable,
 * history
 * @param X0, the start state
 * @param tFinal the end time
 * @param t0 the start time
 */
void gillespie(const stateT& X0, double tFinal, double t0=0){
  nCalls=0; //reset nCalls for the history gathering
  //Fire up the random number generator
  default_random_engine generator; //NOTE: We are seeded!
  generator.seed((Clock::now() - beginning).count());
  cout << "running iteration with seed " << (Clock::now() - beginning).count() 
    <<endl;
  uniform_real_distribution<double> distribution(0.0,1.0);
  //Initialization
  copy(X0, X0+N_SPECIES, X);
  t = t0;
  double ns = t0; // next sample time

  //Algorithm
  while(t < tFinal){
    //cout << "iterating t = " << t << endl;
    double rand1 = distribution(generator);
    double* prop = calcProp(X);
    int mu = wrs(prop, rand1);
    double rand2 = distribution(generator);
    double sum = 0;
    double dt = (1/accumulate(prop, prop+N_SPECIES, sum)) * log(1/rand2);
    applyReaction[mu](X);
    t += dt;
    if(t > ns){
      recordState(X);
      ns += ST;
    }
    delete[] prop;
  }
  cout << "exiting with t = " << t << " and tfinal " << tFinal << endl;
}

void printHistory(){
  for(int i = 0; i < histSize; i++){
    for(int j = 0; j < N_SPECIES; j++){
      cout << history[i][j] << " ";
    }
    cout << endl;
  }
}

int main(){
  applyReaction.push_back(rxn0);
  applyReaction.push_back(rxn1);
  applyReaction.push_back(rxn2);
  for(int i = 0; i < N_RUNS; i++){
    gillespie(X0, ENDT, T0);
    //printHistory();
    writeHis("Out/out" + to_string(i) + ".tsv");
  }
}

