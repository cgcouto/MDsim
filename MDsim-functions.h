#ifndef MDSIM_FUNCTIONS
#define MDSIM_FUNCTIONS

using namespace std;

// - - - - - SIM PARAMETERS - - - - -

extern const int PARTICLE_DIAM;
extern const int NUM_PARTICLES;
extern const long double WIDTH;
extern const long double HEIGHT;

// - - - - - SIM FUNCTIONS - - - - -

long double** importData(const char* fileName);
long double ** doOneFrame(long double ** particles); 

#endif