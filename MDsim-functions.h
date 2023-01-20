#ifndef MDSIM_FUNCTIONS
#define MDSIM_FUNCTIONS

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <random>
#include <vector>
#include <array>
#include <chrono>

using namespace std;

// - - - - - SIM PARAMETERS - - - - -

extern const int PARTICLE_DIAM;
extern const int NUM_COLUMNS;
extern const int NUM_PARTICLES;
extern const long double WIDTH;
extern const long double HEIGHT;

// - - - - - IMPORTANT CONSTANTS - - - - -

extern const long double STEP_SIZE;
extern const long double VISCOSITY;
extern const int TEMPERATURE;
extern const long double BOLTZ_CONSTANT;
extern const long double SCALING_FACTOR;
extern const long double  SQUIG;
extern const long double NEIGHB_THRESHOLD;

extern const int MU;
extern const long double SIGMA;
extern random_device rd;
extern unsigned int seed;
// extern mt19937 gen(seed);
// extern normal_distribution <float> d(MU,SIGMA);

extern const int TABLE_WIDTH;
extern const int TABLE_HEIGHT;
extern const int TABLE_SIZE;


extern const long double CELL_WIDTH;
extern const long double CELL_HEIGHT;


long double** importData(const char* fileName);

void exportData(long double** particles, int frameNum, const char* saveName);

void exportNeighbors(array<vector<int>,5700> neighbors, const char* saveName);

void exportMapping(int** neighborMapping, const char* saveName);

int ** buildNeighborMapping();

array<vector<int>,5700> getNeighborsCell(long double** particles, int** neighborMapping);

long double ** resolveCollisionsCell(long double** particles, array<vector<int>,5700> neighbors, int** neighborMapping);

long double ** runSimCell(long double ** particles, int numFrames, int savingFrequency, const char * saveName); 

long double ** doOneFrame(long double ** particles, const char * saveName); 

#endif