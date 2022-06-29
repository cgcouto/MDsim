#define _USE_MATH_DEFINES

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

// packages i will probably need in the future

// OpenCV for visualization
// some kind of graphing package for making plots
// a package for delaunay triangulation/voronoi


// things i will need to write in the future

// makeHCP
// initializeGrains
// getNeighborAndPsi6
// getImage
// getPsi6Image


using namespace std; // so we don't have to put std:: in front of certain things


// - - - - - SIM PARAMETERS - - - - -

const int PARTICLE_DIAM = 10;
const int NUM_COLUMNS = 2;
const int NUM_PARTICLES = 23892; // these three won't be hardcoded once lattice generation is ported
const long double WIDTH = 1918.6;
const long double HEIGHT = 1211.742744975187;

// - - - - - IMPORTANT CONSTANTS - - - - -

const long double STEP_SIZE = 0.0050;
const long double VISCOSITY = 1.99e-3;
const int TEMPERATURE = 300;
const long double BOLTZ_CONSTANT = 1.381e-23;
const long double SCALING_FACTOR = PARTICLE_DIAM/(1.3e-6);
const long double  SQUIG = (6*M_PI*(PARTICLE_DIAM/2))/SCALING_FACTOR;
const long double NEIGHB_THRESHOLD = 2;

const int MU = 0;
const long double SIGMA = sqrt((2*BOLTZ_CONSTANT*TEMPERATURE)/(SQUIG * STEP_SIZE));
auto urbg = mt19937(123); // generator(seed) - mersenne twister based around 2^(19937)-1 
auto norm = normal_distribution <long double> (MU,SIGMA);

/*
    Translates a MATLAB 2D particle array (that is saved via writematrix)
    into a C++ array
    INPUTS:
            fileName: points to the name of your desired file. must include file extension
    OUTPUTS:
            initialParticles: the particle positions in a 2d array format 
*/
long double** importData(const char* fileName) {
    ifstream file(fileName);

    // make space for the particle positions on the heap
    long double** initialParticles = new long double*[NUM_PARTICLES];
    for (int i = 0; i < NUM_PARTICLES; i++) {
        initialParticles[i] = new long double[NUM_COLUMNS];
    }

    // get all the data into our array
    if (file.is_open()) {
        string data = "";
        int count = 0;
        while(getline(file, data)) {   // get a whole line
            stringstream ss(data);
            while(getline(ss, data, ',')) { // get each element in the line
                initialParticles[count/NUM_COLUMNS][count % NUM_COLUMNS] = stold(data);
                // cout << to_string(initialParticles[count/NUM_COLUMNS][count % NUM_COLUMNS]) << endl;
                count++;
            }
        }
    }
    file.close();

    return initialParticles;
}

/*
    Converts the given C++ array into a txt file with the given name
    INPUTS:
            particles: 2d array containing particle positions
            saveName: points to the name of the file you want to save to. must include file extension
*/
void exportData(long double** particles, const char* saveName) {
    // open desired file for saving
    ofstream myfile (saveName);

    if (myfile.is_open()) {
        // add the array data to the file in the proper format
        for (int i = 0; i < NUM_PARTICLES; ++i) {
            for (int j = 0; j < NUM_COLUMNS-1; ++j) {
                myfile << to_string(particles[i][j]);
                myfile << " "; 
            }
            myfile << std::fixed << particles[i][NUM_COLUMNS-1];
            myfile << "\n";
        }
        myfile.close();
  }
    // clear particle positions from heap
    for (int i = 0; i < NUM_PARTICLES; i++) {
        delete[] particles[i];
    }
    delete[] particles;
}

/*
    Creates a neighbor list for the particles using simple means (i.e. not delaunay/voronoi)
    INPUTS:
            particles: 2d array containing particle positions
    OUTPUTS:
            neighbors: an array of vectors containing relevant array indexes from particles
*/
array<vector<int>,NUM_PARTICLES> getNeighborsSimple(long double** particles) {
    array<vector<int>,NUM_PARTICLES> neighbors;
    // loop through each row in neighbors and check distances
    for (int i = 0; i < NUM_PARTICLES; ++i) {
        long double currentParticle [2] = {particles[i][0], particles[i][1]};

        // form the relevant delta around currentParticle
        long double xMin = currentParticle[0]-PARTICLE_DIAM*NEIGHB_THRESHOLD;
        long double xMax = currentParticle[0]+PARTICLE_DIAM*NEIGHB_THRESHOLD;
        long double yMin = currentParticle[1]-PARTICLE_DIAM*NEIGHB_THRESHOLD;
        long double yMax = currentParticle[1]+PARTICLE_DIAM*NEIGHB_THRESHOLD;

        if (xMin < 0) {
            xMin += WIDTH;
        }

        if (xMax > WIDTH) {
            xMax -= WIDTH;
        }

        if (yMin < 0) {
            yMin += HEIGHT;
        }

        if (yMax > HEIGHT) {
            yMax -= HEIGHT;
        }

        for (int j = 0; j < NUM_PARTICLES; ++j) {
            if (i != j) {
                long double testParticle [2] = {particles[j][0], particles[j][1]};
                if (testParticle[0] > xMin && testParticle[0] < xMax && testParticle[1] > yMin && testParticle[1] < yMax) {
                    neighbors[i].push_back(j);
                }
            }
        }
        // if (neighbors[i].size() == 0) {
        //     cout << particles[i][0] << " " << particles[i][1] << endl;
        // }
        // cout << neighbors[i].size() << endl;
    }

    // cout << PARTICLE_DIAM*NEIGHB_THRESHOLD << endl;
    // cout << WIDTH-PARTICLE_DIAM*NEIGHB_THRESHOLD << endl;
    // cout << HEIGHT-PARTICLE_DIAM*NEIGHB_THRESHOLD << endl;

    return neighbors;
}

/*
    Uses the neighbor list to resolve particle collisions
    INPUTS:
            particles: 2d array containing particle positions
            neighbors: array of vectors containing particle indices
    OUTPUTS:
            particles: same as input, just with collisions resolved
*/
long double ** resolveCollisions(long double** particles, array<vector<int>,NUM_PARTICLES> neighbors) {

    // loop through the neighbor list to get relevant particle inds
    for (int i = 0; i < NUM_PARTICLES; ++i) {
        int centralInd = i;
        long double centralParticle [2] = {particles[i][0], particles[i][1]};
        for (int j = 0; j < neighbors[i].size(); ++j) {
            int neighborInd = neighbors[i][j];
            long double neighborParticle [2] = {particles[neighborInd][0], particles[neighborInd][1]};

            // adjust particles based on wraparound boundary conditions
            if (neighborParticle[0] < PARTICLE_DIAM && centralParticle[0] > WIDTH-PARTICLE_DIAM) {
                neighborParticle[0] += WIDTH;
            } else if (neighborParticle[0] > WIDTH-PARTICLE_DIAM && centralParticle[0] < PARTICLE_DIAM) {
                neighborParticle[0] -= WIDTH;
            }

            if (neighborParticle[1] < PARTICLE_DIAM && centralParticle[1] > HEIGHT-PARTICLE_DIAM) {
                neighborParticle[1] += HEIGHT;
            } else if (neighborParticle[1] > HEIGHT-PARTICLE_DIAM && centralParticle[1] < PARTICLE_DIAM) {
                neighborParticle[1] -= HEIGHT;
            }

            // get distance between points
            long double dist = sqrt(pow(centralParticle[0]-neighborParticle[0],2) + pow(centralParticle[1]-neighborParticle[1],2));

            if (dist < PARTICLE_DIAM) {
                // PUSH THEM BACK
                long double r [2] = {((PARTICLE_DIAM-dist)/2)*(1/dist)*(centralParticle[0]-neighborParticle[0]), ((PARTICLE_DIAM-dist)/2)*(1/dist)*(centralParticle[1]-neighborParticle[1])};

                // cout << r[0] << " " << r[1] << endl;

                particles[centralInd][0] += r[0];
                particles[centralInd][1] += r[1];
                particles[neighborInd][0] -= r[0];
                particles[neighborInd][1] -= r[1];

            }
        }

    }

    // enforce wraparound boundary conditions
    for (int i = 0; i < NUM_PARTICLES; ++i) {
        if (particles[i][0] < 0) {
            particles[i][0] += WIDTH;
        } else if (particles[i][0] > WIDTH) {
            particles[i][0] -= WIDTH;
        }

        if (particles[i][1] < 0) {
            particles[i][1] += HEIGHT;
        } else if (particles[i][1] > HEIGHT) {
            particles[i][1] -= HEIGHT;
        }

    }
    return particles;
}

/*
    Runs the molecular dynamics simulation
    INPUTS:
            particles: 2d array containing particle positions
            numFrames: how long you want the simulation to run for
    OUTPUTS:
            particles: the final particle positions
*/
long double ** runSim(long double ** particles, int numFrames) {
    int frameCount = 0;

    while (frameCount < numFrames) {
        // get neighbor lists 
        array<vector<int>,NUM_PARTICLES> neighbors = getNeighborsSimple(particles);
        
        for (int i = 0; i < 10; ++i) {
            // move all the particles via brownian motion
            for (int m = 0; m < NUM_PARTICLES; ++m) {
                for (int n = 0; n < NUM_COLUMNS; ++n) {
                    long double test = norm(urbg)*SCALING_FACTOR;
                    particles[m][n] += test;
                    // cout << test <<endl;
                }
            }
            // do collision resolution
            particles = resolveCollisions(particles, neighbors);
        }

        // report the frame as done
        frameCount++;
        cout << "frame " + to_string(frameCount) + " done" << endl;
    }

    return particles;
}


int main() {
    // import particle positions 
    long double** particles = importData("initial_particles.txt");

    // run the sim
    particles = runSim(particles, 10);

    // export particle positions
    exportData(particles, "test.txt");

    return 0;
}