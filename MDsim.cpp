#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
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


// - - - - - IMPORTANT CONSTANTS - - - - -

const int PARTICLE_DIAM = 10;
const int NUM_COLUMNS = 2;
const int NUM_PARTICLES = 23892;
const long double WIDTH = 1918.6;
const long double HEIGHT = 1211.742744975187;

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
*/
long double** importData(const char* fileName) {
    ifstream file(fileName);

    // make space for the particle positions on the heap
    long double** initialParticles = new long double*[NUM_PARTICLES];
    for (int i = 0; i < NUM_PARTICLES; i++) {
        initialParticles[i] = new long double[NUM_COLUMNS];
    }

    // string* initialParticles = new string[23892][3];
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

// we need an array of vectors for the neighbor lists boi
array<vector<int>,NUM_PARTICLES> getNeighborsSimple(long double** particles) {
    array<vector<int>,NUM_PARTICLES> neighbors;
    // loop through each row in neighbors and check distances
    for (int i = 0; i < NUM_PARTICLES; ++i) {
        long double currentParticle [2] = {particles[i][0], particles[i][1]};

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
// is broken for some reason 
long double ** resolveCollisions(long double** particles, array<vector<int>,NUM_PARTICLES> neighbors) {

    for (int i = 0; i < NUM_PARTICLES; ++i) {
        int centralInd = i;
        long double centralParticle [2] = {particles[i][0], particles[i][1]};
        for (int j = 0; j < neighbors[i].size(); ++j) {
            int neighborInd = neighbors[i][j];
            long double neighborParticle [2] = {particles[neighborInd][0], particles[neighborInd][1]};

            // DO WRAPAROUND BOUNDARY CONDITIONS
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
                // cout << "here" << endl;
                long double r [2] = {((PARTICLE_DIAM-dist)/2)*(1/dist)*(centralParticle[0]-neighborParticle[0]), ((PARTICLE_DIAM-dist)/2)*(1/dist)*(centralParticle[1]-neighborParticle[1])};

                cout << r[0] << " " << r[1] << endl;

                particles[centralInd][0] += r[0];
                particles[centralInd][1] += r[1];
                particles[neighborInd][0] -= r[0];
                particles[neighborInd][1] -= r[1];

            }
        }

    }

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

long double ** runSim(long double ** particles, int numFrames) {
    int frameCount = 0;

    while (frameCount < numFrames) {
        array<vector<int>,NUM_PARTICLES> neighbors = getNeighborsSimple(particles);
        for (int i = 0; i < 10; ++i) {
            for (int m = 0; m < NUM_PARTICLES; ++m) {
                for (int n = 0; n < NUM_COLUMNS; ++n) {
                    particles[m][n] += norm(urbg)*SCALING_FACTOR;
                }
            }

            particles = resolveCollisions(particles, neighbors);
        }

        frameCount++;
        cout << "frame " + to_string(frameCount) + " done" << endl;
    }

    return particles;
}


int main() {
    // import particle positions 
    long double** particles = importData("initial_particles.txt");

    // run the sim
    particles = runSim(particles, 1);

    // export particle positions
    exportData(particles, "test.txt");

    return 0;
}