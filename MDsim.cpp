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

// things i will need to fix in the future

// the sig fig problem with exportData

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
const long double  SQUIG = (6*M_PI*VISCOSITY*(PARTICLE_DIAM/2))/SCALING_FACTOR;
const long double NEIGHB_THRESHOLD = 2;

const int MU = 0;
const long double SIGMA = sqrt((2*BOLTZ_CONSTANT*TEMPERATURE)/SQUIG * STEP_SIZE);
random_device rd; // creates random integer seed values each time we run
unsigned int seed = rd();
mt19937 gen(seed); // generator(seed) - mersenne twister based around 2^(19937)-1 
normal_distribution <float> d(MU,SIGMA);

// sizes the neighbor cells variably (you'll get slightly different result depending on the sim width/height)
// should consult with sharon next week before finalizing this...
const int TABLE_WIDTH = static_cast<int> (floor(WIDTH/(PARTICLE_DIAM*NEIGHB_THRESHOLD)));
const int TABLE_HEIGHT = static_cast<int> (floor(HEIGHT/(PARTICLE_DIAM*NEIGHB_THRESHOLD)));
const int TABLE_SIZE = 95*60;

const long double CELL_WIDTH = WIDTH/TABLE_WIDTH;
const long double CELL_HEIGHT = HEIGHT/TABLE_HEIGHT;


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
    // for (int i = 0; i < NUM_PARTICLES; i++) {
    //     delete[] particles[i];
    // }
    // delete[] particles;
}


int ** buildNeighborMapping() {
    int** neighborMapping = new int*[TABLE_SIZE];
    for (int i = 0; i < TABLE_SIZE; ++i) {
        neighborMapping[i] = new int[8];
    }

    // loop through all the possible cell coordinates
    for (int m = 0; m < TABLE_WIDTH; ++m) {
        for (int n = 0; n < TABLE_HEIGHT; ++n) {

            int tableIndex = m + TABLE_WIDTH * n;
            // do our mapping and encode it in indices
            // the ninth neighbor (itself) is implied by the tableIndex
            neighborMapping[tableIndex][0] = ((TABLE_WIDTH + ((m-1))) % TABLE_WIDTH) + TABLE_WIDTH*n;
            neighborMapping[tableIndex][1] = (((TABLE_WIDTH + ((m+1))) % TABLE_WIDTH) % TABLE_WIDTH) + TABLE_WIDTH*n;
            neighborMapping[tableIndex][2] = m + TABLE_WIDTH*(((TABLE_HEIGHT + ((n-1))) % TABLE_HEIGHT) % TABLE_HEIGHT);
            neighborMapping[tableIndex][3] = m + TABLE_WIDTH*(((TABLE_HEIGHT + ((n+1))) % TABLE_HEIGHT) % TABLE_HEIGHT);
            neighborMapping[tableIndex][4] = (((TABLE_WIDTH + ((m-1))) % TABLE_WIDTH) % TABLE_WIDTH) + TABLE_WIDTH*(((TABLE_HEIGHT + ((n-1))) % TABLE_HEIGHT) % TABLE_HEIGHT);
            neighborMapping[tableIndex][5] = (((TABLE_WIDTH + ((m-1))) % TABLE_WIDTH) % TABLE_WIDTH) + TABLE_WIDTH*(((TABLE_HEIGHT + ((n+1))) % TABLE_HEIGHT) % TABLE_HEIGHT);
            neighborMapping[tableIndex][6] = (((TABLE_WIDTH + ((m+1))) % TABLE_WIDTH) % TABLE_WIDTH) + TABLE_WIDTH*(((TABLE_HEIGHT + ((n-1))) % TABLE_HEIGHT) % TABLE_HEIGHT);
            neighborMapping[tableIndex][7] = (((TABLE_WIDTH + ((m+1))) % TABLE_WIDTH) % TABLE_WIDTH) + TABLE_WIDTH*(((TABLE_HEIGHT + ((n+1))) % TABLE_HEIGHT) % TABLE_HEIGHT);

            int test = 1000;
            if (tableIndex == test) {
                std::cout << neighborMapping[test][0] << endl;
                std::cout << neighborMapping[test][1] << endl;
                std::cout << neighborMapping[test][2] << endl;
                std::cout << neighborMapping[test][3] << endl;
                std::cout << neighborMapping[test][4] << endl;
                std::cout << neighborMapping[test][5] << endl;
                std::cout << neighborMapping[test][6] << endl;
                std::cout << neighborMapping[test][7] << endl;
            }
        }

    }

    return neighborMapping;
}

array<vector<int>,TABLE_SIZE> getNeighborsCell(long double** particles, int** neighborMapping) {
    array<vector<int>,TABLE_SIZE> neighbors;

    // int arrMax = 0;
    // int sumMax = 0;
    for (int i = 0; i < NUM_PARTICLES; ++i) {
        // get the position of the cell the particle is in
        int tableX = static_cast<int> (floor(particles[i][0] / CELL_WIDTH));
        int tableY = static_cast<int> (floor(particles[i][1] / CELL_HEIGHT));
        
        // translate position to index and add to vector
        // cout << tableX + TABLE_HEIGHT*tableY << " " << particles[i][0] << " " << particles[i][1] << endl;
        neighbors[tableX + TABLE_WIDTH*tableY].push_back(i);

        // if (tableX + TABLE_WIDTH*tableY > arrMax) {
        //     arrMax = tableX + TABLE_WIDTH*tableY;
        // }
        // if (tableX + tableY > sumMax) {
        //     sumMax = tableX+tableY;
        // }

    }

    // cout << TABLE_SIZE << endl;
    // cout << arrMax << endl;
    // cout << sumMax << endl;
    return neighbors;
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
        for (int j = 0; j < NUM_PARTICLES; ++j) {
            if (i != j) {
                long double currentParticle [2] = {particles[i][0], particles[i][1]};
                long double testParticle [2] = {particles[j][0], particles[j][1]};

                // adjust particles based on wraparound boundary conditions
                if (testParticle[0] < PARTICLE_DIAM*NEIGHB_THRESHOLD && currentParticle[0] > WIDTH-PARTICLE_DIAM*NEIGHB_THRESHOLD ) {
                    testParticle[0] += WIDTH;
                } else if (testParticle[0] > WIDTH-PARTICLE_DIAM*NEIGHB_THRESHOLD  && currentParticle[0] < PARTICLE_DIAM*NEIGHB_THRESHOLD ) {
                    testParticle[0] -= WIDTH;
                }

                if (testParticle[1] < PARTICLE_DIAM*NEIGHB_THRESHOLD  && currentParticle[1] > HEIGHT-PARTICLE_DIAM*NEIGHB_THRESHOLD ) {
                    testParticle[1] += HEIGHT;
                } else if (testParticle[1] > HEIGHT-PARTICLE_DIAM*NEIGHB_THRESHOLD  && currentParticle[1] < PARTICLE_DIAM*NEIGHB_THRESHOLD ) {
                    testParticle[1] -= HEIGHT;
                }

                long double dist = sqrt(pow(currentParticle[0]-testParticle[0],2) + pow(currentParticle[1]-testParticle[1],2));

                if (dist < PARTICLE_DIAM*NEIGHB_THRESHOLD) {
                    neighbors[i].push_back(j);
                }
                    
            }
    
        }
    }

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
long double ** resolveCollisionsSimple(long double** particles, array<vector<int>,NUM_PARTICLES> neighbors) {

    int collisionCount = 0;
    // loop through the neighbor list to get relevant particle inds
    for (int i = 0; i < NUM_PARTICLES; ++i) {
        int centralInd = i;
        long double centralParticle [2] = {particles[i][0], particles[i][1]};
        cout << neighbors[i].size() << endl;
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
                collisionCount++;
                // PUSH THEM BACK
                long double r [2] = {((PARTICLE_DIAM-dist)/2)*(1/dist)*(centralParticle[0]-neighborParticle[0]), ((PARTICLE_DIAM-dist)/2)*(1/dist)*(centralParticle[1]-neighborParticle[1])};

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
    cout << collisionCount << endl;
    return particles;
}

long double ** resolveCollisionsCell(long double** particles, array<vector<int>,TABLE_SIZE> neighbors, int** neighborMapping) {
    
    // go through all the indices in the neighbor list (size is equal to NUM_PARTICLES)

    int collisionCount = 0;
    for (int i = 0; i < TABLE_SIZE; ++i) {
        for (int j = 0; j < neighbors[i].size(); ++j) {
            int centralInd = neighbors[i][j];
            long double centralParticle [2] = {particles[centralInd][0], particles[centralInd][1]};

            // use the neighbor mapping to get all the neighbor particles in the appropriate cells
            for (int m = 0; m < 9; ++m) {
                // make sure we catch the cell the particle is in as well!
                int currentCell;
                if (m == 8) {
                    currentCell = i;
                } else {
                    currentCell = neighborMapping[i][m];
                    
                }

                for (int n = 0; n < neighbors[currentCell].size(); ++n) {
                    int neighborInd = neighbors[currentCell][n];

                    if (neighborInd != centralInd) {
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
                        // std::cout << dist << endl;
                        if (dist < PARTICLE_DIAM) {
                            // cout << "actually found a collision" << endl;
                             collisionCount++;
                            // PUSH THEM BACK
                            long double r [2] = {((PARTICLE_DIAM-dist)/2)*(1/dist)*(centralParticle[0]-neighborParticle[0]), ((PARTICLE_DIAM-dist)/2)*(1/dist)*(centralParticle[1]-neighborParticle[1])};

                            particles[centralInd][0] += r[0];
                            particles[centralInd][1] += r[1];
                            particles[neighborInd][0] -= r[0];
                            particles[neighborInd][1] -= r[1];
                        }
                    }
                }
            }
            // cout << "- - - -" << endl;
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

    // cout << collisionCount << endl;
    return particles;
}

void exportNeighbors(array<vector<int>,TABLE_SIZE> neighbors, const char* saveName) {
    // open desired file for saving
    ofstream myfile (saveName);

    if (myfile.is_open()) {
        // add the array data to the file in the proper format
        for (int i = 0; i < TABLE_SIZE; ++i) {
            for (int j = 0; j < neighbors[i].size(); ++j) {
                myfile << to_string(neighbors[i][j]+1);
                myfile << " "; 
            }
            for (int k = 0; k < 19 - neighbors[i].size(); ++k) {
                myfile << "0 ";
            }
            myfile << std::fixed << "0";
            myfile << "\n";
        }
        myfile.close();
  }
    // // clear particle positions from heap
    // for (int i = 0; i < TABLE_SIZE; i++) {
    //     delete[] neighbors[i];
    // }
    // delete[] neighbors;
}

void exportMapping(int** neighborMapping, const char* saveName) {
    // open desired file for saving
    ofstream myfile (saveName);

    if (myfile.is_open()) {
        // add the array data to the file in the proper format
        for (int i = 0; i < TABLE_SIZE-1; ++i) {
            for (int j = 0; j < 7; ++j) {
                myfile << to_string(neighborMapping[i][j]+1);
                myfile << " "; 
            }
            myfile << std::fixed << neighborMapping[i][8]+1;
            myfile << "\n";
        }
        myfile.close();
  }
    // clear neighbor mapping from heap
    for (int i = 0; i < TABLE_SIZE; i++) {
        delete[] neighborMapping[i];
    }
    delete[] neighborMapping;
}

/*
    Runs the molecular dynamics simulation
    INPUTS:
            particles: 2d array containing particle positions
            numFrames: how long you want the simulation to run for
    OUTPUTS:
            particles: the final particle positions
*/
long double ** runSimCell(long double ** particles, int numFrames) {
    int frameCount = 0;

    int ** neighborMap = buildNeighborMapping();
    array<vector<int>,TABLE_SIZE> neighbors;
    while (frameCount < numFrames) {
        // get neighbor lists 
        neighbors = getNeighborsCell(particles, neighborMap);
        std::cout << "neighbors done" << endl;
        // array<vector<int>,NUM_PARTICLES> neighbors = getNeighborsSimple(particles);
        
        for (int i = 0; i < 10; ++i) {
           
            
            // move all the particles via brownian motion
            for (int m = 0; m < NUM_PARTICLES; ++m) {
                for (int n = 0; n < NUM_COLUMNS; ++n) {
                    particles[m][n] += d(gen)*SCALING_FACTOR;
                }
            }
            // do collision resolution
            particles = resolveCollisionsCell(particles, neighbors, neighborMap);
        }
        std::cout << "movement and collisions done" << endl;

        // report the frame as done
        frameCount++;
        std::cout << "frame " + to_string(frameCount) + " done" << endl;

        if (frameCount == numFrames) {
            exportNeighbors(neighbors, "neighbors.txt");
            exportMapping(neighborMap, "mapping.txt");
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
long double ** runSimSimple(long double ** particles, int numFrames) {
    int frameCount = 0;

    while (frameCount < numFrames) {
        // get neighbor lists 
        
        array<vector<int>,NUM_PARTICLES> neighbors = getNeighborsSimple(particles);
        std::cout << "neighbors done" << endl;
        
        for (int i = 0; i < 10; ++i) {
            // move all the particles via brownian motion
            for (int m = 0; m < NUM_PARTICLES; ++m) {
                for (int n = 0; n < NUM_COLUMNS; ++n) {
                    particles[m][n] += d(gen)*SCALING_FACTOR;
                }
            }

            // do collision resolution
            particles = resolveCollisionsSimple(particles, neighbors);
        }
        std::cout << "movement and collisions done" << endl;

        // report the frame as done
        frameCount++;
        std::cout << "frame " + to_string(frameCount) + " done" << endl;
    }

    return particles;
}


int main() {
    // import particle positions 
    long double** particles = importData("initial_particles.txt");

    cout << CELL_WIDTH << endl;
    cout << CELL_HEIGHT << endl;
    // std::cout << -1 % 95 << endl;
    // run the sim
    particles = runSimCell(particles, 1000);


    exportData(particles, "final_particles.txt");

    return 0;
}