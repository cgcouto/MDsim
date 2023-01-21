#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <array>
#include <random>

#include "MDsim-functions.h"

// - - - - - SIM PARAMETERS - - - - -

const int PARTICLE_DIAM = 10; // diameter of each simulated particles
const int NUM_COLUMNS = 2; // 
const int NUM_PARTICLES = 23892; // the number of particles we are simulating
const long double WIDTH = 1918.6; // exact width of the sim window
const long double HEIGHT = 1211.742744975187; // exact height of the sim window

const long double STEP_SIZE = 0.0050; // how much time elapsed in each frame (in seconds)
const long double VISCOSITY = 1.99e-3; // dynamic viscosity of our simulated colloid mixture
const int TEMPERATURE = 300; // temperature of the simulation (in kelvin)
const long double BOLTZ_CONSTANT = 1.381e-23; // boltzmann constant (in m^2*kg / s^2*K)
const long double SCALING_FACTOR = PARTICLE_DIAM/(1.27e-6); // size difference between simulated and real particles
const long double  SQUIG = (6*M_PI*VISCOSITY*(PARTICLE_DIAM/2))/SCALING_FACTOR; // Fd/v, derived from Stokes' law
const long double NEIGHB_THRESHOLD = 2; // shapes the size of the neighbor cells

const int MU = 0; // a variable in our normal distribution
const long double SIGMA = sqrt((2*BOLTZ_CONSTANT*TEMPERATURE)/SQUIG * STEP_SIZE);
random_device rd; // creates random integer seed values each time we run
unsigned int seed = rd(); // puts the seed into an int so we can store it if we wish
mt19937 gen(seed); // generator(seed) - mersenne twister based around 2^(19937)-1 
normal_distribution <float> d(MU,SIGMA); // plugs random numbers into normal distribution to move particles

// sizes the neighbor cells variably (you'll get slightly different result depending on the sim width/height)
const int TABLE_WIDTH = static_cast<int> (floor(WIDTH/(PARTICLE_DIAM*NEIGHB_THRESHOLD))); // no. of cells horizontally
const int TABLE_HEIGHT = static_cast<int> (floor(HEIGHT/(PARTICLE_DIAM*NEIGHB_THRESHOLD))); // no. of cells vertically

 // just TABLE_WIDTH * TABLE_HEIGHT, but must be valued at compile time to size arrays of vectors
const int TABLE_SIZE = 95*60;

// - - - - - SIM FUNCTIONS - - - - -

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
    Converts the given C++ particle array into a txt file with the given name
    INPUTS:
            particles: 2d array containing particle positions
            frameNum: the frame number of the data we're saving
            saveName: points to the name of the file you want to save to. must include file extension
*/
void exportData(long double** particles, int frameNum, const char* saveName) {
    // open desired file for appending data
    ofstream myfile;
    myfile.open(saveName, ios_base::app);

    if (myfile.is_open()) {
        // add the array data to the file in the proper format
        for (int i = 0; i < NUM_PARTICLES; ++i) {
            for (int j = 0; j < NUM_COLUMNS; ++j) {
                myfile << to_string(particles[i][j]); // positions in first and second cols
                myfile << " "; 
            }
            myfile << std::fixed << frameNum; // frame number in third column - needed for pos_to_cnt_cell
            myfile << "\n";
        }
        myfile.close();
  }
}

/*
    Saves the given neighbor array as a txt file with dimensions TABLE_SIZE x 20
    Extra zeros are added to each row to get it to 20 total entries
    INPUTS:
            neighbors: array of vectors that contains the indices of the particles within each cell
            saveName: points to the name of the file you want to save to. must include file extension
*/
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
}

/*
    Saves the given cell mapping as a txt file with dimensions TABLE_SIZE x 8
    INPUTS:
            neighbors: 2d array of ints that contains the indices of the cells that each cell is neighbors with
            saveName: points to the name of the file you want to save to. must include file extension
*/
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
    Creates and stores the neighbors of each cell so we can access their indices in constant time
    Imagine a checkerboard where each square has an id - for each square we want to log the id's of its eight
    surrounding squares
    OUTPUTS:
            neighborMapping: a 2d array of ints that stores the neighbors of each cell
*/
int ** buildNeighborMapping() {
    // create the 2d array of ints
    int** neighborMapping = new int*[TABLE_SIZE];
    for (int i = 0; i < TABLE_SIZE; ++i) {
        neighborMapping[i] = new int[8];
    }

    // loop through all the possible cell coordinates
    for (int m = 0; m < TABLE_WIDTH; ++m) {
        for (int n = 0; n < TABLE_HEIGHT; ++n) {

            int tableIndex = m + TABLE_WIDTH * n; // index of current cell

            // do our mapping and encode it in indices
            // the ninth neighbor (itself) is implied by the tableIndex
            neighborMapping[tableIndex][0] = ((TABLE_WIDTH + ((m-1))) % TABLE_WIDTH) + TABLE_WIDTH*n;
            neighborMapping[tableIndex][1] = (((TABLE_WIDTH + ((m+1))) % TABLE_WIDTH) % TABLE_WIDTH) + TABLE_WIDTH*n;
            neighborMapping[tableIndex][2] = m + TABLE_WIDTH*(((TABLE_HEIGHT + ((n-1))) % TABLE_HEIGHT) % TABLE_HEIGHT);
            neighborMapping[tableIndex][3] = m + TABLE_WIDTH*(((TABLE_HEIGHT + ((n+1))) % TABLE_HEIGHT) % TABLE_HEIGHT);
            neighborMapping[tableIndex][4] = (((TABLE_WIDTH + ((m-1))) % TABLE_WIDTH) % TABLE_WIDTH) 
                                            + TABLE_WIDTH*(((TABLE_HEIGHT + ((n-1))) % TABLE_HEIGHT) % TABLE_HEIGHT);
            neighborMapping[tableIndex][5] = (((TABLE_WIDTH + ((m-1))) % TABLE_WIDTH) % TABLE_WIDTH) 
                                            + TABLE_WIDTH*(((TABLE_HEIGHT + ((n+1))) % TABLE_HEIGHT) % TABLE_HEIGHT);
            neighborMapping[tableIndex][6] = (((TABLE_WIDTH + ((m+1))) % TABLE_WIDTH) % TABLE_WIDTH) 
                                            + TABLE_WIDTH*(((TABLE_HEIGHT + ((n-1))) % TABLE_HEIGHT) % TABLE_HEIGHT);
            neighborMapping[tableIndex][7] = (((TABLE_WIDTH + ((m+1))) % TABLE_WIDTH) % TABLE_WIDTH) 
                                            + TABLE_WIDTH*(((TABLE_HEIGHT + ((n+1))) % TABLE_HEIGHT) % TABLE_HEIGHT);
        }
    }

    return neighborMapping;
}

/*
    Creates a neighbor list for the particles using the cell list method
    Find and store which cell each particle is in
    INPUTS:
            particles: 2d array containing particle positions
            neighborMapping: 2d array storing the neighbors of each cell
    OUTPUTS:
            neighbors: an array of vectors containing relevant array indexes from particles that fit in each cell
*/
array<vector<int>,TABLE_SIZE> getNeighbors(long double** particles, int** neighborMapping) {
    array<vector<int>,TABLE_SIZE> neighbors;

    for (int i = 0; i < NUM_PARTICLES; ++i) {
        // get the position of the cell the particle is in
        int tableX = static_cast<int> (floor(particles[i][0] / (WIDTH/TABLE_WIDTH)));
        int tableY = static_cast<int> (floor(particles[i][1] / (HEIGHT/TABLE_HEIGHT)));
        
        // translate position to index and add to vector
        neighbors[tableX + TABLE_WIDTH*tableY].push_back(i);

    }

    return neighbors;
}

/*
    Uses the neighbor list to resolve particle collisions
    INPUTS:
            particles: 2d array containing particle positions
            neighbors: array of vectors containing particle indices that lie in each cell
    OUTPUTS:
            particles: same as input, just with collisions resolved
*/
long double ** resolveCollisions(long double** particles, array<vector<int>,TABLE_SIZE> neighbors, int** neighborMapping) {
    
    // go through all the indices in the neighbor list (size is equal to NUM_PARTICLES)
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

                    // we don't want a particle to consider itself as its own neighbor
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
                        long double dist = sqrt(pow(centralParticle[0]-neighborParticle[0],2) 
                                               + pow(centralParticle[1]-neighborParticle[1],2));

                        // if the particles are overlapping, repel them so they no longer intersect
                        if (dist < PARTICLE_DIAM) {

                            long double r [2] = {((PARTICLE_DIAM-dist)/2)*(1/dist)*(centralParticle[0]-neighborParticle[0]), 
                                                 ((PARTICLE_DIAM-dist)/2)*(1/dist)*(centralParticle[1]-neighborParticle[1])};

                            particles[centralInd][0] += r[0];
                            particles[centralInd][1] += r[1];
                            particles[neighborInd][0] -= r[0];
                            particles[neighborInd][1] -= r[1];
                        }
                    }
                }
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
    Runs the molecular dynamics simulation with cell lists
    INPUTS:
            particles: 2d array containing particle positions
            numFrames: how long you want the simulation to run for
            savingFrequency: how often you want to save particle positions
            saveName: the name of the file you want to save into
    OUTPUTS:
            particles: the final particle positions
*/
long double ** runSim(long double ** particles, int numFrames, int savingFrequency, const char * saveName) {

    // clear the given txt file of any existing data
    ofstream ofs;
    ofs.open(saveName, std::ofstream::out | std::ofstream::trunc);
    ofs.close();

    int frameCount = 0;

    // store our mapping of cells to their neighbors
    int ** neighborMap = buildNeighborMapping();
    array<vector<int>,TABLE_SIZE> neighbors;

    while (frameCount < numFrames) {
        // get neighbor lists 
        neighbors = getNeighbors(particles, neighborMap);
        
        for (int i = 0; i < 10; ++i) {
           
            // move all the particles via brownian motion
            for (int m = 0; m < NUM_PARTICLES; ++m) {
                for (int n = 0; n < NUM_COLUMNS; ++n) {
                    particles[m][n] += d(gen)*SCALING_FACTOR;
                }
            }

            // do collision resolution
            particles = resolveCollisions(particles, neighbors, neighborMap);
        }

        // report the frame as done
        frameCount++;


        // save the frame data if it aligns with the saving frequency
        if (frameCount % savingFrequency == 0) {
            exportData(particles, frameCount, saveName);
            cout << "saved frame " + to_string(frameCount) << endl;
        } else {
            cout << "frame " + to_string(frameCount) + " done" << endl;
        }

        

        // save the final neighbor list and cell mapping
        if (frameCount == numFrames) {
            exportNeighbors(neighbors, "neighbors.txt");
            exportMapping(neighborMap, "mapping.txt");
        }
    }


    return particles;
}

/*
    Runs one frame of the molecular dynamics simulation
    Ideal for real-time visualization in a render loop
    INPUTS:
            particles: 2d array containing particle positions
    OUTPUTS:
            particles: the updated particle positions
*/
long double ** doOneFrame(long double ** particles) {

    // store our mapping of cells to their neighbors
    int ** neighborMap = buildNeighborMapping();
    array<vector<int>,TABLE_SIZE> neighbors;


    // get neighbor lists 
     neighbors = getNeighbors(particles, neighborMap);
        
    for (int i = 0; i < 10; ++i) {
           
        // move all the particles via brownian motion
        for (int m = 0; m < NUM_PARTICLES; ++m) {
            for (int n = 0; n < NUM_COLUMNS; ++n) {
                particles[m][n] += d(gen)*SCALING_FACTOR;
            }
        }

        // do collision resolution
        particles = resolveCollisions(particles, neighbors, neighborMap);
    }

    return particles;
}