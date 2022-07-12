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
#include <chrono>


// things i will need to fix in the future

// the sig fig problem with exportData

using namespace std; // so we don't have to put std:: in front of certain things


// - - - - - SIM PARAMETERS - - - - - REPLACE ME WHEN YOU ARE RUNNING SIM

const string INITIAL_FILE [48] = {"initial_particles_theta5tightness1.06extra_forces0flr_01.txt","initial_particles_theta5tightness1.06extra_forces0flr_02.txt","initial_particles_theta5tightness1.06extra_forces0flr_03.txt","initial_particles_theta5tightness1.06extra_forces0flr_04.txt","initial_particles_theta5tightness1.06extra_forces0flr_05.txt","initial_particles_theta5tightness1.06extra_forces0flr_06.txt","initial_particles_theta5tightness1.06extra_forces0flr_07.txt","initial_particles_theta5tightness1.06extra_forces0flr_08.txt","initial_particles_theta5tightness1.06extra_forces0flr_09.txt","initial_particles_theta5tightness1.06extra_forces0flr_10.txt","initial_particles_theta5tightness1.06extra_forces0flr_11.txt","initial_particles_theta5tightness1.06extra_forces0flr_12.txt","initial_particles_theta10tightness1.06extra_forces0flr_01.txt","initial_particles_theta10tightness1.06extra_forces0flr_02.txt","initial_particles_theta10tightness1.06extra_forces0flr_03.txt","initial_particles_theta10tightness1.06extra_forces0flr_04.txt","initial_particles_theta10tightness1.06extra_forces0flr_05.txt","initial_particles_theta10tightness1.06extra_forces0flr_06.txt","initial_particles_theta10tightness1.06extra_forces0flr_07.txt","initial_particles_theta10tightness1.06extra_forces0flr_08.txt","initial_particles_theta10tightness1.06extra_forces0flr_09.txt","initial_particles_theta10tightness1.06extra_forces0flr_10.txt","initial_particles_theta10tightness1.06extra_forces0flr_11.txt","initial_particles_theta10tightness1.06extra_forces0flr_12.txt","initial_particles_theta15tightness1.06extra_forces0flr_01.txt","initial_particles_theta15tightness1.06extra_forces0flr_02.txt","initial_particles_theta15tightness1.06extra_forces0flr_03.txt","initial_particles_theta15tightness1.06extra_forces0flr_04.txt","initial_particles_theta15tightness1.06extra_forces0flr_05.txt","initial_particles_theta15tightness1.06extra_forces0flr_06.txt","initial_particles_theta15tightness1.06extra_forces0flr_07.txt","initial_particles_theta15tightness1.06extra_forces0flr_08.txt","initial_particles_theta15tightness1.06extra_forces0flr_09.txt","initial_particles_theta15tightness1.06extra_forces0flr_10.txt","initial_particles_theta15tightness1.06extra_forces0flr_11.txt","initial_particles_theta15tightness1.06extra_forces0flr_12.txt","initial_particles_theta20tightness1.06extra_forces0flr_01.txt","initial_particles_theta20tightness1.06extra_forces0flr_02.txt","initial_particles_theta20tightness1.06extra_forces0flr_03.txt","initial_particles_theta20tightness1.06extra_forces0flr_04.txt","initial_particles_theta20tightness1.06extra_forces0flr_05.txt","initial_particles_theta20tightness1.06extra_forces0flr_06.txt","initial_particles_theta20tightness1.06extra_forces0flr_07.txt","initial_particles_theta20tightness1.06extra_forces0flr_08.txt","initial_particles_theta20tightness1.06extra_forces0flr_09.txt","initial_particles_theta20tightness1.06extra_forces0flr_10.txt","initial_particles_theta20tightness1.06extra_forces0flr_11.txt","initial_particles_theta20tightness1.06extra_forces0flr_12.txt"};
const string FINAL_FILE [48] = {"plist_theta5tightness1.06extra_forces0flr_01.txt","plist_theta5tightness1.06extra_forces0flr_02.txt","plist_theta5tightness1.06extra_forces0flr_03.txt","plist_theta5tightness1.06extra_forces0flr_04.txt","plist_theta5tightness1.06extra_forces0flr_05.txt","plist_theta5tightness1.06extra_forces0flr_06.txt","plist_theta5tightness1.06extra_forces0flr_07.txt","plist_theta5tightness1.06extra_forces0flr_08.txt","plist_theta5tightness1.06extra_forces0flr_09.txt","plist_theta5tightness1.06extra_forces0flr_10.txt","plist_theta5tightness1.06extra_forces0flr_11.txt","plist_theta5tightness1.06extra_forces0flr_12.txt","plist_theta10tightness1.06extra_forces0flr_01.txt","plist_theta10tightness1.06extra_forces0flr_02.txt","plist_theta10tightness1.06extra_forces0flr_03.txt","plist_theta10tightness1.06extra_forces0flr_04.txt","plist_theta10tightness1.06extra_forces0flr_05.txt","plist_theta10tightness1.06extra_forces0flr_06.txt","plist_theta10tightness1.06extra_forces0flr_07.txt","plist_theta10tightness1.06extra_forces0flr_08.txt","plist_theta10tightness1.06extra_forces0flr_09.txt","plist_theta10tightness1.06extra_forces0flr_10.txt","plist_theta10tightness1.06extra_forces0flr_11.txt","plist_theta10tightness1.06extra_forces0flr_12.txt","plist_theta15tightness1.06extra_forces0flr_01.txt","plist_theta15tightness1.06extra_forces0flr_02.txt","plist_theta15tightness1.06extra_forces0flr_03.txt","plist_theta15tightness1.06extra_forces0flr_04.txt","plist_theta15tightness1.06extra_forces0flr_05.txt","plist_theta15tightness1.06extra_forces0flr_06.txt","plist_theta15tightness1.06extra_forces0flr_07.txt","plist_theta15tightness1.06extra_forces0flr_08.txt","plist_theta15tightness1.06extra_forces0flr_09.txt","plist_theta15tightness1.06extra_forces0flr_10.txt","plist_theta15tightness1.06extra_forces0flr_11.txt","plist_theta15tightness1.06extra_forces0flr_12.txt","plist_theta20tightness1.06extra_forces0flr_01.txt","plist_theta20tightness1.06extra_forces0flr_02.txt","plist_theta20tightness1.06extra_forces0flr_03.txt","plist_theta20tightness1.06extra_forces0flr_04.txt","plist_theta20tightness1.06extra_forces0flr_05.txt","plist_theta20tightness1.06extra_forces0flr_06.txt","plist_theta20tightness1.06extra_forces0flr_07.txt","plist_theta20tightness1.06extra_forces0flr_08.txt","plist_theta20tightness1.06extra_forces0flr_09.txt","plist_theta20tightness1.06extra_forces0flr_10.txt","plist_theta20tightness1.06extra_forces0flr_11.txt","plist_theta20tightness1.06extra_forces0flr_12.txt"};
const int PARTICLE_DIAM [48] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
const int NUM_PARTICLES [48] = {23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892, 23892};
const long double SIM_WIDTH [48] = {1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6, 1918.6};
const long double SIM_HEIGHT [48] = {1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427, 1211.7427};
const int NEIGHB_THRESHOLD [48] = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
const int TABLE_WIDTH [48] = {95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95};
const int TABLE_HEIGHT [48] = {60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60};
const int TABLE_SIZE [48] = {5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700, 5700};
const int NUM_FRAMES [48] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const int SAVING_FREQUENCY [48] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};


// - - - - - IMPORTANT CONSTANTS - - - - -


/*
    Translates a MATLAB 2D particle array (that is saved via writematrix)
    into a C++ array
    INPUTS:
            fileName: points to the name of your desired file. must include file extension
    OUTPUTS:
            initialParticles: the particle positions in a 2d array format 
*/
long double** importData(const char* fileName, int trialNum) {
    ifstream file(fileName);

    // make space for the particle positions on the heap
    long double** initialParticles = new long double*[NUM_PARTICLES[trialNum]];
    for (int i = 0; i < NUM_PARTICLES[trialNum]; i++) {
        initialParticles[i] = new long double[2];
    }

    // get all the data into our array
    if (file.is_open()) {
        string data = "";
        int count = 0;
        while(getline(file, data)) {   // get a whole line
            stringstream ss(data);
            while(getline(ss, data, ' ')) { // get each element in the line
                initialParticles[count/2][count % 2] = stold(data);
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
void exportData(long double** particles, int frameNum, const char* saveName, int trialNum) {
    // open desired file for appending data
    ofstream myfile;
    myfile.open(saveName, ios_base::app);

    if (myfile.is_open()) {
        // add the array data to the file in the proper format
        for (int i = 0; i < NUM_PARTICLES[trialNum]; ++i) {
            for (int j = 0; j < 2; ++j) {
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
void exportNeighbors(vector<vector<int>> neighbors, const char* saveName, int trialNum) {
    // open desired file for saving
    ofstream myfile (saveName);

    if (myfile.is_open()) {
        // add the array data to the file in the proper format
        for (int i = 0; i < TABLE_SIZE[trialNum]; ++i) {
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
void exportMapping(int** neighborMapping, const char* saveName, int trialNum) {
    // open desired file for saving
    ofstream myfile (saveName);

    if (myfile.is_open()) {
        // add the array data to the file in the proper format
        for (int i = 0; i < TABLE_SIZE[trialNum]; ++i) {
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
    for (int i = 0; i < TABLE_SIZE[trialNum]; i++) {
        delete[] neighborMapping[i];
    }
    delete[] neighborMapping;
}

/*
    Creates and stores the neighbors of each cell so we can access their indices in constant time
    OUTPUTS:
            neighborMapping: a 2d array of ints that stores the neighbors of each cell
*/
int ** buildNeighborMapping(int trialNum) {
    // create the 2d array of ints
    int** neighborMapping = new int*[TABLE_SIZE[trialNum]];
    for (int i = 0; i < TABLE_SIZE[trialNum]; ++i) {
        neighborMapping[i] = new int[8];
    }

    // loop through all the possible cell coordinates
    for (int m = 0; m < TABLE_WIDTH[trialNum]; ++m) {
        for (int n = 0; n < TABLE_HEIGHT[trialNum]; ++n) {

            int tableIndex = m + TABLE_WIDTH[trialNum] * n;
            // do our mapping and encode it in indices
            // the ninth neighbor (itself) is implied by the tableIndex
            neighborMapping[tableIndex][0] = ((TABLE_WIDTH[trialNum] + ((m-1))) % TABLE_WIDTH[trialNum]) + TABLE_WIDTH[trialNum]*n;
            neighborMapping[tableIndex][1] = (((TABLE_WIDTH[trialNum] + ((m+1))) % TABLE_WIDTH[trialNum]) % TABLE_WIDTH[trialNum]) + TABLE_WIDTH[trialNum]*n;
            neighborMapping[tableIndex][2] = m + TABLE_WIDTH[trialNum]*(((TABLE_HEIGHT[trialNum] + ((n-1))) % TABLE_HEIGHT[trialNum]) % TABLE_HEIGHT[trialNum]);
            neighborMapping[tableIndex][3] = m + TABLE_WIDTH[trialNum]*(((TABLE_HEIGHT[trialNum] + ((n+1))) % TABLE_HEIGHT[trialNum]) % TABLE_HEIGHT[trialNum]);
            neighborMapping[tableIndex][4] = (((TABLE_WIDTH[trialNum] + ((m-1))) % TABLE_WIDTH[trialNum]) % TABLE_WIDTH[trialNum]) + TABLE_WIDTH[trialNum]*(((TABLE_HEIGHT[trialNum] + ((n-1))) % TABLE_HEIGHT[trialNum]) % TABLE_HEIGHT[trialNum]);
            neighborMapping[tableIndex][5] = (((TABLE_WIDTH[trialNum] + ((m-1))) % TABLE_WIDTH[trialNum]) % TABLE_WIDTH[trialNum]) + TABLE_WIDTH[trialNum]*(((TABLE_HEIGHT[trialNum] + ((n+1))) % TABLE_HEIGHT[trialNum]) % TABLE_HEIGHT[trialNum]);
            neighborMapping[tableIndex][6] = (((TABLE_WIDTH[trialNum] + ((m+1))) % TABLE_WIDTH[trialNum]) % TABLE_WIDTH[trialNum]) + TABLE_WIDTH[trialNum]*(((TABLE_HEIGHT[trialNum] + ((n-1))) % TABLE_HEIGHT[trialNum]) % TABLE_HEIGHT[trialNum]);
            neighborMapping[tableIndex][7] = (((TABLE_WIDTH[trialNum] + ((m+1))) % TABLE_WIDTH[trialNum]) % TABLE_WIDTH[trialNum]) + TABLE_WIDTH[trialNum]*(((TABLE_HEIGHT[trialNum] + ((n+1))) % TABLE_HEIGHT[trialNum]) % TABLE_HEIGHT[trialNum]);
        }
    }

    return neighborMapping;
}

/*
    Creates a neighbor list for the particles using the cell list method
    INPUTS:
            particles: 2d array containing particle positions
            neighborMapping: 2d array storing the neighbors of each cell
    OUTPUTS:
            neighbors: an array of vectors containing relevant array indexes from particles that fit in each cell
*/
vector<vector<int>> getNeighbors(long double** particles, int** neighborMapping, int trialNum) {
    vector<vector<int>> neighbors(TABLE_SIZE[trialNum]);

    for (int i = 0; i < NUM_PARTICLES[trialNum]; ++i) {
        // get the position of the cell the particle is in
        int tableX = static_cast<int> (floor(particles[i][0] / (SIM_WIDTH[trialNum]/TABLE_WIDTH[trialNum])));
        int tableY = static_cast<int> (floor(particles[i][1] / (SIM_HEIGHT[trialNum]/TABLE_HEIGHT[trialNum])));
        
        // translate position to index and add to vector
        neighbors[tableX + TABLE_WIDTH[trialNum]*tableY].push_back(i);

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
long double ** resolveCollisions(long double** particles, vector<vector<int>> neighbors, int** neighborMapping, int trialNum) {
    
    // go through all the indices in the neighbor list (size is equal to NUM_PARTICLES)
    for (int i = 0; i < TABLE_SIZE[trialNum]; ++i) {
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
                        if (neighborParticle[0] < PARTICLE_DIAM[trialNum] && centralParticle[0] > SIM_WIDTH[trialNum]-PARTICLE_DIAM[trialNum]) {
                            neighborParticle[0] += SIM_WIDTH[trialNum];
                        } else if (neighborParticle[0] > SIM_WIDTH[trialNum]-PARTICLE_DIAM[trialNum] && centralParticle[0] < PARTICLE_DIAM[trialNum]) {
                            neighborParticle[0] -= SIM_WIDTH[trialNum];
                        }

                        if (neighborParticle[1] < PARTICLE_DIAM[trialNum] && centralParticle[1] > SIM_HEIGHT[trialNum]-PARTICLE_DIAM[trialNum]) {
                            neighborParticle[1] += SIM_HEIGHT[trialNum];
                        } else if (neighborParticle[1] > SIM_HEIGHT[trialNum]-PARTICLE_DIAM[trialNum] && centralParticle[1] < PARTICLE_DIAM[trialNum]) {
                            neighborParticle[1] -= SIM_HEIGHT[trialNum];
                        }

                        // get distance between points
                        long double dist = sqrt(pow(centralParticle[0]-neighborParticle[0],2) + pow(centralParticle[1]-neighborParticle[1],2));

                        if (dist < PARTICLE_DIAM[trialNum]) {

                            // PUSH THEM BACK
                            long double r [2] = {((PARTICLE_DIAM[trialNum]-dist)/2)*(1/dist)*(centralParticle[0]-neighborParticle[0]), ((PARTICLE_DIAM[trialNum]-dist)/2)*(1/dist)*(centralParticle[1]-neighborParticle[1])};

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
    for (int i = 0; i < NUM_PARTICLES[trialNum]; ++i) {
        if (particles[i][0] < 0) {
            particles[i][0] += SIM_WIDTH[trialNum];
        } else if (particles[i][0] > SIM_WIDTH[trialNum]) {
            particles[i][0] -= SIM_WIDTH[trialNum];
        }

        if (particles[i][1] < 0) {
            particles[i][1] += SIM_HEIGHT[trialNum];
        } else if (particles[i][1] > SIM_HEIGHT[trialNum]) {
            particles[i][1] -= SIM_HEIGHT[trialNum];
        }
    }

   return particles;
}

/*
    Runs the molecular dynamics simulation with cell lists - O(N)
    INPUTS:
            particles: 2d array containing particle positions
            numFrames: how long you want the simulation to run for
            savingFrequency: how often you want to save particle positions
            saveName: the name of the file you want to save into
            trialNum: the number of the trial so we know where to pull sim details from
    OUTPUTS:
            particles: the final particle positions
*/
long double ** runSim(long double ** particles, int numFrames, int savingFrequency, const char * saveName, int trialNum) {

    // clear the given txt file of any existing data
    ofstream ofs;
    ofs.open(saveName, std::ofstream::out | std::ofstream::trunc);
    ofs.close();

    int frameCount = 0;

    const long double STEP_SIZE = 0.0050;
    const long double VISCOSITY = 1.99e-3;
    const int TEMPERATURE = 300;
    const long double BOLTZ_CONSTANT = 1.381e-23;
    const long double SCALING_FACTOR = PARTICLE_DIAM[trialNum]/(1.3e-6);
    const long double  SQUIG = (6*M_PI*VISCOSITY*(PARTICLE_DIAM[trialNum]/2))/SCALING_FACTOR;

    const int MU = 0;
    const long double SIGMA = sqrt((2*BOLTZ_CONSTANT*TEMPERATURE)/SQUIG * STEP_SIZE);
    random_device rd; // creates random integer seed values each time we run
    unsigned int seed = rd(); // put the seed into an int so we can store it if we wish
    mt19937 gen(seed); // generator(seed) - mersenne twister based around 2^(19937)-1 
    normal_distribution <float> d(MU,SIGMA);

    // exportData(particles, frameCount, saveName);

    // store our mapping of cells to their neighbors
    int ** neighborMap = buildNeighborMapping(trialNum);
    vector<vector<int>> neighbors;

    while (frameCount < numFrames) {
        // get neighbor lists 
        neighbors = getNeighbors(particles, neighborMap, trialNum);
        
        for (int i = 0; i < 10; ++i) {
           
            // move all the particles via brownian motion
            for (int m = 0; m < NUM_PARTICLES[trialNum]; ++m) {
                for (int n = 0; n < 2; ++n) {
                    particles[m][n] += d(gen)*SCALING_FACTOR;
                }
            }


            // do collision resolution
            particles = resolveCollisions(particles, neighbors, neighborMap, trialNum);
        }

        // report the frame as done
        frameCount++;

        

        // save the frame data if it aligns with the saving frequency
        if (frameCount % savingFrequency == 0) {
            exportData(particles, frameCount, saveName, trialNum);
            cout << "saved frame " + to_string(frameCount) << endl;
        } else {
            cout << "frame " + to_string(frameCount) + " done" << endl;
        }
    }


    return particles;
} 


int main() {
    #pragma omp parallel 
    {
        # pragma omp for
            for (int i = 0; i < 48; ++i) {
                // start the timer
                std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
                
                // import particle positions 
                long double** particles = importData(INITIAL_FILE[i].c_str(), i);

                // run the sim
                particles = runSim(particles, NUM_FRAMES[i], SAVING_FREQUENCY[i], FINAL_FILE[i].c_str(), i);

                // end the timer and print the time elapsed
                std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
            }

    }
    return 0;
}