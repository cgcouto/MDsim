#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace std; // so we don't have to put std:: in front of certain things


// - - - - - IMPORTANT CONSTANTS - - - - -

const int PARTICLE_SIZE = 10;
const int NUM_COLUMNS = 3;
const int NUM_PARTICLES = 23892;
const long double WIDTH = 1920;
const long double HEIGHT = 1215;

const long double VISCOSITY = 1.99e-3;
const int TEMPERATURE = 300;
const long double BOLTZ_CONSTANT = 1.381e-23;
const long double SCALING_FACTOR = PARTICLE_SIZE/(1.3e-6);
const long double  SQUIG = (6*M_PI*(PARTICLE_SIZE/2))/SCALING_FACTOR;


/*
    Translates a MATLAB 2D particle array (that is saved via writematrix)
    into a C++ array
*/
long double** importData(const char* fileName)
{
    ifstream file(fileName);

    // make space for the particle positions on the heap
    long double** initialParticles = new long double*[NUM_PARTICLES];
    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        initialParticles[i] = new long double[NUM_COLUMNS];
    }

    // string* initialParticles = new string[23892][3];
    if (file.is_open()) 
    {
        string data = "";
        int count = 0;
        while(getline(file, data)){   // get a whole line
            stringstream ss(data);
            while(getline(ss, data, ',')){ // get each element in the line
                initialParticles[count/NUM_COLUMNS][count % NUM_COLUMNS] = stold(data);
                cout << to_string(initialParticles[count/NUM_COLUMNS][count % NUM_COLUMNS]) << endl;
                count++;
            }
        }
    }
    file.close();

    return initialParticles;
}

void exportData(long double** particles, const char* saveName)
{
    ofstream myfile (saveName);
    if (myfile.is_open())
    {
        for (int i = 0; i < NUM_PARTICLES; ++i) 
        {
            for (int j = 0; j < NUM_COLUMNS-1; ++j)
            {
                myfile << to_string(particles[i][j]);
                myfile << ","; 
            }
            myfile << std::fixed << particles[i][NUM_COLUMNS-1];
            myfile << "\n";
        }
        myfile.close();
  }
    // clear particle positions from heap
    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        delete[] particles[i];
    }
    delete[] particles;
}

int getNeighborsSimple(long double** particles)
{
    // loop through each row in neighbors and check distances
    return 0;
}
// make sure you pass references into these functions!
int resolveCollisions(long double** particles, double** neighbors)
{
    return 0;
}

// exportData

int runSim(long double particles, int numFrames)
{
    int frameCount = 0;

    while (frameCount < numFrames)
    {
        // neighbors = getNeighborsSimple(particles);
        for (int i = 0; i < 10; ++i)
        {
            // do brownian motion and collision checks
            // use random numbers on a normal distribution
            // particles = resolveCollisions(particles, neighbors);
        }
    }

    return 0;
}


int main()
{
    long double** particles = importData("initial_particles.txt");

    exportData(particles, "test.txt");

    // cout << "frame " + to_string(5) + " done" << endl;

    return 0;
}