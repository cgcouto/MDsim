#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std; // so we don't have to put std:: in front of certain things

const int PARTICLE_SIZE = 10;
const int NUM_PARTICLES = 23892;

/*
    Translates a MATLAB 2D particle array (that is saved via writematrix)
    into a C++ array
*/
double** importData(const char* fileName, const int numColumns)
{
    const int c = 3;
    ifstream file(fileName);

    double** initialParticles = new double*[NUM_PARTICLES];
    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        initialParticles[i] = new double[numColumns];
    }

    // string* initialParticles = new string[23892][3];
    if (file.is_open()) 
    {
        string data = "";
        int count = 0;
        while(getline(file, data)){   // get a whole line
            stringstream ss(data);
            while(getline(ss, data, ',')){ // get each element in the line
                initialParticles[count/numColumns][count % numColumns] = stod(data);
                cout << to_string(initialParticles[count/numColumns][count % numColumns]) << endl;
                count++;
            }
        }
    }
    file.close();

    return initialParticles;
}

int exportData(double** particles, const char* saveName)
{
    return 0;
}

int getNeighborsSimple(double** particles)
{
    // loop through each row in neighbors and check distances
    return 0;
}
// make sure you pass references into these functions!
int resolveCollisions(double** particles, double** neighbors)
{
    return 0;
}

// exportData

int runSim(double particles, int numFrames)
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
    double** particles = importData("initial_particles.txt", 3);

    // clear particle positions from heap
    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        delete[] particles[i];
    }
    delete[] particles;

    return 0;
}