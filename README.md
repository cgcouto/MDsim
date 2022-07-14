# MDsim

## Introduction

Welcome to MDsim! The goal of this project was to translate the original MDsim code from MATLAB to C++ to improve sim performance. It was written during the summer of 2022. These simulations function similarly to MDsim in MATLAB, but the implementations differ in a few choice areas.

## Overview

Let's start by walking through how a single frame of sim is formed. The code starts by generating neighbor lists. A neighbor list is a 2D array where each row correspond to the index of a particle in the simulation and the row entries are the particles in the general vicinity of the row particle. These neighbor lists are vital because they make resolving collisions much more efficient. Instead of checking your particle against every other particle in the sim, you only check it against the particles in its neighbor list.

Once we have neighbor lists, we will do ten steps of Brownian motion and collision resolution. Doing the Brownian motion involves adding some random displacement to each of the particle positions, with the displacement drawn from a normal distribution set by sim parameters. Resolving collisions involves going through the neighbor lists, retrieving particle positions, and checking distances. If the distance is less than a particle diameter, the particles must be overlapping, and we give them a slight nudge away from each other to resolve this.

With out ten steps completed, our frame of sim is done! The last bit of work depends on the saving frequency. If we're due to save our particle positions (determined by whether the current frame is a multiple of the saving frequency) we will append to our plist file and add our current particle positions. If we're not due to save, we just progress to the next frame of sim!

## Using This Code

To start, you're going to want to run the 'importIntoC.m' script over in MATLAB. It will create MDsim objects in labelled folders and place all the initial particle positions into txt files at the top-level. 

From here, copy the lines of C++ code produced by the MATLAB script and open MDsim.cpp with your IDE of choice. Once you're there, replace the lines at the top of the code with what you've copied. This is all you need to change in the code! Now you should compile the C++ code in the terminal with the command g++ -o  mdsim MDsim.cpp -fopenmp . From there, it's as simple as running the executable with ./mdsim.

Finally, to read in the plist txt files and create cnt cells in the labelled folders. Make sure that the sim info at the top of the script matches that from 'importIntoC.m'. 

## Differences from the MATLAB Version

### Building Neighbor Lists

As discussed in the overview, building neighbor lists is a crucial part of these molecular dynamics simulations. In MATLAB we loop through each of the particles and use logical indexing to get all the particles neighbor threshold. We do not have logical indexing in C++, so building neighbor lists efficiently demands a different approach. 

As a result, we've implemented cell lists. The Wikipedia page does a good job explaining how it works, but the gist is that we partition the 2D sim space into cells. Before we run the sim, we need to Having this stored in an array allows us to access this information quickly and not waste time calculating it again. 

### Data Structures

MATLAB makes it really simple to generate and work with 2D arrays. In C++ it is much Once they are created, working with 2D arrays in C++ is similar to MATLAB. Just remember that your indices start at 0, not 1! The neighbor information is stored in a vector of vectors. Ideally we would use an array of vectors - . We only use a bit more memory with this vector of vectors approach, the performance is the same.

### Generating Random Numbers
For these simulations, the Brownian motion is dictated by random numbers drawn on a normal distribution. In MATLAB, we can do this with the function normrnd. In C++, however, getting random numbers takes a bit more work.

### Running Multiple Trials

In MATLAB, parallelizing our loops is really easy. We just use parfor and tweak the parallel pool settings as needed.  In C++, the solutions give us much more control. We could have used pthreads - the built-in solution for multi-threading - but we decided to use OpenMP because it's more flexible and more intuitive to implement. If you want to know more about OpenMP, there's a book in the lab called 'The OpenMP Common Core'. It's a good read - the first two chapters explain the fundamentals of parallel computing really well.

## Closing Thoughts

So that's the code base! I hope that it is useful to you. If you have any questions, feel free to shoot me an email at ccouto@g.hmc.edu. 

Chris Couto (HMC '23)
