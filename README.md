# MDsim

## Introduction

Welcome to MDsim! The goal of this project was to translate the original MDsim code from MATLAB to C++ to improve sim performance. It was written during the summer of 2022. These simulations function similarly to MDsim in MATLAB, but the implementations differ in a few choice areas.

## Overview

Let's start by walking through how a single frame of sim is formed. The code starts by generating neighbor lists. A neighbor list is a 2D array where each row correspond to the index of a particle in the simulation and the row entries are the particles in the general vicinity of the row particle. These neighbor lists are vital because they make resolving collisions much more efficient. Instead of checking your particle against every other particle in the sim, you only check it against the particles in its neighbor list.

Once we have neighbor lists, we will do ten steps of Brownian motion and collision resolution. Resolving collisions involves going through the neighbor lists, retrieving particle positions, and checking distances. If the distance is less than a particle diameter, the particles must be overlapping, and we give them a slight nudge away from each other to resolve this.

## Differences from the MATLAB Version

### Building Neighbor Lists

As discussed in the overview, building neighbor lists is a crucial part of these molecular dynamics simulations. In MATLAB we loop through each of the particles use logical indexing to . We do not have logical indexing in C++, so building neighbor lists efficiently demands a different approach. 

As a result, we've implemented cell lists. The Wikipedia page does a good job explaining how it works, but the gist is that we partition the 2D sim space into cells. Before we run the sim, we need to Having this stored in an array allows us to access this information quickly and not waste time calculating it again. 

### Data Structures

MATLAB makes it really simple to generate and work with 2D arrays. In C++ it is much Once they are created, working with 2D arrays in C++ is similar to MATLAB. Just remember that your indices start at 0, not 1!

### Generating Random Numbers
For these simulations, the Brownian motion is dictated by random numbers drawn on a normal distribution. In MATLAB, we can do this with the function normrnd. In C++, however, getting random numbers takes a bit more work.


### Running Multiple Trials

In MATLAB, parallelizing our loops is really easy. We just use parfor and tweak the parallel pool settings as needed. In C++, we need to We could have used pthreads - the built-in solution for multi-threading - but we decided to use OpenMP because it's more flexible and more intuitive to implement. If you want to know more about OpenMP, there's a book in the lab called 'The OpenMP Common Core'. It's a good read - the first two chapters explain the fundamentals of parallel computing really well.

## Closing Thoughts

So that's the code base! I hope that it is useful to you. If you have any questions, feel free to shoot me an email at ccouto@g.hmc.edu. 

Chris Couto (HMC '23)
