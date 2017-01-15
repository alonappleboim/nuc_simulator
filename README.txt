Overview Of The Project Codebase:
=================================
The project contains code for both running different simulations, reviewing the results of those simulations and visualizing the results as well. The code is split into three folders:
1. simulators - containing the main code for running and reviewing simulations.
2. clustering - containing the code for running the simulations on the cluster.
3. visualization - containing code for making movies out of simulations.

Each of the folders has their own README file inside, specifing how to run that part of the project.


The General Kinds Of Simulations In The Project:
================================================
There are two main kinds of simulations that were used so far in the project - static simulations and dynamic simulations. 
The static simulation is the basic simulation. the assembly, eviction and sliding rates are generated from the DNA sequence of the given gene, and then those rates are fed into the simulator, which returns a nucleosome position distribution which can be compared to the experimental data.
The dynamic simulation is used to compare the change in time of the nucleosome positions of a certain gene. A gene is fed into the simulation, along with its best parameters (that were found earlier), and the simulation runs a static simulation and then the final state is fed as an initial state into a new static simulation with null parameters (representing a cell with no RSC). The resulting state history can be used to see how the nucleosome positions change over time when RSC is damaged.