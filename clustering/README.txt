Clustering Codebase Summary:
============================
The code in this folder allows running simulations on the SULFUR cluster, and reviewing results of simulations that were run on the cluster.
The clustering scripts are based on two main scripts - genome_sim.sh (for static simulations) and genome_dynamic.sh (for dynamic simulations).
genome_sim then calls gene_sim for every gene we want to simulate, each getting its own job on the cluster.
genome_dynamic does the same by calling gene_dynamic.
The "experiment_data" folder contains the experiment data from the different timestamps of the experiment, and are loaded in different scripts the are used in the simulator.


genome_sim.sh:
==============
genome_sim gets 5 parameters:
1. lower bound gene
2. upper bound gene
3. number of sims per gene (calculated from create_params - the product of the vector lengths)
4. the path of the data containing the .mat file with the experimental data
5. the output path of the result files (MUST END IN "/")

example for genome_sim for the genes 20-53:
genome_sim.sh 20 53 120 /cs/bd/Daniel/nuc_simulator/clustering/experiment_data/wt_centers.mat /cs/bd/Daniel/simulations/OUTPUT/

The only thing we need to touch in order to run our simulation, is the parameters in the create_params function, and to make sure that that is the function that is called in cluster_sim.m and cluster_reduce.m
A full explanation of running the simulation on the cluster will be presented in a moment.


genome_dynamic.sh:
==================
genome_dynamic get 3 parameters:
1. lower bound gene
2. upper bound gene
3. the output path of the result files (MUST END IN "/")

example for genome_dynamic for the genes 1-45:
genome_dynamic.sh 1 45 /cs/bd/Daniel/simulations/output_dynamic/

The only thing we need to change here in addition to the parameters in the create_params function, is to make sure the gene_params.mat file has the updated mapping of the genes we want to simulate and their optimal parameters for t=0 minutes. A full explanation will follow soon.


=====|||*|||=====


Running genome_sim.sh On The Cluster:
=====================================
The cluster scripts currently are able to test different combinations of parameters and return the ones that best fit the data. All possible combinations of the given parameters are calculated - the script isn't very intelligent at the moment...
In order to run the simulation on the cluster, there are a few things that you need to do:
1. Set the parameters that you want to try
2. Push the changes in git and then pull from the CS computer
3. Run the simulation
4. Copy the result file to your machine and look at the results

Setting The Parameters:
=======================
The parameters are set in the "create_params" script. In it there are a few vectors that are defined, each representing all options for a parameter in the simulations. You need to set those to be what you want. You can create a script of your choosing that sets the required variables, just make sure that the cluster_sim and cluster_reduce functions call your "create_params" script and not some other script.
When you've finished setting the parameters - multiply all of the vector lengths to get the amount of simulations that will run for each gene. You will need to enter that number as the third parameter of the genome_sim.sh script.
Note that the genes indices you give the script are not the gene ids, and that is to avoid NaN genes. What you are giving the script is the gene indices in the "genes" vector (which is created in the create_params script). This gene contains only the genes that are not NaN.

Running The Simulation:
=======================
After connecting to sulfur-gw and setting the parameters you want in the create_params script, all that's left to do is run genome_sim.sh with the parameters (an example can be seen above). The simulations will run and the results will be saved to the path you specified. The result files that are of interest to you are the ones that are named "results_<gene_id>".
If something is wrong and no results are saved (or the results are bad), then you can check the logs in the log path that is specified in the genome_sim.sh file itself (when calling sbatch).

Examining The Results:
======================
Once you have the results file, you can import it to your matlab and then run one of the Review scripts, which will give you a couple of plots of the simulation that had the best likelihood, so you can see what how good the fit was and so on. There are a few Review scripts for different simulations that I ran during the project - you can create your own review scripts by looking at the ones I made - They usually just import the results from the result files in a loop and make a few plots...


=====|||*|||=====


Running genome_dynamic.sh on the cluster:
=========================================
Running this script is not so different from genome_sim.sh - we need to make the same preperations as genome_sim.sh, but also have the gene_params.mat file contain a matrix of the genes you want to simulate and their optimal parameters for t=0 minutes.
The matrix should be called gene_params and be of size(X,2), where X is the number of genes you are simulating. The first column contains the gene ID, and the second contains the optimal parameter index in the create_params script.
Running the simulation and examining the results is done in the same way...