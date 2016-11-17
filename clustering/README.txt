The clustering scripts are based on one main script - genome_sim.
genome_sim then calls gene_sim for every gene we want to simulate, each getting its own job on the cluster.

genome_sim gets 5 parameters:
1. lower bound gene
2. upper bound gene
3. number of sims per gene (calculated from create_params - the product of the vector lengths)
4. the path of the data containing the .mat file with the experimental data
5. the output path of the result files (MUST END IN "/")

example for genome_sim for the genes 20-53:
genome_sim.sh 20 53 120 /cs/bd/Daniel/nuc_simulator/clustering/experiment_data/wt_centers.mat /cs/bd/Daniel/simulations/OUTPUT/

The only thing we need to touch in order to run our simulation, is the parameters in the create_params function, and to make sure that that is the function that is called in cluster_sim.m and cluster_reduce.m