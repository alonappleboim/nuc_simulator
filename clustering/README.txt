The clustering scripts are based on two steps - map and reduce.

The map step is run through main_map.sh, and it simulataniously runs all the simulations, each saves a .m file in the output folder in /cs/bd/Daniel/simulations/output.
The reduce step goes over all those .m files and finds the best feature (the best parameter set), and then saves a final "result" .m file in the same folder, which we can take and extract the parameters from that file.

In order to change the parameter set, we need to change the following files:
1. update the actual parameters we want to check in create_full_params.
2. change the sim_params indices to reflect the new params in cluster_sim.
3. update the number of iterations that are needed (one for every parameter set) in main_map.sh
4. update the number of iterations and vector lengths in cluster_reduce.

In order to change the gene that is tested, we need to change the gene number in the following files:
1. map.sh
2. main_reduce.sh