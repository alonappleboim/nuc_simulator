#!/bin/bash

for ((j=0 ; j<=59 ; j++));
do

	for ((i=0 ; i<=9 ; i++));
	do
		let "k=$j*10+$i+1" 
		matlab -nodesktop -nosplash -nodisplay -nojvm -r "load('/cs/bd/Daniel/nuc_simulator/clustering/experiment_data/wt_centers.mat'); load('/cs/bd/Daniel/experiment_data/sequences_structure.mat'); addpath(genpath('/cs/bd/Daniel/nflab_scripts')); addpath(genpath('/cs/bd/Daniel/nuc_simulator')); full_cluster_sim($k,$1,sequences_structure,wt_3h); quit;" &
	done

	wait
	
done


