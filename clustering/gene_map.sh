#!/bin/bash

for ((j=0 ; j<=9 ; j++));
do

	for ((i=0 ; i<=9 ; i++));
	do
		let "k=$j*10+$i+1"
		matlab -nodesktop -nosplash -nodisplay -nojvm -r "load('/cs/bd/Daniel/nuc_simulator/clustering/experiment_data/sth1_6h_centers.mat'); load('/cs/bd/Daniel/experiment_data/sequences_structure.mat'); addpath(genpath('/cs/bd/Daniel/nflab_scripts')); addpath(genpath('/cs/bd/Daniel/nuc_simulator')); full_cluster_sim($k,$1,sequences_structure,sth1_6h); quit;" &
	done

	wait
	
done


