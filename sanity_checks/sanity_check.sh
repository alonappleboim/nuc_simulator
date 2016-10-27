#!/bin/bash

for ((i=1 ; i<=10 ; i++));
do

	for ((j=1 ; j<=10 ; j++));
	do
		matlab -nodesktop -nosplash -nodisplay -nojvm -r "load('/cs/bd/Daniel/nuc_simulator/clustering/experiment_data/wt_centers.mat'); load('/cs/bd/Daniel/experiment_data/sequences_structure.mat'); addpath(genpath('/cs/bd/Daniel/nflab_scripts')); addpath(genpath('/cs/bd/Daniel/nuc_simulator')); sanity_check($1,sequences_structure,wt_3h,$i,$j); quit;" &
	done

	wait
	
done


