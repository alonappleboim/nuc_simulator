#!/bin/bash

for ((i=1 ; i<=$2 ; i++));
do
	matlab -nodesktop -nosplash -nodisplay -nojvm -r "load('/cs/bd/Daniel/experiment_data/wt_centers.mat'); load('/cs/bd/Daniel/experiment_data/sequences_structure.mat'); addpath(genpath('/cs/bd/Daniel/nflab_scripts')); addpath(genpath('/cs/bd/Daniel/nuc_simulator')); [nuc_sum,feature]=cluster_sim($i,$1,sequences_structure,wt_3h); quit;" &

done

wait
