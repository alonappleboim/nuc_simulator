#!/bin/bash

matlab -nodesktop -nosplash -nodisplay -nojvm -r "load('/cs/bd/Daniel/experiment_data/wt_centers.mat'); load('/cs/bd/Daniel/experiment_data/sequences_structure.mat'); addpath(genpath('/cs/bd/Daniel/nflab_scripts')); addpath(genpath('/cs/bd/Daniel/nuc_simulator')); [nuc_sum,feature]=cluster_sim($1,6641,sequences_structure,wt_3h); quit;" &

wait


