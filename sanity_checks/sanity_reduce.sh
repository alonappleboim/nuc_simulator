#!/bin/bash

matlab -nodesktop -nosplash -nodisplay -nojvm -r "load('/cs/bd/Daniel/nuc_simulator/clustering/experiment_data/wt_centers.mat'); load('/cs/bd/Daniel/experiment_data/sequences_structure.mat'); addpath(genpath('/cs/bd/Daniel/nflab_scripts')); addpath(genpath('/cs/bd/Daniel/nuc_simulator')); sanity_reduce($1); quit;"
