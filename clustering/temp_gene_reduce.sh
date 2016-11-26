#!/bin/bash

matlab -nodesktop -nosplash -nodisplay -nojvm -r "load('$3'); load('/cs/bd/Daniel/experiment_data/sequences_structure.mat'); addpath(genpath('/cs/bd/Daniel/nflab_scripts')); addpath(genpath('/cs/bd/Daniel/nuc_simulator')); cluster_reduce_global_search($1,data,'$4');  quit;"