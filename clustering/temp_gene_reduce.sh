#!/bin/bash

matlab -nodesktop -nosplash -nodisplay -nojvm -r "load('$3'); load('/cs/bd/Daniel/experiment_data/sequences_structure.mat'); addpath(genpath('/cs/bd/Daniel/nflab_scripts')); addpath(genpath('/cs/bd/Daniel/nuc_simulator')); cluster_reduce($1,data,'$4');  quit;"