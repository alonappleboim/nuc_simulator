#!/bin/bash

matlab -nodesktop -nosplash -nodisplay -nojvm -r "load('$2'); load('/cs/bd/Daniel/experiment_data/sequences_structure.mat'); addpath(genpath('/cs/bd/Daniel/nflab_scripts')); addpath(genpath('/cs/bd/Daniel/nuc_simulator')); cluster_null($1,sequences_structure,data,'$3'); quit;"
