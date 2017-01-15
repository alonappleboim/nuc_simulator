#!/bin/bash

matlab -nodesktop -nosplash -nodisplay -nojvm -r "load('/cs/bd/Daniel/experiment_data/sequences_structure.mat'); addpath(genpath('/cs/bd/Daniel/nflab_scripts')); addpath(genpath('/cs/bd/Daniel/nuc_simulator')); cluster_dynamic_simulation($1,'$2'); quit;"