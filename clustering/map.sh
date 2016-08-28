#!/bin/bash

matlab -nodesktop -nosplash -nodisplay -nojvm -r "load('/cs/bd/Daniel/experiment_data/wt_centers.mat');
load('/cs/bd/Daniel/experiment_data/sequences_structure.mat');
[nuc_sum, feature] = cluster_sim($1, 1001, sequences_structure, wt_3h);  quit;"