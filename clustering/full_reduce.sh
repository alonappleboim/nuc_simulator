#!/bin/bash

matlab -nodesktop -nosplash -nodisplay -nojvm -r "addpath(genpath('/cs/bd/Daniel/nflab_scripts')); addpath(genpath('/cs/bd/Daniel/nuc_simulator')); full_cluster_reduce($1,$2);  quit;"
