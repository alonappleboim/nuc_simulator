#!/bin/bash

matlab -nodesktop -nosplash -nodisplay -nojvm -r "addpath(genpath('/cs/bd/Daniel/nflab_scripts')); addpath(genpath('/cs/bd/Daniel/nuc_simulator')); cluster_reduce(3333);  quit;"
