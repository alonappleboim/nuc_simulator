#!/bin/bash

for ((j=0 ; j<=$2 ; j++));
do

	for ((i=0 ; i<=9 ; i++));
	do
		let "k=$j*10+$i+1" 
		matlab -nodesktop -nosplash -nodisplay -nojvm -r "load('$3'); load('/cs/bd/Daniel/experiment_data/sequences_structure.mat'); addpath(genpath('/cs/bd/Daniel/nflab_scripts')); addpath(genpath('/cs/bd/Daniel/nuc_simulator')); cluster_map($1,$k,sequences_structure,data,'$4'); quit;" &
	done

	wait
	
done

matlab -nodesktop -nosplash -nodisplay -nojvm -r "load('$3'); load('/cs/bd/Daniel/experiment_data/sequences_structure.mat'); addpath(genpath('/cs/bd/Daniel/nflab_scripts')); addpath(genpath('/cs/bd/Daniel/nuc_simulator')); cluster_reduce($1,data,'$4');  quit;"