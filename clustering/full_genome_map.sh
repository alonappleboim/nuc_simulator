#!/bin/bash

for ((i=1 ; i<=500 ; i++));
do
	sbatch -D /cs/bd/Daniel/simulations/logs -e log --mem 10000 --priority 2 "/cs/bd/Daniel/nuc_simulator/clustering/gene_map.sh" $i
done
