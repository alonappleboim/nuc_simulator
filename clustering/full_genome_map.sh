#!/bin/bash

for ((i=1 ; i<=6648 ; i++));
do
	sbatch -D /cs/bd/Daniel/simulations/logs -e log --mem 2000 --priority 2 "/cs/bd/Daniel/nuc_simulator/clustering/gene_map.sh" $1 $i
done
