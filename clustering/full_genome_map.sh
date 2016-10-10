#!/bin/bash

for ((i=1 ; i<=6648 ; i++));
do
	sbatch -D /cs/bd/Daniel/simulations/logs -e log --mem 24000 --priority 2 "/cs/bd/Daniel/nuc_simulator/clustering/gene_map.sh" $i
done
