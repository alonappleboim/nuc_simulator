#!/bin/bash

for ((i=$1 ; i<=$2 ; i++));
do
	sbatch -D /cs/bd/Daniel/simulations/logs -e log$i --mem 10001 -t 200:00:00 "/cs/bd/Daniel/nuc_simulator/clustering/gene_sim_null.sh" $i $3 $4
done