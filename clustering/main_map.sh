#!/bin/bash

echo $1
echo $2

for (i=1 ; i<=$2 ; i++);
do
	sbatch -D /cs/bd/Daniel/simulations/logs -e log --mem 2000 "/cs/bd/Daniel/nuc_simulator/clustering/map.sh" $i $1
done