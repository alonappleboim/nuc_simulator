#!/bin/bash

for i in {1..1800}
do
	sbatch -D /cs/bd/Daniel/simulations/logs -e log --mem 2000 "/cs/bd/Daniel/nuc_simulator/clustering/map.sh" $i
done