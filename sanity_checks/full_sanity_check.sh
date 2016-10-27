#!/bin/bash

for ((i=1 ; i<=153 ; i++));
do
	sbatch -D /cs/bd/Daniel/simulations/logs -e log --mem 10001 "/cs/bd/Daniel/nuc_simulator/sanity_checks/sanity_check.sh" $i
done
