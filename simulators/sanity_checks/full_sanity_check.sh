#!/bin/bash

for ((i=1 ; i<=80 ; i++));
do
	sbatch -D /cs/bd/Daniel/simulations/logs -e log --mem 10001 -t 23:00:00 "/cs/bd/Daniel/nuc_simulator/sanity_checks/sanity_check.sh" $i
done
