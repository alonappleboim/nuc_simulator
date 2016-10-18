#!/bin/bash

for ((i=1 ; i<=153 ; i++));
do
	/cs/bd/Daniel/nuc_simulator/clustering/full_main_reduce.sh $i $1 &
done

wait
