#!/bin/bash

for ((i=1 ; i<=80 ; i++));
do
	/cs/bd/Daniel/nuc_simulator/clustering/full_main_reduce.sh $i &
done

wait
