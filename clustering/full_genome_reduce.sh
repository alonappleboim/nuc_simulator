#!/bin/bash

for ((i=1 ; i<=200 ; i++));
do
	/cs/bd/Daniel/nuc_simulator/clustering/full_main_reduce.sh $i &
done

wait
