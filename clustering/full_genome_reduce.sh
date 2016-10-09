#!/bin/bash

for ((i=1 ; i<=6648 ; i++));
do
	/cs/bd/Daniel/nuc_simulator/clustering/main_reduce.sh $i $1 &
done

wait
