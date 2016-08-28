#!/bin/bash

for i in {1..1800}
do
	sbatch -D /cs/bd/Daniel/simulations -e log --mem 2000 "map.sh" $i
done

exit 0