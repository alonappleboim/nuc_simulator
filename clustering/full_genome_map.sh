#!/bin/bash

for ((i=1 ; i<=6648 ; i++));
do
	sbatch -D /cs/bd/Daniel/simulations/logs -e log --mem 40000 --priority 2 --mail-type FAIL --mail-user daniel.gissin@mail.huji.ac.il "/cs/bd/Daniel/nuc_simulator/clustering/gene_map.sh" $i $1
done

### memory amount, mail notification on failure

 --mail-type END --mail-user daniel.gissin@mail.huji.ac.il