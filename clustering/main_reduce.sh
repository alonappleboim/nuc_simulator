#!/bin/bash

sbatch -D /cs/bd/Daniel/simulations/logs -e log --mem 2000 "/cs/bd/Daniel/nuc_simulator/clustering/reduce.sh" $1 $2