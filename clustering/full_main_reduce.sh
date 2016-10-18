#!/bin/bash

sbatch -D /cs/bd/Daniel/simulations/logs -e log --mem 200 "/cs/bd/Daniel/nuc_simulator/clustering/full_reduce.sh" $1 $2