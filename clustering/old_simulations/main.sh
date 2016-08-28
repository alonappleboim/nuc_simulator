#!/bin/bash

matlab -nodesktop -nosplash -nodisplay -nojvm -r "a=1; save('/cs/bd/Daniel/simulations/output/mat.mat','a'); quit;"

for i in {1..1800}
do
	sbatch -D /cs/bd/Daniel/simulations -e log --mem 2000 "simulate.sh" $i
done

exit 0