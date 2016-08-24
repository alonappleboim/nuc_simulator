#!/bin/bash
cd /cs/bd/Daniel/simulations
for i in {1..10}
do
   matlab -nodesktop -nosplash -nodisplay -nojvm -r "a=$i; save('mat$i.mat','a'); quit;" &
done