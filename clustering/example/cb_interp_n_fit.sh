#!/bin/bash
cd ~/Dropbox/MATLAB/chip_analysis/dynamics/cluster_run/
for i in {1..67}
do
   matlab -nodesktop -nosplash -nodisplay -nojvm -r "bi=$i;BS=1000; run cluster_interpandfit.m ; quit;" &
done

