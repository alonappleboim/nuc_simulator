#!/bin/bash
cd ~/Dropbox/MATLAB/chip_analysis/dynamics/cluster_run/
matlab -nodesktop -nosplash -nodisplay -nojvm -r "bi=$SGE_TASK_ID;BS=1000; run cluster_interpandfit.m ; quit;"
