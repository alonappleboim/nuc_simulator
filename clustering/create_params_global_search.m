genes = [97,93,90,89,86,84,83,82,80,79,77,76,74,73,61,56,46,40,38,30,29,28,24,23,22,16,9,195,193,190,189,188,187,185,184,183,180,178,177,176,173,170,169,168,161,160,155,152,150,149,148,144,138,117];

tf_evic_eff = [0.0001];
RSC_evic_length = [60, 70, 80, 100, 120]; % this is the eviction length of RSC - the sliding length will be twice this
RSC_evic_slide_length_ratio = [1, 1.5, 2, 2.5, 3];
rsc_evic_eff = [0.02, 0.04, 0.07, 0.1, 0.2, 0.4];
rsc_evic_slide_eff_ratio = [2, 4, 8, 32]; % the ratio between slide and eviction

params = combvec(tf_evic_eff, RSC_evic_length, RSC_evic_slide_length_ratio, rsc_evic_eff, rsc_evic_slide_eff_ratio);