genes = find(~isnan(data(:,1)));

tf_evic_eff = [0.0001];
RSC_evic_length = [40, 50, 60, 70, 80]; % this is the eviction length of RSC - the sliding length will be twice this
RSC_evic_slide_length_ratio = [1.5, 2];
rsc_evic_eff = [0.001, 0.003, 0.005, 0.01, 0.02, 0.04];
rsc_evic_slide_eff_ratio = [1, 2, 4, 8, 16]; % the ratio between slide and eviction

params = combvec(tf_evic_eff, RSC_evic_length, RSC_evic_slide_length_ratio, rsc_evic_eff, rsc_evic_slide_eff_ratio);