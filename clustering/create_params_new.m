tf_evic_eff = [0.0001];
RSC_evic_length = [60, 70, 80]; % this is the eviction length of RSC - the sliding length will be twice this
RSC_evic_slide_length_ratio = [0.5, 1, 1.5, 2, 2.5];
rsc_evic_eff = [0.02, 0.04, 0.07, 0.1, 0.2, 0.4];
rsc_evic_slide_eff_ratio = [1, 4, 32]; % the ratio between slide and eviction

params = combvec(poly_rates, tf_rates, tf_evic_eff, rsc_length, rsc_evic_eff, rsc_slide_eff);