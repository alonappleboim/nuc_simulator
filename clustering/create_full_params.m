poly_rates = [0];
tf_rates = [0.0001]; % this is the TF assembly rate, and the eviction will be half of this
tf_evic_eff = [0.0001];
rsc_length = [10, 20, 40, 60, 80]; % this is the eviction length of RSC - the sliding length will be twice this
rsc_evic_eff = [0.0001, 0.01, 0.05, 0.15];
rsc_slide_eff = [0.0001, 0.5, 1, 2, 4];

params = combvec(poly_rates, tf_rates, tf_evic_eff, rsc_length, rsc_evic_eff, rsc_slide_eff);