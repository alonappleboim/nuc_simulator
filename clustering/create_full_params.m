poly_rates = [0];
tf_rates = [0.0001]; % this is the TF assembly rate, and the eviction will be half of this
tf_evic_eff = [0.0001];
rsc_length = [10, 15, 20]; % this is the eviction length of RSC - the sliding length will be twice this
rsc_evic_eff = [0.0001, 0.1, 0.2];
rsc_slide_eff = [0.0001, 1, 2, 4];

params = combvec(poly_rates, tf_rates, tf_evic_eff, rsc_length, rsc_evic_eff, rsc_slide_eff);