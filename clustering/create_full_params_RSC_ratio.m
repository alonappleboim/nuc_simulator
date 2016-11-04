poly_rates = [0];
tf_rates = [0.0001]; % this is the TF assembly rate, and the eviction will be half of this
tf_evic_eff = [0.0001];
rsc_length = [20, 40, 45, 50, 65]; % this is the eviction length of RSC - the sliding length will be twice this
rsc_evic_eff = [0.02, 0.04, 0.07, 0.1, 0.2];
rsc_slide_eff = [10, 20, 40, 70, 100]; % the ratio between slide and eviction

params = combvec(poly_rates, tf_rates, tf_evic_eff, rsc_length, rsc_evic_eff, rsc_slide_eff);