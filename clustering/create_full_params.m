poly_rates = [0,0.5,1];
tf_rates = [0.0001]; % this is the TF assembly rate, and the eviction will be half of this
tf_evic_eff = [0.0001];
rsc_length = [10, 15, 20, 25, 30, 40, 50, 60]; % this is the eviction length of RSC - the sliding length will be twice this
rsc_evic_eff = [0.0001, 0.1, 0.2, 0.4, 0.6];
rsc_slide_eff = [0.0001, 1, 2, 4, 6];

params = combvec(poly_rates, tf_rates, tf_evic_eff, rsc_length, rsc_evic_eff, rsc_slide_eff);