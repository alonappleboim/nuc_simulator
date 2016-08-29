poly_rates = [0];
tf_rates = [0.0001]; % make the evic be half of the assem
tf_evic_eff = [0.0001];
rsc_length = [10, 15, 20, 25, 30, 40]; % make the slide twice as long as the evic
rsc_evic = [0.0001, 0.1, 0.2];
rsc_slide = [0.0001, 1, 2, 4, 6, 10, 15];

params = combvec(poly_rates, tf_rates, tf_evic_eff, rsc_length, rsc_evic, rsc_slide);