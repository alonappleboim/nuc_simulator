poly_rates = [0];
tf_rates = [0.001]; % make the evic be half of the assem
tf_evic_eff = [0.001];
rsc_length = [10, 20, 30, 40, 50]; % make the slide twice as long as the evic
rsc_evic = [0.001];
rsc_slide = [0.001, 1, 2, 4, 6, 10];

params = combvec(poly_rates, tf_rates, tf_evic_eff, rsc_length, rsc_evic, rsc_slide);