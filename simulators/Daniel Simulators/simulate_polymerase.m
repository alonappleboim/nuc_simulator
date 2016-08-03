function [peak_distance] = simulate_polymerase(l_rate)

%% %%% PARAMETERS %%%%%

n_steps = 20000;
params = struct();
gen_len = 10000;
params.nuc_width = 147;
params.slide_len = 10;
params.a_rate = ones(1,gen_len); 
params.e_rate = ones(1,gen_len);
params.r_rate = 0.1*ones(1,gen_len); 
params.l_rate = 0.1*ones(1,gen_len); 
params.nuc_footprint = ones(1,(params.nuc_width.*2) - 1);
params.linker_len = 5;

% NFRs
params.r_rate(5050:5200) = 30;
params.l_rate(5000:5050) = 30;

% Polymerase Parameters
params.l_rate(5200:end) = l_rate;

%%

[time, s_hist] = gillespie(params, 'n_steps', n_steps, 's0', zeros(1,gen_len));

%% %%% FINAL CALCULATIONS %%%%%

% get the number of times each bp had a nuceosome center on it
centers_vector = sum(s_hist(:,4161:6660));
peak_distance = get_average_peak_distance(centers_vector,25);