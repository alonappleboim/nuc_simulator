function [peak_distance, peak_ratio, simulation] = simulate_polymerase(poly_rate)

%% %%% PARAMETERS %%%%%

n_steps = 20000;
params = struct();
gen_len = 6000;
params.nuc_width = 147;
params.slide_len = 10;
params.a_rate = ones(1,gen_len); 
params.e_rate = ones(1,gen_len);
params.r_rate = 0.1*ones(1,gen_len); 
params.l_rate = 0.1*ones(1,gen_len); 
params.nuc_footprint = ones(1,(params.nuc_width.*2) - 1);
params.linker_len = 5;

% NFRs
params.r_rate(1050:1200) = 30;
params.l_rate(1000:1050) = 30;

% Polymerase Parameters
params.l_rate(1200:3200) = poly_rate;

%%

[time, s_hist] = gillespie(params, 'n_steps', n_steps, 's0', zeros(1,gen_len));

%% %%% FINAL CALCULATIONS %%%%%

% get the number of times each bp had a nuceosome center on it
centers_vector = sum(s_hist(:,161:2660));
smooth_vector = ksdensity(1000:length(centers_vector)-1,1000:length(centers_vector)-1,'weights',double(centers_vector(1000:end-1)),'width',20);
smooth_vector_ratio = ksdensity(1000:1800,1000:1800,'weights',double(centers_vector(1000:1800)),'width',20);

% extract the desired parameters from the simulation:
peak_distance = get_average_peak_distance(smooth_vector);
peak_ratio = get_peak_decline_ratio(smooth_vector_ratio);
simulation = smooth_vector;