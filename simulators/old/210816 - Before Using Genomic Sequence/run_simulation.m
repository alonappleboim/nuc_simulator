function [s_hist, s_hist_coverage, peak_distance_1_2, peak_distance_2_3] = run_simulation(varargin)

% the function runs a single simulation and returns:
% 	s_hist - the state history of the simulation
% 	s_hist_coverage - a summation of the state history, smoothed with a gaussian of width 20 and normalised (ksdensity)
%	peak_distance_1_2 - the distance between the +1 nuc and the +2 nuc (the distance between the two peaks)
%	peak_distance_2_3 - the distance between the +2 nuc and the +3 nuc (the distance between the two peaks)
%
% the functions possible parameters are:
%	n_steps - the number of steps in the simulation
%	gen_len - the length of the genome in the simulation
%	linker_len - the linker_len parameter (the width of the sigmoid added to the ends of the nucleosome footprint)
%	nuc_width - the width of the nucleosome (before adding the sigmoids)
%	slide_len - the amount of bps in a single slide
% 	poly_rate - the rate of the left sliding that is due to the polymerase
%	poly_pos - a vector indicating the positions of the polymerase (for example: 1200:2800)
%	NFR_pos - a vector indicating the position of the NFR (for example: 1000:1200)


%% %%% PARAMETERS %%%%%

defaults = struct('n_steps', 20000,...
                  'gen_len', 3500,...
				  'nuc_width', 147, ...
				  'slide_len', 10, ...
                  'poly_rate', 3, ...
                  'linker_len', 5, ...
				  'NFR_pos', [1000:1200], ...
				  'poly_pos', [1200:2800] ...
				  );
inputs = parse_namevalue_pairs(defaults, varargin);

n_steps = inputs.n_steps;
gen_len = inputs.gen_len;

params = struct();
params.nuc_width = inputs.nuc_width;
params.slide_len = inputs.slide_len;
params.linker_len = inputs.linker_len;
params.a_rate = ones(1,gen_len); 
params.e_rate = ones(1,gen_len);
params.r_rate = 0.1*ones(1,gen_len); 
params.l_rate = 0.1*ones(1,gen_len); 
params.nuc_footprint = ones(1,(params.nuc_width.*2) - 1);

% make NFRs
params.l_rate(inputs.NFR_pos(1) : inputs.NFR_pos((inputs.NFR_pos(end)-inputs.NFR_pos(1)) / 2)) = 10;
params.r_rate(inputs.NFR_pos((inputs.NFR_pos(end)-inputs.NFR_pos(1)) / 2) : inputs.NFR_pos(end)) = 10;
params.e_rate(inputs.NFR_pos(20) : inputs.NFR_pos(end-20)) = 1;

% make a Polymerase
params.l_rate(inputs.poly_pos(1) : inputs.poly_pos(end)) = inputs.poly_rate;

%% Run the Simulation

[time, s_hist] = gillespie(params, 'n_steps', n_steps);%, 's0', zeros(1,gen_len));

%% %%% FINAL CALCULATIONS %%%%%

% save the state histories:
centers_vector = sum(s_hist(:,:));
s_hist_coverage = ksdensity(1:length(centers_vector),1:length(centers_vector),'weights',double(centers_vector(1:end)),'width',20);

% get the peak distances
[peak_distance_1_2, peak_distance_2_3, tmp1, tmp2] = get_peak_distances(s_hist, inputs.NFR_pos);

