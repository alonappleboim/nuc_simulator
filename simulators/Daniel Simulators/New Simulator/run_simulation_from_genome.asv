function [full_s_hist, nuc_s_hist, peak_distance_1_2, peak_distance_2_3, NFR_width] = run_simulation_from_genome(genome, varargin)
%run_simulation_from_genome the function for running a single simulation from genomic data.
%
% the function runs a single simulation and returns:
% 	full_s_hist - the full state history of the simulation (includes trans factors)
% 	nuc_s_hist - the nucleosome state history.
%	peak_distance_1_2 - the distance between the +1 nuc and the +2 nuc (the distance between the two peaks)
%	peak_distance_2_3 - the distance between the +2 nuc and the +3 nuc (the distance between the two peaks)
%   NFR_width - the average width of the NFR along the simulation
%
% the function accepts a 2501-bp-long genome string, and additional
% optional parameters.
% the functions optional parameters are:
%	n_steps - the number of steps in the simulation
%	gen_len - the length of the genome in the simulation
%	linker_len - the linker_len parameter (the width of the sigmoid added to the ends of the nucleosome footprint)
%	nuc_width - the width of the nucleosome (before adding the sigmoids)
%   REB1_width - the width of the REB1 trans factor
%   ABF1_width - the width of the ABF1 trans factor
%   RAP1_width - the width of the RAP1 trans factor
%	slide_len - the amount of bps in a single slide
% 	poly_rate - the rate of the left sliding that is due to the polymerase
%	poly_pos - a vector indicating the positions of the polymerase (for example: 1200:2800)
%	NFR_pos - a vector indicating the position of the NFR (for example: 1000:1200)
%   REB1_a_rate - the REB1 assembly rate
%   REB1_e_rate - the REB1 eviction rate
%   ABF1_a_rate - the ABF1 assembly rate
%   ABF1_e_rate - the ABF1 eviction rate
%   RAP1_a_rate - the RAP1 assembly rate
%   RAP1_e_rate - the RAP1 eviction rate

%% %%% PARAMETERS %%%%%

defaults = struct('n_steps', 20000,...
                  'gen_len', 3500,...
				  'nuc_width', 147, ...
                  'REB1_width', 10, ...
				  'ABF1_width', 10, ...
				  'RAP1_width', 10, ...
				  'slide_len', 10, ...
                  'linker_len', 5, ...
                  'poly_rate', 0, ...
				  'poly_pos', 1000:2600, ...
                  'REB1_a_rate', 0, ...
                  'REB1_e_rate', 0, ...
                  'ABF1_a_rate', 0, ...
                  'ABF1_e_rate', 0, ...
                  'RAP1_a_rate', 0, ...
                  'RAP1_e_rate', 0 ...
				  );
inputs = parse_namevalue_pairs(defaults, varargin);

n_steps = inputs.n_steps;
gen_len = inputs.gen_len;

params = struct();
params.nuc_width = inputs.nuc_width;
params.slide_len = inputs.slide_len;
params.linker_len = inputs.linker_len;
params.a_rate = ones(1,gen_len); 
params.nuc_footprint = ones(1,(params.nuc_width.*2) - 1);
params.REB1_width = inputs.REB1_width;
params.ABF1_width = inputs.ABF1_width;
params.RAP1_width = inputs.RAP1_width;

% extract the relevant sites from the genome:
[ PolyA_sites, PolyT_sites, REB1_sites, ABF1_sites, RAP1_sites ] = ...
    Extract_Sites_From_Gene(genome);

% set the trans factors rate vectors:
params.REB1_a_rate = REB1_sites .* inputs.REB1_a_rate;
REB1_sites(REB1_sites > 0) = 1; % make all sites equal 1 for equal evic rate
params.REB1_e_rate = REB1_sites .* inputs.REB1_e_rate;
params.ABF1_a_rate = ABF1_sites .* inputs.ABF1_a_rate;
ABF1_sites(ABF1_sites > 0) = 1; % make all sites equal 1 for equal evic rate
params.ABF1_e_rate = ABF1_sites .* inputs.ABF1_e_rate;
params.RAP1_a_rate = RAP1_sites .* inputs.RAP1_a_rate;
RAP1_sites(RAP1_sites > 0) = 1; % make all sites equal 1 for equal evic rate
params.RAP1_e_rate = RAP1_sites .* inputs.RAP1_e_rate;

%%% ADD THE EXTRA PARAMETERS FOR THE FUNCTION
% set the nucleosome rate vectors:
[params.e_rate, params.r_rate, params.l_rate] = ...
    generate_nuc_rates_from_sites(PolyA_sites, PolyT_sites, REB1_sites, ABF1_sites, RAP1_sites);

% make the Polymerase
%{
params.l_rate(inputs.poly_pos(1) : inputs.poly_pos(end)) = ...
    params.l_rate(inputs.poly_pos(1) : inputs.poly_pos(end)) + inputs.poly_rate;
%}

%% Run the Simulation

[time, full_s_hist, nuc_s_hist] = Simulate(params, 'n_steps', n_steps, 's0', zeros(1,gen_len));

%% %%% EXTACT FEATURES %%%%%

% save the state histories:
centers_vector = sum(nuc_s_hist(:,:));

% get the peak distances
[peak_distance_1_2, peak_distance_2_3, tmp1, tmp2] = get_peak_distances(s_hist, inputs.NFR_pos);

% get the NFR width:
NFR_width

