function [nuc_sum, time, nuc_s_hist, REB1_s_hist, ABF1_s_hist, RAP1_s_hist] = ... 
        run_simulation_from_genome(genome, varargin)
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
%	gen_len - the length of the genome in the simulation (at least 2501)
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
%   nuc_base_a_rate - the global starting assembly rate of the nucs, before changes
%                     from the binding sites
%   nuc_base_e_rate - the global starting eviction rate of the nucs, before changes
%                     from the binding sites
%   nuc_base_r_rate - the global starting right rate of the nucs, before changes
%                     from the binding sites
%   nuc_base_l_rate - the global starting left rate of the nucs, before changes
%                     from the binding sites
%   TF_evic_intensity - the factor that multiplies the convolution of the
%                       TF sites on the nuc eviction rate
%   RSC_evic_intensity - the factor that multiplies the convolution of the
%                       PolyAT sites on the nuc eviction rate
%   RSC_evic_length - the length of the convolution of the
%                       PolyAT sites on the nuc eviction rate
%   RSC_slide_intensity - the factor that multiplies the convolution of the
%                       PolyAT sites on the nuc sliding rates
%   RSC_slide_length - the length of the convolution of the
%                       PolyAT sites on the nuc sliding rates


%% %%% PARAMETERS %%%%%

defaults = struct('report', nan, ...
                  'n_steps', 25000,...
                  'gen_len', 3500,...
				  'nuc_width', 147, ...
                  'REB1_width', 10, ...
				  'ABF1_width', 14, ...
				  'RAP1_width', 15, ...
				  'slide_len', 5, ...
                  'linker_len', 10, ...
                  'poly_rate', 0, ...
				  'poly_pos', 1100:2500, ...
                  'REB1_a_rate', 0, ...
                  'REB1_e_rate', 0, ...
                  'ABF1_a_rate', 0, ...
                  'ABF1_e_rate', 0, ...
                  'RAP1_a_rate', 0, ...
                  'RAP1_e_rate', 0, ...
                  'nuc_base_a_rate', 0.01, ...
                  'nuc_base_e_rate', 0.01, ...
                  'nuc_base_r_rate', 0.1, ...
                  'nuc_base_l_rate', 0.1, ...
                  'TF_evic_intensity', 0, ...
                  'RSC_evic_intensity', 0.1, ...
                  'RSC_evic_length', 20, ...
                  'RSC_slide_intensity', 4, ...
                  'RSC_slide_length', 40);
p = parse_namevalue_pairs(defaults, varargin);

% set the rates as vectors:
p.REB1_a_rate = p.REB1_a_rate .* ones(1, p.gen_len);
p.REB1_e_rate = p.REB1_e_rate .* ones(1, p.gen_len);
p.ABF1_a_rate = p.ABF1_a_rate .* ones(1, p.gen_len);
p.ABF1_e_rate = p.ABF1_e_rate .* ones(1, p.gen_len);
p.RAP1_a_rate = p.RAP1_a_rate .* ones(1, p.gen_len);
p.RAP1_e_rate = p.RAP1_e_rate .* ones(1, p.gen_len);
p.nuc_base_a_rate = p.nuc_base_a_rate .* ones(1, p.gen_len);
p.nuc_base_e_rate = p.nuc_base_e_rate .* ones(1, p.gen_len);
p.nuc_base_r_rate = p.nuc_base_r_rate .* ones(1, p.gen_len);
p.nuc_base_l_rate = p.nuc_base_l_rate .* ones(1, p.gen_len);

p.nuc_footprint = ones(1,(p.nuc_width.*2) - 1);

% extract the relevant sites from the genome:
[ PolyA_sites, PolyT_sites, REB1_sites, ABF1_sites, RAP1_sites ] = ...
    Extract_Sites_From_Gene(genome, p.gen_len);

% get the rate vectors from the binding sites:
[p.nuc_a_rate, p.nuc_e_rate, p.nuc_r_rate, p.nuc_l_rate,...
    p.REB1_a_rate, p.REB1_e_rate, p.ABF1_a_rate, p.ABF1_e_rate, ...
    p.RAP1_a_rate, p.RAP1_e_rate] = ...
    generate_rates_from_sites(PolyA_sites, PolyT_sites, REB1_sites, ABF1_sites, RAP1_sites, ...
    'nuc_width', p.nuc_width, 'REB1_width', p.REB1_width, 'ABF1_width', p.ABF1_width,...
    'RAP1_width', p.RAP1_width, 'poly_rate', p.poly_rate, 'poly_pos', p.poly_pos, ...
    'REB1_a_rate', p.REB1_a_rate, 'REB1_e_rate', p.REB1_e_rate, 'ABF1_a_rate', p.ABF1_a_rate, ...
    'ABF1_e_rate', p.ABF1_e_rate, 'RAP1_a_rate', p.RAP1_a_rate, 'RAP1_e_rate', p.RAP1_e_rate, ... 
    'nuc_base_a_rate', p.nuc_base_a_rate, 'nuc_base_e_rate', p.nuc_base_e_rate, ... 
    'nuc_base_r_rate', p.nuc_base_r_rate, 'nuc_base_l_rate', p.nuc_base_l_rate, ...
    'TF_evic_intensity', p.TF_evic_intensity, 'RSC_evic_intensity', p.RSC_evic_intensity, ...
    'RSC_evic_length', p.RSC_evic_length, 'RSC_slide_intensity', p.RSC_slide_intensity, ... 
    'RSC_slide_length', p.RSC_slide_length);

%% Run the Simulation

[nuc_sum, time, nuc_s_hist, REB1_s_hist, ABF1_s_hist, RAP1_s_hist] = ...
        Simulate(p, 'n_steps', p.n_steps, 's0', zeros(1,p.gen_len), 'report', p.report);

