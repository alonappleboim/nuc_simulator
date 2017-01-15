function [nuc_sum, time, nuc_s_hist, nuc_evics, REB1_s_hist, ABF1_s_hist, RAP1_s_hist] = ... 
        run_null_model_simulation_from_genome(genome, varargin)
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
				  'slide_len', 3, ...
                  'linker_len', 10, ...
                  's0', zeros(1,3500), ...
                  'poly_rate', 0, ...
				  'poly_pos', 1100:2500, ...
                  'REB1_a_rate', 0.0001, ...
                  'REB1_e_rate', 0.0001, ...
                  'ABF1_a_rate', 0.0001, ...
                  'ABF1_e_rate', 0.0001, ...
                  'RAP1_a_rate', 0.0001, ...
                  'RAP1_e_rate', 0.0001, ...
                  'nuc_base_a_rate', 0.01, ...
                  'nuc_base_e_rate', 0.01, ...
                  'nuc_base_r_rate', 0.1, ...
                  'nuc_base_l_rate', 0.1, ...
                  'TF_evic_intensity', 0, ...
                  'RSC_evic_intensity', 0.0001, ...
                  'RSC_evic_length', 2, ...
                  'RSC_slide_intensity', 0.0001, ...
                  'RSC_slide_length', 2);
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

% extract the null model rates from the genome sequence
dinuc_vec = get_xnucleotid_vector(genome, 4);
dinuc_vec(50:end-50) = conv(dinuc_vec(50:end-50), ones(1,100), 'same');
dinuc_vec(1:49) = 14;
dinuc_vec(end-49:end) = 14;
dinuc_vec = create_gene_buffer(dinuc_vec, p.gen_len);
dinuc_vec(dinuc_vec == 0) = 14;
p.nuc_e_rate = (0.0005064.*(dinuc_vec - 14).^2 + 1) .* p.nuc_base_e_rate;
%p.nuc_e_rate = (0.01.*(dinuc_vec - 14).^2 + 1) .* p.nuc_base_e_rate;
p.nuc_a_rate = p.nuc_base_a_rate;
p.nuc_r_rate = p.nuc_base_r_rate; 
p.nuc_l_rate = p.nuc_base_l_rate;

%% Run the Simulation

[nuc_sum, time, nuc_s_hist, nuc_evics, REB1_s_hist, ABF1_s_hist, RAP1_s_hist] = ...
        Simulate(p, 'n_steps', p.n_steps, 's0', p.s0, 'report', p.report);

