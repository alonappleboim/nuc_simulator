function full_cluster_sim(params_index , gene_index, sequences_structure, wt_3h)

load_bad_genes;
if (~isempty(find(bad_genes(bad_genes==gene_index),1)))
    return
end

genlen = 3500;
TSS = fix(genlen/2);
NFR_pos = [TSS-299 : TSS+150];

% load all the necessary data:
addpath(genpath('/cs/bd/Daniel/nflab_scripts'));
addpath(genpath('/cs/bd/Daniel/nuc_simulator'));
seq = sequences_structure(gene_index,:);
wt_data = wt_3h(gene_index,:);

if (isnan(wt_data(1)))
    nuc_sum = 0;
    likelihood = nan;
    plus1_dist = nan;
    minus1_dist = nan;
    peak_num_delta = nan;
    plus1width = nan;
    minus1width = nan; 
    height_ratio = nan;

else
    
    % make the wt data the right length:
    buffer = genlen - 2501;
    right_buffer = fix((buffer-500)/2);
    left_buffer = right_buffer + 500;
    if (right_buffer + left_buffer < buffer)
        left_buffer = left_buffer + 1;
    end
    wt_data = [zeros(1,left_buffer), wt_data, zeros(1,right_buffer)];

    % create the full parameter matrix
    create_full_params_301016;

    % choose this specific sim parameters:
    sim_params = params(: , params_index);

    % the combvec params, for reference:
    % combvec(poly_rates, tf_rates, tf_evic_eff, rsc_length, rsc_evic, rsc_slide)

    % create ten simulations that will be combined:
	nuc_sum = zeros(1,genlen);
	for i = 1:20
    nuc_sum1 = run_simulation_from_genome(seq,'report',0, ...
        'poly_rate',sim_params(1), ...
        'REB1_a_rate', sim_params(2), 'REB1_e_rate', (sim_params(2)/2), ...
        'ABF1_a_rate', sim_params(2), 'ABF1_e_rate', (sim_params(2)/2), ...
        'RAP1_a_rate', sim_params(2), 'RAP1_e_rate', (sim_params(2)/2), ...
        'TF_evic_intensity', sim_params(3), ...
        'RSC_evic_length', sim_params(4), 'RSC_slide_length', sim_params(4).*2, ...
        'RSC_evic_intensity', sim_params(5), ...
        'RSC_slide_intensity', sim_params(6), 'slide_len', 3, 'gen_len', genlen, 'n_steps', 10000);
	nuc_sum = nuc_sum + nuc_sum1;
	end
	
    % get the feature of the simulation:
    [likelihood, plus1_dist, minus1_dist, peak_num_delta, plus1width, minus1width, height_ratio] = ...
        Compare_Sum_To_Data(nuc_sum, wt_data, NFR_pos, true);
    
end

% save the data to a .mat file:
save(['/cs/bd/Daniel/simulations/full_output_wt_301016/sim_' num2str(params_index) 'gene_' num2str(gene_index) '.mat'] , ...
	'nuc_sum', 'likelihood', 'plus1_dist', 'minus1_dist','peak_num_delta','plus1width','minus1width','height_ratio');
