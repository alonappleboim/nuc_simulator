function cluster_null( gene_index, sequences_structure, all_data, results_path)

% create the full parameter matrix
create_params_genome;
gene_index = genes(gene_index);

genlen = 3500;
TSS = fix(genlen/2);
NFR_pos = [TSS-299 : TSS+150];

% get the relevant gene data:
seq = sequences_structure(gene_index,:);
exp_data = all_data(gene_index,:);

if (isnan(exp_data(1)))
    nuc_sum = 0;
    likelihood = nan;
    plus1_dist = nan;
    minus1_dist = nan;
    peak_num_delta = nan;
    plus1width = nan;
    minus1width = nan; 
    height_ratio = nan;
    optimum_likelihood = nan;
    bad_likelihood = nan;

else
    
    % make the wt data the right length:
    buffer = genlen - 2501;
    right_buffer = fix((buffer-500)/2);
    left_buffer = right_buffer + 500;
    if (right_buffer + left_buffer < buffer)
        left_buffer = left_buffer + 1;
    end
    exp_data = [zeros(1,left_buffer), exp_data, zeros(1,right_buffer)];

    % create 20 simulations that will be combined:
	nuc_sum = zeros(1,genlen);
	for i = 1:20
        nuc_sum1 = run_null_model_simulation_from_genome(seq,'report',0, ...
        'poly_rate',0, ...
        'slide_len', 3, 'gen_len', genlen, 'n_steps', 10000);
	nuc_sum = nuc_sum + nuc_sum1;
	end
	
    % get the feature of the simulation:
    [best_likelihood, plus1_dist, minus1_dist, peak_num_delta, plus1width, minus1width, height_ratio] = ...
        Compare_Sum_To_Data(nuc_sum, exp_data, NFR_pos, true);
    [optimum_likelihood, ~,~,~,~,~,~] = ...
        Compare_Sum_To_Data(exp_data, exp_data, NFR_pos, true);
    [bad_likelihood,~,~,~,~,~,~] = ...
        Compare_Sum_To_Data(ones(size(exp_data)), exp_data, NFR_pos, true);
    best_ratio = (best_likelihood - bad_likelihood) / (optimum_likelihood - bad_likelihood);
    nuc_sum_likelihood = nuc_sum;    
    
end

% save the data to a .mat file:
save([results_path 'results_' data_type '_' num2str(gene_index)] , ...
	'best_likelihood', 'best_ratio', 'nuc_sum_likelihood');