function cluster_reduce_global_search(sim_index, all_data, results_path)

% create the full parameter matrix
create_params_global_search;
num_of_runs = 6648;

% create a matrix for all the results:
features = zeros(num_of_runs, 1);
likelihoods = zeros(num_of_runs, 1);
optimal_likelihoods = zeros(num_of_runs, 1);
bad_likelihoods = zeros(num_of_runs, 1);

% load all of the results:
for i = 1:num_of_runs
	try
        load([results_path 'sim_' num2str(sim_index) '_gene_' num2str(i) '.mat']);
    catch a
        features(i,1) = 0;
        likelihoods(i,1) = 0;
        bad_likelihoods(i,1) = 0;
        optimal_likelihoods(i,1) = 1;
        continue
    end
            
    % deal with NaNs:
    if (isnan(plus1_dist))
        plus1_dist = 50;
    end
    if (isnan(minus1_dist))
        minus1_dist = 50;
    end
    
    likelihoods(i,1) = likelihood;
    optimal_likelihoods(i,1) = optimum_likelihood;
    bad_likelihoods(i,1) = bad_likelihood;
    features(i, 1) = likelihood - plus1_dist - minus1_dist;
    
end

ratios = (likelihoods - bad_likelihoods) ./ (optimal_likelihoods - bad_likelihoods);
over_3 = find(ratios > 0.3);
over_2 = find(ratios > 0.2);
sum_ratio = sum(ratios(ratios > 0));

% save to a new .mat file:
save([results_path 'results_' num2str(sim_index)] , ...
	'ratios', 'over_3', 'over_2', 'sum_ratio');