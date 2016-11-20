function cluster_reduce(gene_index, all_data, results_path)

% change the gene id for avital's simulation:
testtest;
gene_index = genes(gene_index);

% create the full parameter matrix
create_params_avital;
num_of_runs = length(params(1,:));

% ignore NaN genes:
exp_data = all_data(gene_index,:);
if (isnan(exp_data(1)))
    return
end

% create a matrix for all the results:
features = zeros(num_of_runs, 1);
likelihoods = zeros(num_of_runs, 1);
nuc_sums = zeros(num_of_runs, 3500);

% load all of the results:
for i = 1:num_of_runs
	try
        load([results_path 'sim_' num2str(i) '_gene_' num2str(gene_index) '.mat']);
    catch a
        features(i,1) = nan;
        likelihoods(i,1) = nan;
        nuc_sums(i,:) = zeros(1,3500);
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
    features(i, 1) = likelihood - plus1_dist - minus1_dist;
	nuc_sums(i, :) = nuc_sum;
    
    % TODO - make the features be the ratio between the sim and the uniform
    % distribution?
end

% find the best result:
[best_sim_feature, best_sim_index] = max(features);
[best_likelihood, best_likelihood_index] = max(likelihoods);
nuc_sum_feature = nuc_sums(best_sim_index, :);
nuc_sum_likelihood = nuc_sums(best_likelihood_index, :);

% save to a new .mat file:
save([results_path 'results_' num2str(gene_index)] , ...
	'best_sim_feature', 'best_sim_index', 'best_likelihood', 'best_likelihood_index', 'features', 'nuc_sum_feature', 'nuc_sum_likelihood');