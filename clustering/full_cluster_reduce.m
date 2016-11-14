function full_cluster_reduce(gene_index)

% load all the necessary data:

load('/cs/bd/Daniel/nuc_simulator/clustering/experiment_data/wt_centers.mat');
load('/cs/bd/Daniel/experiment_data/sequences_structure.mat');
addpath(genpath('/cs/bd/Daniel/nflab_scripts'));
addpath(genpath('/cs/bd/Daniel/nuc_simulator'));

wt_data = wt_3h(gene_index,:);
if (isnan(wt_data(1)))
    return
end

% create the full parameter matrix
create_full_params_RSC_ratio;
num_of_runs = length(params(1,:));

% create a matrix for all the results:
features = zeros(num_of_runs, 1);
likelihoods = zeros(num_of_runs, 1);
nuc_sums = zeros(num_of_runs, 3500);

% load all of the results:
for i = 1:num_of_runs
	try
        load(['/cs/bd/Daniel/simulations/full_output_RSC_ratio/sim_' num2str(i) 'gene_' num2str(gene_index) '.mat']);
    catch a
        features(i,1) = nan;
        continue
    end
        
    if (length(nuc_sum) == 1) % in case the gene was a NaN gene...
        return
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
save(['/cs/bd/Daniel/simulations/full_output_RSC_ratio/results_' num2str(gene_index)] , ...
	'best_sim_feature', 'best_sim_index', 'best_likelihood', 'best_likelihood_index', 'features', 'nuc_sum_feature', 'nuc_sum_likelihood');