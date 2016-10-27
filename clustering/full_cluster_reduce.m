function full_cluster_reduce(gene_index)

% load all the necessary data:
load('/cs/bd/Daniel/nuc_simulator/clustering/experiment_data/wt_centers.mat');
load('/cs/bd/Daniel/experiment_data/sequences_structure.mat');
addpath(genpath('/cs/bd/Daniel/nflab_scripts'));
addpath(genpath('/cs/bd/Daniel/nuc_simulator'));

% create the full parameter matrix
create_full_params_271016;
num_of_runs = length(params(1,:));

% create a matrix for all the results:
features = zeros(num_of_runs, 1);
nuc_sums = zeros(num_of_runs, 3500);

% load all of the results:
for i = 1:num_of_runs
	try
        load(['/cs/bd/Daniel/simulations/full_output_6h/sim_' num2str(i) 'gene_' num2str(gene_index) '.mat']);
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
    
    features(i, 1) = likelihood - plus1_dist - minus1_dist;
	nuc_sums(i, :) = nuc_sum;
end

% find the best result:
[best_sim_feature, best_sim_index] = max(features);

% save to a new .mat file:
save(['/cs/bd/Daniel/simulations/full_output_wt_271016/results_' num2str(gene_index)] , ...
	'best_sim_feature', 'best_sim_index', 'features', 'nuc_sums');