function cluster_reduce(gene_index)

% load all the necessary data:
load('/cs/bd/Daniel/experiment_data/wt_centers.mat');
load('/cs/bd/Daniel/experiment_data/sequences_structure.mat');
addpath(genpath('/cs/bd/Daniel/nflab_scripts'));
addpath(genpath('/cs/bd/Daniel/nuc_simulator'));

% create the full parameter matrix
create_full_params;

% create a matrix for all the results:
features = zeros(1800, 1);
nuc_sums = zeros(1800, 3500);

% load all of the results:
for i = 1:1800
	load(['/cs/bd/Daniel/simulations/output/sim_' num2str(i) 'gene_' num2str(gene_index) '.mat']);
	features(i) = feature;
	nuc_sums(i, :) = nuc_sum;
end

% find the best result:
[best_sim_feature, best_sim_index] = max(features);

% save to a new .mat file:
save(['/cs/bd/Daniel/simulations/output/results_' num2str(gene_index)] , ...
	'best_sim_feature', 'best_sim_index', 'features', 'nuc_sums');