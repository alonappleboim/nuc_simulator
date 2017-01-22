%%% This script will be used to see how many genes we simulate
%%% successfully, and choose the genes with which to test the dynamics
%%% simulations.

load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\sequences_structure.mat')
load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_0m_centers.mat')

genlen = 3500;
TSS = fix(genlen/2);
NFR_pos = [TSS-299 : TSS+150];

create_params_genome;
genes = [1:6648]';

optimal_param_indices = zeros(size(genes));
optimal_likelihoods = zeros(size(genes));
optimal_ratios = zeros(size(genes));

% extract the best parameters for each gene:
for i=1:length(genes)

    gene_id = genes(i);
    
    try
        load(['C:\Users\Daniel\Documents\MATLAB\Friedman Lab\results\genome_results\results_0m_' num2str(gene_id) '.mat'])
    catch a
        optimal_param_indices(i) = nan;
        optimal_likelihoods(i) = nan;
        optimal_ratios(i) = nan;
        continue
    end
    
    optimal_param_indices(i) = best_likelihood_index;
    optimal_likelihoods(i) = best_likelihood;
    optimal_ratios(i) = best_ratio;
        
end

% find the genes that are simulated well (gene + optimal_params):
over_90 = [genes(optimal_ratios > 0.9) , optimal_param_indices(optimal_ratios > 0.9)];
over_80 = [genes(optimal_ratios > 0.8) , optimal_param_indices(optimal_ratios > 0.8)];
over_70 = [genes(optimal_ratios > 0.7) , optimal_param_indices(optimal_ratios > 0.7)];
over_60 = [genes(optimal_ratios > 0.6) , optimal_param_indices(optimal_ratios > 0.6)];
over_50 = [genes(optimal_ratios > 0.5) , optimal_param_indices(optimal_ratios > 0.5)];
over_40 = [genes(optimal_ratios > 0.4) , optimal_param_indices(optimal_ratios > 0.4)];
over_30 = [genes(optimal_ratios > 0.3) , optimal_param_indices(optimal_ratios > 0.3)];
over_20 = [genes(optimal_ratios > 0.2) , optimal_param_indices(optimal_ratios > 0.2)];
over_10 = [genes(optimal_ratios > 0.1) , optimal_param_indices(optimal_ratios > 0.1)];