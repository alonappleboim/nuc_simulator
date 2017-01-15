function sanity_reduce(gene_index)

% load all the necessary data:
load('/cs/bd/Daniel/nuc_simulator/clustering/experiment_data/wt_centers.mat');
load('/cs/bd/Daniel/experiment_data/sequences_structure.mat');
addpath(genpath('/cs/bd/Daniel/nflab_scripts'));
addpath(genpath('/cs/bd/Daniel/nuc_simulator'));

features = zeros(10,10);

% load all of the results:
for i = 1:10
    for l = 1:10
        try
            load(['/cs/bd/Daniel/simulations/sanity_checks/sim_' num2str(i) '_' num2str(l) 'gene_' num2str(gene_index) '.mat']);
        catch a
            features(i,l) = nan;
            continue
        end

        if (length(nuc_sum) == 1) % in case the gene was a NaN gene...
            return
        end

        features(i,l) = likelihood;
    end
end

% save to a new .mat file:
save(['/cs/bd/Daniel/simulations/sanity_checks/sanity_results_' num2str(gene_index)] , ...
	'features');