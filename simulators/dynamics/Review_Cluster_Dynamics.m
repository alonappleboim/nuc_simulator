%%% This script is used to review the dynamic simulations that were run on
%%% the cluster, for genes that are simulated well at 0 minutes and their 2
%%% hour state is different than theit 0 minute state.

load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\sequences_structure.mat')
load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\gene_params.mat')

genlen = 3500;
TSS = fix(genlen/2);
NFR_pos = [TSS-299 : TSS+150];

timeCoefficients = zeros(1,length(gene_params));
for i=1:length(gene_params)

    gene_id = gene_params(i,1);
    load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_0m_centers.mat')
    
    try
        load(['C:\Users\Daniel\Documents\MATLAB\Friedman Lab\results\dynamic_genome_results\dynamic_results_' num2str(gene_id) '.mat'])
    catch a
        continue
    end
    
    %%% Keep track of the time expansions:
    timeCoefficients(i) = timeExpansion;

    %%% Test whether the time expansions converge to a value:
    %{
    figure;
    plot(timeExpansions, costPerTimeExpansion, 'sb')
    xlabel('Time Expansion Coefficient')
    ylabel('Total Cost')
    title(['Time Expansion Convergence - Gene ' num2str(gene_id)])
    %}
    
    %%% Look at the state at 0m and at 2h for each gene:
    seq = sequences_structure(gene_id,:);
    [ PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites ] = ...
        Extract_Sites_From_Gene(seq, genlen);
    PolyA_Sites = create_gene_buffer(PolyA_Sites, genlen);
    PolyT_Sites = create_gene_buffer(PolyT_Sites, genlen);
    
    figure;
    subplot(2,1,1)
    plot(data_matrix(1,:),'b')
    hold on
    plot(nuc_sum_matrix(1,:) .* sum(data_matrix(1,:)) ./ sum(nuc_sum_matrix(1,:)),'r')
    plot(PolyA_Sites(NFR_pos), 'k')
    plot(PolyT_Sites(NFR_pos), 'm')
    legend('experiment at t=0','simulation','PolyA', 'PolyT')
    title(['sim-experiment states for 0 minutes and 2 hours - gene ' num2str(gene_id)])
    subplot(2,1,2)
    plot(data_matrix(end,:),'b')
    hold on
    plot(nuc_sum_matrix(end,:) .* sum(data_matrix(end,:)) ./ sum(nuc_sum_matrix(end,:)),'r')
    plot(PolyA_Sites(NFR_pos), 'k')
    plot(PolyT_Sites(NFR_pos), 'm')
    legend('experiment at t=2h','simulation','PolyA', 'PolyT')
    
end

% Look at the time coefficients:
figure;
hist(timeCoefficients)
title('Time Coefficients Calculated For The Genes')
xlabel('Time Coefficient')
ylabel('# Of Genes')


