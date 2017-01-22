%%% This script tests the metagene of the experimental data and the
%%% simulation data, to see how well the simulator works when we look at
%%% all the genes together

load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_0m_centers.mat')
genes = find(~isnan(data(:,1)));

genlen = 3500;
TSS = round(genlen/2);
NFR_pos = [TSS-299 : TSS+150];

simulation_data = zeros(1,genlen);
experimental_data = zeros(6648, genlen);

for i = 1:length(data(:,1))
    experimental_data(i,:) = create_gene_buffer(data(i,:), genlen);
end

exp_data = sum(experimental_data(genes,:));

for gene = genes'
    
    try
        load(['C:\Users\Daniel\Documents\MATLAB\Friedman Lab\results\genome_results\results_0m_' num2str(gene) '.mat'])
    catch a
        continue
    end
    simulation_data = simulation_data + nuc_sum_likelihood;
end

% plot the metagene
figure;
plot(exp_data(1200:2200),'b')
hold on
plot(simulation_data(1200:2200) .* sum(exp_data(1200:2200)) ./ sum(simulation_data(1200:2200)),'r')
xlabel('position (TSS at 1000)')
ylabel('intensity')
title('Metage Comparisen Between Simulation And Experiment')
legend('Experiment','Simulation')