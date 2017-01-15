function cluster_dynamic_simulation( Gene_id , results_path )
%cluster_dynamic_simulation run a dynamic simulation for a specific gene on
%the cluster, saving the data about the time expansion constant

% load the relevant parameters:
load('/cs/bd/Daniel/nuc_simulator/clustering/gene_params.mat');
load('/cs/bd/Daniel/experiment_data/sequences_structure.mat')
param_index = gene_params(Gene_id, 2);
Gene_id = gene_params(Gene_id,1);
create_params_genome;

%%

iterations = 3;
first_n_steps = 10000;
second_n_steps = 20000;

genlen = 3500;
TSS = round(genlen/2);
NFR_pos = [TSS-299 : TSS+150];

seq = sequences_structure(Gene_id,:);

nuc_history = zeros(first_n_steps + second_n_steps + 2, genlen);
mat_time = zeros(first_n_steps + second_n_steps + 2, 1);

%%

% run the simulation:
for i = 1:iterations
    
    sim_params = params(: , param_index);
    [~, time, nuc_s_hist, ~, ~, ~, ~] = ...
        run_simulation_from_genome(seq, 'poly_rate',0, ...
        'REB1_a_rate', 0.0001, 'REB1_e_rate', 0.0001, ...
        'ABF1_a_rate', 0.0001, 'ABF1_e_rate', 0.0001, ...
        'RAP1_a_rate', 0.0001, 'RAP1_e_rate', 0.0001, ...
        'TF_evic_intensity', sim_params(1), ...
        'RSC_evic_length', sim_params(2), 'RSC_slide_length', sim_params(2).*sim_params(3), ...
        'RSC_evic_intensity', sim_params(4), ...
        'RSC_slide_intensity', sim_params(4)*sim_params(5), 'slide_len', 3, 'gen_len', genlen, 'n_steps', first_n_steps);
    
    mat_time_1 = repmat(time, 1, genlen);
    wt_nuc_hist = nuc_s_hist .* mat_time_1;
    
    % second simulation with null model params:
    [~, time_2, nuc_s_hist_2, ~, ~, ~, ~] = ...
        run_null_model_simulation_from_genome(seq, 's0', nuc_s_hist(end,:), ...
        'slide_len', 3, 'gen_len', genlen, 'n_steps', second_n_steps);

    mat_time_2 = repmat(time_2, 1, genlen);
    sth1_nuc_hist = nuc_s_hist_2 .* mat_time_2;

    % combine the two state histories into one big state history:
    temp = [mat_time_1 ; mat_time_2];
    mat_time = mat_time + cumsum(temp(:, 1));
    nuc_history = nuc_history + (temp .* [wt_nuc_hist ; sth1_nuc_hist]);
end
mat_time = mat_time ./ iterations;

%%% TODO - use the mat_time and not just the simulation steps themselves

%%

% create a data matrix for the gene:
dat_matrix = zeros(7,genlen);
load('/cs/bd/Daniel/nuc_simulator/clustering/experiment_data/sth1_0m_centers.mat');
temp = create_gene_buffer(data(Gene_id,:),genlen);
dat_matrix(1, :) = temp;
load('/cs/bd/Daniel/nuc_simulator/clustering/experiment_data/sth1_10m_centers.mat');
temp = create_gene_buffer(data(Gene_id,:),genlen);
dat_matrix(2, :) = temp;
load('/cs/bd/Daniel/nuc_simulator/clustering/experiment_data/sth1_20m_centers.mat');
temp = create_gene_buffer(data(Gene_id,:),genlen);
dat_matrix(3, :) = temp;
load('/cs/bd/Daniel/nuc_simulator/clustering/experiment_data/sth1_30m_centers.mat');
temp = create_gene_buffer(data(Gene_id,:),genlen);
dat_matrix(4, :) = temp;
load('/cs/bd/Daniel/nuc_simulator/clustering/experiment_data/sth1_45m_centers.mat');
temp = create_gene_buffer(data(Gene_id,:),genlen);
dat_matrix(5, :) = temp;
load('/cs/bd/Daniel/nuc_simulator/clustering/experiment_data/sth1_60m_centers.mat');
temp = create_gene_buffer(data(Gene_id,:),genlen);
dat_matrix(6, :) = temp;
load('/cs/bd/Daniel/nuc_simulator/clustering/experiment_data/sth1_2h_centers.mat');
temp = create_gene_buffer(data(Gene_id,:),genlen);
dat_matrix(7, :) = temp;

%%

jump_len = 50;
window_size = 1000;
until_when = first_n_steps + second_n_steps - 5*window_size;

% create the Simulation Matrix:
sim_matrix = zeros(floor((until_when + window_size - first_n_steps) / jump_len), genlen);
for i = first_n_steps-window_size : jump_len : until_when
    sim_matrix(((i - first_n_steps+window_size) / jump_len) + 1,:) = sum(nuc_history(i:i+window_size , :));
end

%%

% get the best time expansion coefficient:
[timeExpansion, timeExpansions, costPerTimeExpansion] = getTimeExpansion(dat_matrix, sim_matrix);

%%

exp_time = [0, 10, 20, 30, 45, 60, 120];

% make a nuc_sum matrix for the different times:
nuc_sum_matrix = zeros(length(exp_time), length(NFR_pos));
data_matrix = zeros(length(exp_time), length(NFR_pos));
for i = 1:length(exp_time)
    data_matrix(i,:) = conv(dat_matrix(1, NFR_pos), gausswin(10)./sum(gausswin(10)), 'same');
    nuc_sum_matrix(i,:) = conv(sim_matrix(round(exp_time(1) * timeExpansion)+1, NFR_pos), gausswin(10)./sum(gausswin(10)), 'same');
end

% save everything for the review:
save([results_path 'dynamic_results_' num2str(Gene_id)] , ...
	'data_matrix', 'nuc_sum_matrix', 'timeExpansion', 'timeExpansions', 'costPerTimeExpansion');

end