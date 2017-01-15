%%% this was used before I make the cost function that finds the best time
%%% expasion coefficiennt

%%

iterations = 100;
first_n_steps = 10000;
second_n_steps = 30000;

gene_index = 7;

load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\sequences_structure.mat')

genes = [3903, 168, 169, 3655, 4694, 86, 4513];
wt_param_indices = [251, 314, 329, 344, 321, 38, 319];
sth1_6h_param_indices = [126,125,123,189,205,305,302];

Gene_id = genes(gene_index);
wt_param = wt_param_indices(gene_index);
sth1_param = sth1_6h_param_indices(gene_index);
create_params_sth1;

genlen = 3500;
TSS = round(genlen/2);
NFR_pos = [TSS-299 : TSS+150];

% Starting Variables:
seq = sequences_structure(Gene_id,:);

[ PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites ] = ...
    Extract_Sites_From_Gene(seq, genlen);

nuc_history = zeros(first_n_steps + second_n_steps + 2, genlen);
mat_time = zeros(first_n_steps + second_n_steps + 2, 1);
% run the simulation:
for i = 1:iterations
    
    % first simulation with wt params:
    sim_params = params(: , wt_param);
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
    
    % second simulation with sth1 params and s0 as the last state from the
    % wt simulation:
    sim_params = params(: , sth1_param);
    [~, time_2, nuc_s_hist_2, ~, ~, ~, ~] = ...
        run_simulation_from_genome(seq, 'poly_rate',0, 's0', nuc_s_hist(end,:), ...
        'REB1_a_rate', 0.0001, 'REB1_e_rate', 0.0001, ...
        'ABF1_a_rate', 0.0001, 'ABF1_e_rate', 0.0001, ...
        'RAP1_a_rate', 0.0001, 'RAP1_e_rate', 0.0001, ...
        'TF_evic_intensity', sim_params(1), ...
        'RSC_evic_length', sim_params(2), 'RSC_slide_length', sim_params(2).*sim_params(3), ...
        'RSC_evic_intensity', sim_params(4), ...
        'RSC_slide_intensity', sim_params(4)*sim_params(5), 'slide_len', 3, 'gen_len', genlen, 'n_steps', second_n_steps);

    mat_time_2 = repmat(time_2, 1, genlen);
    sth1_nuc_hist = nuc_s_hist_2 .* mat_time_2;

    % combine the two state histories into one big state history:
    temp = [mat_time_1 ; mat_time_2];
    mat_time = mat_time + cumsum(temp(:, 1));
    nuc_history = nuc_history + (temp .* [wt_nuc_hist ; sth1_nuc_hist]);
end
mat_time = mat_time ./ iterations;

%%

% create a data matrix for the gene:
data_matrix = zeros(10,genlen);
load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_0m_centers.mat')
temp = create_gene_buffer(data(Gene_id,:),genlen);
data_matrix(1, :) = temp;
load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_10m_centers.mat')
temp = create_gene_buffer(data(Gene_id,:),genlen);
data_matrix(2, :) = temp;
load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_20m_centers.mat')
temp = create_gene_buffer(data(Gene_id,:),genlen);
data_matrix(3, :) = temp;
load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_30m_centers.mat')
temp = create_gene_buffer(data(Gene_id,:),genlen);
data_matrix(4, :) = temp;
load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_45m_centers.mat')
temp = create_gene_buffer(data(Gene_id,:),genlen);
data_matrix(5, :) = temp;
load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_60m_centers.mat')
temp = create_gene_buffer(data(Gene_id,:),genlen);
data_matrix(6, :) = temp;
load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_2h_centers.mat')
temp = create_gene_buffer(data(Gene_id,:),genlen);
data_matrix(7, :) = temp;
load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_3h_centers.mat')
temp = create_gene_buffer(data(Gene_id,:),genlen);
data_matrix(8, :) = temp;
load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_4_5h_centers.mat')
temp = create_gene_buffer(data(Gene_id,:),genlen);
data_matrix(9, :) = temp;
load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_6h_centers.mat')
temp = create_gene_buffer(data(Gene_id,:),genlen);
data_matrix(10, :) = temp;

%%

jump_len = 300;
window_size = 9000;
until_when = 10000 + 0.7*window_size;
% find the best index for every data point:
likelihoods = zeros(length(10000-window_size:jump_len:until_when),10);
conv_likelihoods = zeros(size(likelihoods));
for i = 1:10
    for j=10000-window_size:jump_len:until_when
        nuc_sum = sum(nuc_history(j:j+window_size , :));
        likelihood = Compare_Sum_To_Data(nuc_sum, data_matrix(i,:), NFR_pos, true);
        likelihoods(((j-10000+window_size)/jump_len)+1, i) = likelihood;
    end
end

% look at convolutions:
for i = 1:10
    conv_likelihoods(:,i) = conv(likelihoods(:,i), ones(3,1), 'same');
end

figure;
plot(mat_time(jump_len.*(cumsum(ones(size(likelihoods(:,1))))-1)+1),likelihoods(:,1),'color',[1,0,1],'LineWidth',1)
hold on
plot(mat_time(jump_len.*(cumsum(ones(size(likelihoods(:,1))))-1)+1),likelihoods(:,2),'color',[0.9,0,0.9],'LineWidth',1)
plot(mat_time(jump_len.*(cumsum(ones(size(likelihoods(:,1))))-1)+1),likelihoods(:,3),'color',[0.8,0,0.8],'LineWidth',1)
plot(mat_time(jump_len.*(cumsum(ones(size(likelihoods(:,1))))-1)+1),likelihoods(:,4),'color',[0.7,0,0.7],'LineWidth',1)
plot(mat_time(jump_len.*(cumsum(ones(size(likelihoods(:,1))))-1)+1),likelihoods(:,5),'color',[0.6,0,0.6],'LineWidth',1)
plot(mat_time(jump_len.*(cumsum(ones(size(likelihoods(:,1))))-1)+1),likelihoods(:,6),'color',[0.5,0,0.5],'LineWidth',1)
plot(mat_time(jump_len.*(cumsum(ones(size(likelihoods(:,1))))-1)+1),likelihoods(:,7),'color',[0.4,0,0.4],'LineWidth',1)
plot(mat_time(jump_len.*(cumsum(ones(size(likelihoods(:,1))))-1)+1),likelihoods(:,8),'color',[0.3,0,0.3],'LineWidth',1)
plot(mat_time(jump_len.*(cumsum(ones(size(likelihoods(:,1))))-1)+1),likelihoods(:,9),'color',[0.2,0,0.2],'LineWidth',1)
plot(mat_time(jump_len.*(cumsum(ones(size(likelihoods(:,1))))-1)+1),likelihoods(:,10),'color',[0.1,0,0.1],'LineWidth',1)
legend('0m','10m','20m','30m','45m','60m','2h','3h','4.5h','6h')
title(['Dynamic Likelihoods Of Simulation Compared To Experimental Datapoints' char(10) 'Gene ' num2str(Gene_id)])
ylabel('Likelihood Score')
xlabel('Simulation Time')

[~ , best_indices] = max(likelihoods);
figure;
plot([0 10 20 30 45 60 120 180 270 360],mat_time(best_indices.*jump_len), 'x')
title(['Optimal Likelihood Simulation Time For Each Experimental Datapoint' char(10) 'Gene ' num2str(Gene_id)])
xlabel('Experiment Time [minutes]')
ylabel('Simulation Time')

%%

nuc_dynamics_movie(nuc_history(),mat_time,['C:\Users\Daniel\Documents\MATLAB\nuc_simulator\visualization\movies\dynamics_test_gene' num2str(Gene_id) '.mp4'],Gene_id,wt_param,sth1_param);
