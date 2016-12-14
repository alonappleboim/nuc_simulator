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

until_when = length(mat_time)-29000;
jump_len = 200;
% find the best index for every data point:
likelihoods = zeros(length(500:jump_len:until_when),10);
conv_likelihoods = zeros(length(500:jump_len:until_when),10);
for i = 1:10
    for j=500:jump_len:until_when
        nuc_sum = sum(nuc_history(j:j+9500 , :));
        likelihood = Compare_Sum_To_Data(nuc_sum, data_matrix(i,:), NFR_pos, true);
        likelihoods(((j-500)/jump_len)+1, i) = likelihood;
    end
end

% look at convolutions:
for i = 1:10
    conv_likelihoods(:,i) = conv(likelihoods(:,i), ones(3,1), 'same');
end

figure;
plot(likelihoods(:,1))
hold on
plot(likelihoods(:,2))
plot(likelihoods(:,3))
plot(likelihoods(:,4))
plot(likelihoods(:,5))
plot(likelihoods(:,6))
plot(likelihoods(:,7))
plot(likelihoods(:,8))
plot(likelihoods(:,9))
plot(likelihoods(:,10))
legend('1','2','3','4','5','6','7','8','9','10')

figure;
plot(conv_likelihoods(3:end,1))
hold on
plot(conv_likelihoods(3:end-3,2))
plot(conv_likelihoods(3:end-3,3))
plot(conv_likelihoods(3:end-3,4))
plot(conv_likelihoods(3:end-3,5))
plot(conv_likelihoods(3:end-3,6))
plot(conv_likelihoods(3:end-3,7))
plot(conv_likelihoods(3:end-3,8))
plot(conv_likelihoods(3:end-3,9))
plot(conv_likelihoods(3:end-3,10))
legend('1','2','3','4','5','6','7','8','9','10')

[~ , best_indices] = max(likelihoods);
[~ , conv_best_indices] = max(conv_likelihoods(3:end-3, :));
figure;
plot([0 10 20 30 45 60 120 180 270 360],best_indices, 'x')
hold on
plot([0 10 20 30 45 60 120 180 270 360],conv_best_indices, 'x')
title('maximum positions')
legend('reg','conv')

for i=1:10
    [~ , peaks] = findpeaks(likelihoods(:,i));%conv_likelihoods(:,i));
    if isempty(peaks)
        first_peaks_conv(i) = 1;
    else
        first_peaks_conv(i) = peaks(1);
    end
end
figure;
plot([0 10 20 30 45 60 120 180 270 360],first_peaks_conv, 'x')
title('first peaks from convolution')

%%

nuc_dynamics_movie(nuc_history,mat_time,'C:\Users\Daniel\Documents\MATLAB\nuc_simulator\visualization\movies\dynamics_test.mp4',Gene_id,wt_param,sth1_param);

%{

s_hist_coverage = ksdensity(1:length(centers_vector),1:length(centers_vector),'weights',double(centers_vector(1:end)),'width',5);

% plot wild type gene with PolyAT and the simulation (wt smoothed):
smoothed_wt = ksdensity([1:length(FRS2_wt)],[1:length(FRS2_wt)],'weights',double(FRS2_wt),'width',5);

[likelihood, plus1, minus1, peak_num_delta, plus1width, minus1width, ratio] ...
    = Compare_Sum_To_Data(centers_vector, FRS2_wt, NFR_pos, true)

figure;
plot(smoothed_wt(1:end-1),'b')
%plot(FRS2_wt(1:2500) ./ sum(FRS2_wt),'b')
hold on
plot(s_hist_coverage,'r')
%plot(centers_vector(1:2500) ./ sum(centers_vector(1:2500)),'r')
plot(PolyA_Sites .* mean(smoothed_wt),'k')
plot(PolyT_Sites .* mean(smoothed_wt),'m')
plot(REB1_Sites .* 4 .* mean(smoothed_wt), 'g')
plot(ABF1_Sites .* 4 .* mean(smoothed_wt), 'c')
plot(RAP1_Sites .* 4 .* mean(smoothed_wt), 'y')
legend('wild-type','simulation','PolyA (right)','PolyT (left)', 'REB1', 'ABF1', 'RAP1')
xlabel(['Position (TSS at ' num2str(fix(genlen/2)) ')'])
ylabel('Intensity')

[ a_rate, e_rate, r_rate, l_rate, ... 
            REB1_a, REB1_e, ABF1_a, ABF1_e, RAP1_a, RAP1_e] = ...
    generate_rates_from_sites( PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites);

figure
plot(conv(centers_vector(NFR_pos), gausswin(10), 'same'), 'b')
hold on
plot(conv(evics(NFR_pos), gausswin(10), 'same') .* sum(conv(centers_vector(NFR_pos), gausswin(10), 'same')) ./ sum(conv(evics(NFR_pos), gausswin(10), 'same')), 'r')
plot(e_rate(NFR_pos) .* sum(conv(centers_vector(NFR_pos), gausswin(10), 'same')) ./ sum(e_rate(NFR_pos)), 'k')
legend('nucleosome coverage', 'eviction count', 'eviction rates')

%{
% plot the NFR smoothed:
s_hist_coverage = ksdensity(1:length(FRS2_wt),1:length(FRS2_wt),'weights',double(centers_vector(1:length(FRS2_wt))),'width',5);
smoothed_wt = ksdensity([1:length(FRS2_wt)],[1:length(FRS2_wt)],'weights',double(FRS2_wt),'width',5);
figure;
plot(smoothed_wt(NFR_pos) .* sum(FRS2_wt(NFR_pos)),'g')
hold on
plot(s_hist_coverage(NFR_pos) .* sum(FRS2_wt(NFR_pos)), 'r')
%}

% plot the NFR normal:
centers_vector = conv(centers_vector,gausswin(5)./sum(gausswin(5)),'same');
centers_vector = centers_vector(NFR_pos) ./ sum(centers_vector(NFR_pos));
centers_vector = centers_vector .* sum(FRS2_wt(NFR_pos));
FRS2_wt = conv(FRS2_wt,gausswin(5)./sum(gausswin(5)),'same');
figure;
plot(FRS2_wt(NFR_pos),'b')
hold on
plot(centers_vector, 'r')
plot(evics(NFR_pos), 'g')
legend('wild-type', 'simulation', 'evictions')
%}