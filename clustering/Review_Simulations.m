load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\wt_centers.mat')
load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\sequences_structure.mat')

genlen = 3500;
TSS = fix(genlen/2);
NFR_pos = [TSS-299 : TSS+150];

create_full_params_RSC_ratio;

for gene_id=101:200

    try
        load(['C:\Users\Daniel\Documents\MATLAB\Friedman Lab\results\wt_RSC_ratio_161116\results_' num2str(gene_id) '.mat'])
    catch a
        continue
    end
    
    seq = sequences_structure(gene_id,:);
    wt =  wt_3h(gene_id,:);

    % make the wt data the right length:
    buffer = genlen - 2501;
    right_buffer = fix((buffer-500)/2);
    left_buffer = right_buffer + 500;
    if (right_buffer + left_buffer < buffer)
        left_buffer = left_buffer + 1;
    end
    wt = [zeros(1,left_buffer), wt, zeros(1,right_buffer)];

    [ PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites ] = ...
        Extract_Sites_From_Gene_new(seq, genlen, NFR_pos);

    nuc_sum = nuc_sum_likelihood; %nuc_sums(best_sim_index, :);
    best_sim = Compare_Sum_To_Data(nuc_sum, wt, NFR_pos, true); %%%%%

    data = conv(wt(NFR_pos), gausswin(10)./sum(gausswin(10)), 'same');
    nuc_sum = conv(nuc_sum, gausswin(10)./sum(gausswin(10)), 'same');
    nuc_sum = nuc_sum(NFR_pos) .* sum(data) ./ sum(nuc_sum(NFR_pos));
    
    [optimum_likelihood, ~,~,~,~,~,~] = ...
        Compare_Sum_To_Data(wt, wt, NFR_pos, true);
    [bad_likelihood,~,~,~,~,~,~] = ...
        Compare_Sum_To_Data(ones(size(wt)), wt, NFR_pos, true);
    ratio = (best_sim - bad_likelihood) / (optimum_likelihood - bad_likelihood);

    figure;
    %subplot(2,1,1);
    plot(data, 'b')
    hold on
    plot(nuc_sum, 'r')
    plot(PolyA_Sites(NFR_pos), 'k')
    plot(PolyT_Sites(NFR_pos), 'm')
    %plot(REB1_Sites(NFR_pos) .* 10,'k')
    %plot(ABF1_Sites(NFR_pos) .* 10,'k')
    %plot(RAP1_Sites(NFR_pos) .* 10,'k')
    legend('wild-type','simulation','Poly A', 'Poly T')
    xlabel('Position (TSS at 300)')
    ylabel('Nucleosome Intensity')
    %title(['Gene ' num2str(gene_id) char(10) 'Likelihood = ' num2str(features(best_sim_index))])
    title(['Gene ' num2str(gene_id) char(10) num2str(ratio)])
    %{
    textbox_string = [...
    'RSC Eviction Effect = ' num2str(params(5, best_sim_index)) char(10) ...
    'RSC Eviction Length = ' num2str(params(4, best_sim_index)) char(10) ...
    'RSC Slide Effect = ' num2str(params(6, best_sim_index)) char(10) ...
    'RSC Slide Length = ' num2str(params(4, best_sim_index) * 2) char(10) ...
    'TF Eviction Effect = ' num2str(params(3, best_sim_index)) ...
    ];
    annotation('textbox','String',textbox_string,'FitBoxToText','on','FontSize',10,'Position',[0.4,0.75,0.1,0.1])
    %}
    %%%%==================
    %{
    create_full_params;
    
    try
        load(['C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\results\wt_101016\results_' num2str(gene_id) '.mat'])
    catch a
        continue
    end
    
    seq = sequences_structure(gene_id,:);
    wt =  wt_3h(gene_id,:);

    % make the wt data the right length:
    buffer = genlen - 2501;
    right_buffer = fix((buffer-500)/2);
    left_buffer = right_buffer + 500;
    if (right_buffer + left_buffer < buffer)
        left_buffer = left_buffer + 1;
    end
    wt = [zeros(1,left_buffer), wt, zeros(1,right_buffer)];

    [ PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites ] = ...
        Extract_Sites_From_Gene(seq, genlen);

    nuc_sum = nuc_sums(best_sim_index, :);
    data = conv(wt(NFR_pos), gausswin(10)./sum(gausswin(10)), 'same');
    nuc_sum = conv(nuc_sum, gausswin(10)./sum(gausswin(10)), 'same');
    nuc_sum = nuc_sum(NFR_pos) .* sum(data) ./ sum(nuc_sum(NFR_pos));
    
    subplot(2,1,2);
    plot(data, 'b')
    hold on
    plot(nuc_sum, 'r')
    plot(REB1_Sites(NFR_pos) .* 10,'k')
    legend('wild-type','simulation','REB1')
    xlabel('Position (TSS at 300)')
    ylabel('Nucleosome Intensity')
    title(['Old Simulation'])
    textbox_string = [...
    'RSC Eviction Effect = ' num2str(params(5, best_sim_index)) char(10) ...
    'RSC Eviction Length = ' num2str(params(4, best_sim_index)) char(10) ...
    'RSC Slide Effect = ' num2str(params(6, best_sim_index)) char(10) ...
    'RSC Slide Length = ' num2str(params(4, best_sim_index) * 2) char(10) ...
    'TF Eviction Effect = ' num2str(params(3, best_sim_index)) ...
    ];
    annotation('textbox','String',textbox_string,'FitBoxToText','on','FontSize',10,'Position',[0.4,0.28,0.1,0.1])
    %}
end


%% make the histogram of each parameter to find global parameters

good_genes = [16,22,28,38,53,54,62,65,80];
good_genes_2 = [97,93,90,89,83,82,77,74,73,61,56,46,40,38,30,29,23,22,16,9,195,190,189,187,185,184,183,180,178,176,170,169,168,161,160,155,152,150,149,148,144,138,117];
RSC_evics = zeros(1,length(good_genes_2));
RSC_slides = zeros(1,length(good_genes_2));
RSC_len = zeros(1,length(good_genes_2));
TF_evic = zeros(1,length(good_genes_2));

%{
for i = 1:length(good_genes)
    create_full_params_271016;
    
    try
        load(['C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\results\wt_271016\results_' num2str(good_genes(i)) '.mat'])
    catch a
        continue
    end
    
    RSC_evics(i) = params(5, best_sim_index);
    RSC_slides(i) = params(6, best_sim_index);
    RSC_len(i) = params(4, best_sim_index);
    TF_evic(i) = params(3, best_sim_index);

end
%}
for i = 1:length(good_genes_2)
    
    try
        load(['C:\Users\Daniel\Documents\MATLAB\Friedman Lab\results\wt_RSC_ratio_161116\results_' num2str(good_genes_2(i)) '.mat'])
    catch a
        continue
    end
    
    RSC_evics(i) = params(5, best_sim_index);
    RSC_slides(i) = params(6, best_sim_index);
    RSC_len(i) = params(4, best_sim_index);
    TF_evic(i) = params(3, best_sim_index);

end

RSC_evics = RSC_evics(RSC_evics>0);
RSC_slides = RSC_slides(RSC_slides>0);
RSC_len = RSC_len(RSC_len>0);
TF_evic = TF_evic(TF_evic>0);

%%
figure;
subplot(3,1,1)
hist(RSC_evics(RSC_evics>0))
xlabel('RSC eviction Parameter')
ylabel('simulations')
title('Histograms of each parameter')
subplot(3,1,2)
hist(RSC_slides(RSC_slides>0))
xlabel('RSC slide Parameter')
ylabel('simulations')
subplot(3,1,3)
hist(RSC_len(RSC_len>0))
xlabel('RSC length Parameter')
ylabel('simulations')

%%
figure;
ndhist(RSC_evics(RSC_evics>0), RSC_slides(RSC_slides>0));
colorbar;
xlabel('RSC evics')
ylabel('RSC slides')

figure;
ndhist(RSC_len(RSC_len>0), RSC_slides(RSC_slides>0));
colorbar;
xlabel('RSC length')
ylabel('RSC slides')

figure;
ndhist(RSC_evics(RSC_evics>0), RSC_len(RSC_len>0));
colorbar;
xlabel('RSC evics')
ylabel('RSC length')
%%

