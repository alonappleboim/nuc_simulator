load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\wt_centers.mat')
load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\sequences_structure.mat')

genlen = 3500;
TSS = fix(genlen/2);
NFR_pos = [TSS-299 : TSS+150];

create_params_new;

for gene_id=1:100

    try
        load(['C:\Users\Daniel\Documents\MATLAB\Friedman Lab\results\wt_RSC_length_141116\results_' num2str(gene_id) '.mat'])
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

create_params_avital;

lowAhighT = [1895,1845,1829,1691,1645,1592,1099,1097,1091,1089,1030,951,841,831,814,671,568,559,509,494,407,271,184,38,3100,2878,2594,2476,2351,2298,2031,1949,4105,4097,4064,3910,3877,3849,3788,6501,6340,6283,6232,6190,5805,5744,5712,5656,5624,5497,5491,5431,5237,4968,4931,4779,4666];
low_poly_genes = [1348,1298,1276,1193,1135,886,508,417,173,3622,3575,3449,3317,3277,3132,2976,2954,2739,2659,2542,2425,2348,2330,2096,1870,5330,4289,4117,6621,6587,6545,6543,6387,6110,6109,5684,5506,5485];
highAlowT = [1791,1633,1596,1486,1434,1291,1234,1153,891,879,785,731,630,495,343,259,190,148,73,22,3152,3058,3019,2987,2959,2949,2899,2871,2595,2569,2060,1941,1923,1907,1830,1813,4727,4499,4364,4245,4155,3744,3703,3681,3626,3590,3588,3573,3493,3448,3331,3288,3196,6624,6608,6598,6566,6434,6312,6106,5743,5609,5574,5244,5132,5062,4925];
high_poly_genes = [1920,1822,1576,1266,1233,974,969,967,860,804,762,607,597,586,195,140,40,9,3758,3636,3448,3241,3237,3156,3027,2937,2880,2740,2666,2632,2613,2515,2508,2480,2338,2294,2196,2065,2042,2039,5544,5543,5325,4939,4771,4745,4731,4631,4606,4500,4387,4246,4168,4152,3980,3964,3930,3917,3916,6556,6472,6243,6171,5835,5763,5747,5665,5571];
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
for i = 1:length(high_poly_genes)
    
    try
        load(['C:\Users\Daniel\Documents\MATLAB\Friedman Lab\results\wt_avital\results_' num2str(high_poly_genes(i)) '.mat'])
    catch a
        continue
    end
    
    RSC_evics4(i) = params(4, best_sim_index);
    RSC_slides4(i) = params(5, best_sim_index);
    RSC_len4(i) = params(2, best_sim_index);
    TF_evic4(i) = params(1, best_sim_index);

end

RSC_evics = RSC_evics(RSC_evics>0);
RSC_slides = RSC_slides(RSC_slides>0);
RSC_len = RSC_len(RSC_len>0);
TF_evic = TF_evic(TF_evic>0);

%%
figure;
subplot(4,1,1)
hist(RSC_len1,30)
%xlabel('RSC eviction Parameter')
ylabel('Low A Low T')
title('Histograms of groups of genes and their RSC Eviction Parameter')
subplot(4,1,2)
hist(RSC_len2,30)
%xlabel('RSC slide Parameter')
ylabel('Low A High T')
subplot(4,1,3)
hist(RSC_len3(RSC_len3>0),30)
%xlabel('RSC length Parameter')
ylabel('High A Low T')
subplot(4,1,4)
hist(RSC_len4,30)
xlabel('RSC Eviction Parameter')
ylabel('High A High T')


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

