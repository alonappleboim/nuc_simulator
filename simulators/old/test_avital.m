load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\wt_centers.mat')
%load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\sequences_structure.mat')

NFR = [800:1100];
genes = find(~(sequences_structure(:,1) == ' '));
As = zeros(1,length(genes));
Ts = zeros(1,length(genes));
j=0;

for i = genes'
    j = j+1;
    genome = sequences_structure(i,NFR);
    PolyA_Sites = zeros(1,length(NFR));
    PolyT_Sites = zeros(1,length(NFR));
    
    PolyA_Sites(genome == 'A') = 1;
    PolyA_Sites = conv(PolyA_Sites, ones(1,7), 'same');
    PolyA_Sites(PolyA_Sites < 5) = 0;
    PolyA_Sites(PolyA_Sites > 4) = 1;
    PolyA_Sites = conv(PolyA_Sites,ones(1,5),'same');
    PolyA_Sites(PolyA_Sites > 0) = 1;

    PolyT_Sites(genome == 'T') = 1;
    PolyT_Sites = conv(PolyT_Sites, ones(1,7), 'same');
    PolyT_Sites(PolyT_Sites < 5) = 0;
    PolyT_Sites(PolyT_Sites > 4) = 1;
    PolyT_Sites = conv(PolyT_Sites,ones(1,5),'same');
    PolyT_Sites(PolyT_Sites > 0) = 1;
    
    As(j) = sum(PolyA_Sites);
    Ts(j) = sum(PolyT_Sites);
end

low_poly_genes = find(As+Ts < quantile(As+Ts,0.05));
high_poly_genes = find(As+Ts > quantile(As+Ts,0.95));
high_A_low_T_genes = find(As > quantile(As,0.8) & Ts < quantile(Ts,0.2));
low_A_high_T_genes = find(Ts > quantile(Ts,0.8) & As < quantile(As,0.2));
indices = [low_poly_genes, high_poly_genes, high_A_low_T_genes, low_A_high_T_genes];
%genes = genes(indices);

create_params_avital;

%%

genlen = 3500;
TSS = fix(genlen/2);
NFR_pos = [TSS-299 : TSS+150];
for gene_id = genes(high_poly_genes(226:end))'
    try
        load(['C:\Users\Daniel\Documents\MATLAB\Friedman Lab\results\wt_avital\results_' num2str(gene_id) '.mat'])
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
end

%%
for i = high_poly_genes'
    
end

for i = high_A_low_T_genes'
    
end

for i = low_A_high_T_genes'
    
end
