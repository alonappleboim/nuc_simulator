load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\sequences_structure.mat')

genlen = 3500;
TSS = fix(genlen/2);
NFR_pos = [TSS-299 : TSS+150];

create_params_sth1;

for i=1:19

    gene_id = genes(i);
    
    load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_0m_centers.mat')
    
    try
        load(['C:\Users\Daniel\Documents\MATLAB\Friedman Lab\results\dynamic_results\results_0m_' num2str(gene_id) '.mat'])
    catch a
        continue
    end
    
    seq = sequences_structure(gene_id,:);
    wt =  create_gene_buffer(data(gene_id,:),genlen);
    
    [ PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites ] = ...
        Extract_Sites_From_Gene(seq, genlen);

    nuc_sum = nuc_sum_likelihood;
    best_sim = Compare_Sum_To_Data(nuc_sum, wt, NFR_pos, true);

    data = conv(wt(NFR_pos), gausswin(10)./sum(gausswin(10)), 'same');
    nuc_sum = conv(nuc_sum, gausswin(10)./sum(gausswin(10)), 'same');
    nuc_sum = nuc_sum(NFR_pos) .* sum(data) ./ sum(nuc_sum(NFR_pos));
    
    figure;
    subplot(2,1,1);
    plot(data, 'b')
    hold on
    plot(nuc_sum, 'r')
    plot(PolyA_Sites(NFR_pos), 'k')
    plot(PolyT_Sites(NFR_pos), 'm')
    legend('sth1_0m','simulation','PolyA', 'PolyT')
    xlabel('Position (TSS at 300)')
    ylabel('Nucleosome Intensity (0m)')
    title(['Gene ' num2str(gene_id) ' 0 minutes' char(10) num2str(best_ratio)])
     
    
    load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_6h_centers.mat')
    try
        load(['C:\Users\Daniel\Documents\MATLAB\Friedman Lab\results\dynamic_results\results_6h_' num2str(gene_id) '.mat'])
    catch a
        continue
    end
    
    wt =  create_gene_buffer(data(gene_id,:),genlen);

    [ PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites ] = ...
        Extract_Sites_From_Gene(seq, genlen);

    nuc_sum = nuc_sum_likelihood;
    best_sim = Compare_Sum_To_Data(nuc_sum, wt, NFR_pos, true);

    exp_data = conv(wt(NFR_pos), gausswin(10)./sum(gausswin(10)), 'same');
    nuc_sum = conv(nuc_sum, gausswin(10)./sum(gausswin(10)), 'same');
    nuc_sum = nuc_sum(NFR_pos) .* sum(exp_data) ./ sum(nuc_sum(NFR_pos));
        
    subplot(2,1,2);
    plot(exp_data, 'b')
    hold on
    plot(nuc_sum, 'r')
    plot(PolyA_Sites(NFR_pos), 'k')
    plot(PolyT_Sites(NFR_pos), 'm')
    legend('sth1_6h','simulation','PolyA', 'PolyT')
    xlabel('Position (TSS at 300)')
    ylabel('Nucleosome Intensity (6h)')
    title(['Gene ' num2str(gene_id) ' 6 hours' char(10) num2str(best_ratio)])
end


