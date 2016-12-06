load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\sequences_structure.mat')

genlen = 3500;
TSS = fix(genlen/2);
NFR_pos = [TSS-299 : TSS+150];

genes = [86,169,185,144,168,3903,4284,3051,4694,4513,6210,3227,1841,2327,6178,3655,3656,905,5776];
%genes = randi([1 6648],1,20);
for gene_id=genes
    
    load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_0m_centers.mat')
    
    seq = sequences_structure(gene_id,:);
    if (seq == ' ')
        continue
    end
    wt =  create_gene_buffer(data(gene_id,:),genlen);
    
    [ PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites ] = ...
        Extract_Sites_From_Gene(seq, genlen);

    data = conv(wt(NFR_pos), gausswin(10)./sum(gausswin(10)), 'same');
    
    figure('units','normalized','position',[0 0 1 1]);
    subplot(2,1,1);
    plot(data, 'b')
    hold on
    plot(PolyA_Sites(NFR_pos), 'k')
    plot(PolyT_Sites(NFR_pos), 'm')
    legend('sth1_0m','PolyA', 'PolyT')
    xlabel('Position (TSS at 300)')
    ylabel('Nucleosome Intensity (0m)')
    title(['Gene ' num2str(gene_id) ' 0 minutes'])
     
    
    load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_6h_centers.mat')
    
    wt =  create_gene_buffer(data(gene_id,:),genlen);

    [ PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites ] = ...
        Extract_Sites_From_Gene(seq, genlen);

    exp_data = conv(wt(NFR_pos), gausswin(10)./sum(gausswin(10)), 'same');
        
    subplot(2,1,2);
    plot(exp_data, 'b')
    hold on
    plot(PolyA_Sites(NFR_pos), 'k')
    plot(PolyT_Sites(NFR_pos), 'm')
    legend('sth1_6h','PolyA', 'PolyT')
    xlabel('Position (TSS at 300)')
    ylabel('Nucleosome Intensity (6h)')
    title(['Gene ' num2str(gene_id) ' 6 hours'])
end


