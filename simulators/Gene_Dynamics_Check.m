load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\sequences_structure.mat')

genes = [86,169,185,144,168,3903,4284,3051,4694,4513,6210,3227,1841,2327,6178,3655,3656,905,5776];
TSS = 1000;
NFR_pos = [TSS-299 : TSS+100];

for gene_id = genes
    load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_0m_centers.mat')
    data = data(gene_id,:);
    data = conv(data(NFR_pos), gausswin(10)./sum(gausswin(10)), 'same');

    figure('units','normalized','position',[0 0 1 1]);
    
    subplot(5,1,1);
    plot(data, 'b')
    legend('0m')
    title(num2str(gene_id))
    
    load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_20m_centers.mat')
    data = data(gene_id,:);
    data = conv(data(NFR_pos), gausswin(10)./sum(gausswin(10)), 'same');
    subplot(5,1,2);
    plot(data, 'b')
    legend('20m')

    load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_60m_centers.mat')
    data = data(gene_id,:);
    data = conv(data(NFR_pos), gausswin(10)./sum(gausswin(10)), 'same');
    subplot(5,1,3);
    plot(data, 'b')
    legend('60m')

    load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_4_5h_centers.mat')
    data = data(gene_id,:);
    data = conv(data(NFR_pos), gausswin(10)./sum(gausswin(10)), 'same');
    subplot(5,1,4);
    plot(data, 'b')
    legend('4.5h')

    load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_6h_centers.mat')
    data = data(gene_id,:);
    data = conv(data(NFR_pos), gausswin(10)./sum(gausswin(10)), 'same');
    subplot(5,1,5);
    plot(data, 'b')
    legend('6h')

    
end