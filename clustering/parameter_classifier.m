create_full_params;

very_good_genes = [9,16,19,21,63];
good_genes = [23,24,40,41,48,56,62,73,80,82,86,90,97];
middle_genes = [15,17,22,26,28,29,30,31,35,37,38,39,42,43,45,46,49,50,52,53,54,55,58,61,65,66,67,70,74,76,77,84,87,89,95];
bad_genes = [6,7,14,25,32,33,34,44,59,60,68,72,75,79,83,85,88,91,92,93,96,98,99];

very_good_evic = zeros(size(very_good_genes));
very_good_slide = zeros(size(very_good_genes));
very_good_length = zeros(size(very_good_genes));
good_evic = zeros(size(good_genes));
good_slide = zeros(size(good_genes));
good_length = zeros(size(good_genes));
middle_evic = zeros(size(middle_genes));
middle_slide = zeros(size(middle_genes));
middle_length = zeros(size(middle_genes));
bad_evic = zeros(size(bad_genes));
bad_slide = zeros(size(bad_genes));
bad_length = zeros(size(bad_genes));

for i = 1:length(very_good_genes)
    gene = very_good_genes(i);
    load(['C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\results\wt\results_' num2str(gene) '.mat'])
    parameters = params(: , best_sim_index);
    very_good_length(i) = parameters(4);
    very_good_evic(i) = parameters(5);
    very_good_slide(i) = parameters(6);
end

for i = 1:length(good_genes)
    gene = good_genes(i);
    load(['C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\results\wt\results_' num2str(gene) '.mat'])
    parameters = params(: , best_sim_index);
    good_length(i) = parameters(4);
    good_evic(i) = parameters(5);
    good_slide(i) = parameters(6);
end

for i = 1:length(good_genes)
    gene = good_genes(i);
    load(['C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\results\wt\results_' num2str(gene) '.mat'])
    parameters = params(: , best_sim_index);
    middle_length(i) = parameters(4);
    middle_evic(i) = parameters(5);
    middle_slide(i) = parameters(6);
end

for i = 1:length(bad_genes)
    gene = bad_genes(i);
    load(['C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\results\wt\results_' num2str(gene) '.mat'])
    parameters = params(: , best_sim_index);
    bad_length(i) = parameters(4);
    bad_evic(i) = parameters(5);
    bad_slide(i) = parameters(6);
end

figure
subplot(3,1,1)
plot(very_good_genes,very_good_evic,'^r')
title('Very Good Genes')
ylabel('Eviction Intensity')
xlabel('Gene ID')
subplot(3,1,2)
plot(very_good_genes,very_good_slide,'ob')
ylabel('Slide Intensity')
xlabel('Gene ID')
subplot(3,1,3)
plot(very_good_genes,very_good_length,'xk')
ylabel('RSC Length')
xlabel('Gene ID')

figure
subplot(3,1,1)
plot(good_genes,good_evic,'^r')
title('Good Genes')
ylabel('Eviction Intensity')
xlabel('Gene ID')
subplot(3,1,2)
plot(good_genes,good_slide,'ob')
ylabel('Slide Intensity')
xlabel('Gene ID')
subplot(3,1,3)
plot(good_genes,good_length,'xk')
ylabel('RSC Length')
xlabel('Gene ID')

figure
subplot(3,1,1)
plot(middle_genes,middle_evic,'^r')
title('Middle Genes')
ylabel('Eviction Intensity')
xlabel('Gene ID')
subplot(3,1,2)
plot(middle_genes,middle_slide,'ob')
ylabel('Slide Intensity')
xlabel('Gene ID')
subplot(3,1,3)
plot(middle_genes,middle_length,'xk')
ylabel('RSC Length')
xlabel('Gene ID')


figure
subplot(3,1,1)
plot(bad_genes,bad_evic,'^r')
title('Bad Genes')
ylabel('Eviction Intensity')
xlabel('Gene ID')
subplot(3,1,2)
plot(bad_genes,bad_slide,'ob')
ylabel('Slide Intensity')
xlabel('Gene ID')
subplot(3,1,3)
plot(bad_genes,bad_length,'xk')
ylabel('RSC Length')
xlabel('Gene ID')
