%%% in this script we will understand what the correct NULL MODEL should
%%% be, given that the NULL MODEL only relies on AA/AT/TA/TT dinucleotids.

length = 4;
nuc = 100;

% we use the 6h data, assuming it has no RSC present and so the main
% effects are from eviction and nor from sliding:
load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\sequences_structure.mat')
load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_6h_centers.mat')

% for every gene that isn't null, get its xnucleotid vector:
genes = find(~isnan(data(:,1)));
dinucs = zeros(size(data));
for gene_id = genes'
    dinucs(gene_id, :) = get_xnucleotid_vector(sequences_structure(gene_id,:), length);
end

% for every position, sum up the number of xnucs that are around it:
buffer = ceil(nuc/2); %75;
for gene_id = genes'
    dinucs(gene_id, buffer:end-buffer) = conv(dinucs(gene_id, buffer:end-buffer), ones(1,nuc), 'same');
end

% normalize the data:
for gene_id = genes'
    data(gene_id,buffer:end-buffer) = data(gene_id,buffer:end-buffer) ./ sum(data(gene_id,buffer:end-buffer));
end

% for every position, get the intensity and add it to the vector of
% dinucleotid intensity:
dinuc_intensity = zeros(2,nuc);
for gene_id = genes'
    for index = buffer:max(size(data(1,:)))-buffer
        di_intensity = dinucs(gene_id, index);
        
        if (di_intensity ~= 0)
            dinuc_intensity(1,di_intensity) = dinuc_intensity(1,di_intensity) + data(gene_id, index);
            dinuc_intensity(2,di_intensity) = dinuc_intensity(2,di_intensity) + 1;
        end
    end
end

% normalize the dinuc vector according to the amount of times we saw each
% intensity:
di_intensity_vec = dinuc_intensity(1,:);
normalizer_vec = dinuc_intensity(2,:);
di_intensity_vec(normalizer_vec ~= 0) = di_intensity_vec(normalizer_vec ~= 0) ./ normalizer_vec(normalizer_vec ~= 0);

% plot the distribution of dinucs and their intensities (not normalized):
%{
figure;
plot(dinuc_intensity(1,:))
hold on
plot(dinuc_intensity(2,:) .* sum(dinuc_intensity(1,:)) ./ sum(dinuc_intensity(2,:)))
legend([num2str(length) '-nucleotid Intensity'],[num2str(length) '-nucleotid Count'])
xlabel(['# of ' num2str(length) '-nucleotids'])
ylabel('Intensity [a.u.]')
title(['Comparison between Intensity and Count of ' num2str(length) '-nucleotids'])
%}

% plot the normalized intensities of each number of dinucs:
figure;
plot(di_intensity_vec,'kx')
hold on
plot(dinuc_intensity(2,:) .* sum(di_intensity_vec) ./ sum(dinuc_intensity(2,:)), 'b')
legend('Normalized Intensity','# of Appearences Around a bp')
xlabel(['# of ' num2str(length) '-nucleotids'])
ylabel('Intensity [a.u.]')
title(['The Normalized Intensity of each ' num2str(length) '-nucleotid Count' char(10) num2str(nuc) ' bps around center'])

% find the factors by which to multiply the global eviction rate:
%rates_to_multiply = max(di_intensity_vec(1:50)) ./ di_intensity_vec(1:50);