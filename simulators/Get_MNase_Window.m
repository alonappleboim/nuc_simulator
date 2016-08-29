function [ window ] = Get_MNase_Window( window_size )
%Get_MNase_Window get the normalized window for the convolution in the end of the
%simulation (simulating the MNase digestion). Use window_size=1 for no
%effect...

load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\quantiles_feature_detailed_3.mat')
load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\wt_centers.mat')

% is 7 the right part of the features to take? (the center of 100% of the
% nucs)
centers = center_features_p1(gene_ids, 7);
wt_3h_p1Alligned = wt_3h(gene_ids, [centers-500 : centers+500]);
alligned_sum = sum(wt_3h_p1Alligned);

window = alligned_sum([500-fix(window_size/2) : 500+fix(window_size/2)]) ./ sum(alligned_sum([500-fix(window_size/2) : 500+fix(window_size/2)]));