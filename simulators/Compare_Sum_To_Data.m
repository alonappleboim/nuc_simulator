function [ likelihood, plus1_dist, minus1_dist ] = Compare_Sum_To_Data( nuc_sum , data , NFR_pos, smooth)
%Compare_Sum_To_Data a function for comparison between a simulation and the
%experimental data - this is what we want to maximize when we are fitting
%parameters.
%   The parameters are the nucleosome state history of the simulation, 
%   the experimental data (the reads) and the NFR positions (701:1100), for  
%   example - the two vectors need to be of the
%   same region of the gene (the length should be 2501 bp). Also, the
%   functions accepts a boolean argument (smooth) - True for smoothing the
%   vectors before the comparison and false for keeping them the same.
%   The function uses the s_hist as the probability distribution of the
%   nucleosomes, and then checks what is the chance to get the experimental
%   data, given the probability distribution.
%   Assumptions of the function:
%       - The MNase cuts exactly 147 bps (the exact nucleosomes)
%       - The distribution of the experimental data, given the probability
%         distribution from the simulation, is Poissonal.

% the normalized probability distribution from the simulation:
prob_dist = nuc_sum ./ sum(nuc_sum);

if (smooth == true)
    prob_dist = conv(prob_dist, gausswin(10)./sum(gausswin(10)), 'same');
    prob_dist = prob_dist ./ sum(prob_dist);
    data = round(conv(data, gausswin(10)./sum(gausswin(10)), 'same'));
end

% the number of reads from the experimental data:
read_num = sum(data);

% find the average amount of reads we expect in every bp:
lambda_vector = prob_dist .* read_num;

% get the log-chance that we get the data from the given distribution:
data(data == 0) = 0.0001; % actual zeros make problems in the log...
pois_vec = log_poisspdf(data, lambda_vector);

likelihood = sum(pois_vec(NFR_pos));

% stop using get_NFR_features_data until we fix the problems (temp is
% sometimes of length 0)
%{
[plus1_data , minus1_data, temp] = get_NFR_features_data(data);
[plus1_sim, minus1_sim, temp] = get_NFR_features_data(nuc_sum);

plus1_dist = abs(plus1_sim - plus1_data);
minus1_dist = abs(minus1_sim - minus1_data);
%}

plus1_dist = 100;
minus1_dist = 100;

end

