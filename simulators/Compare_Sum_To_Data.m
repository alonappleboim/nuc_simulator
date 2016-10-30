function [ likelihood, plus1_dist, minus1_dist, peak_num_delta, plus_one_width_delta, minus_one_width_delta, peak_ratio_delta ] ...
    = Compare_Sum_To_Data( nuc_sum , data , NFR_pos, smooth)
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
%prob_dist = nuc_sum ./ sum(nuc_sum);
prob_dist = nuc_sum(NFR_pos) ./ sum(nuc_sum(NFR_pos));

%
data_likelihood = data(NFR_pos);

if (smooth == true)
    prob_dist = conv(prob_dist, gausswin(10)./sum(gausswin(10)), 'same');
    prob_dist = prob_dist ./ sum(prob_dist);
    data_conv = conv(data_likelihood, gausswin(10)./sum(gausswin(10)), 'same');
    data_likelihood = round(data_conv .* (sum(data_likelihood) / sum(data_conv)));
end

% the number of reads from the experimental data:
read_num = sum(data_likelihood);

% find the average amount of reads we expect in every bp:
lambda_vector = prob_dist .* read_num;

% get the log-chance that we get the data from the given distribution:
data_likelihood(data_likelihood == 0) = (min(data_likelihood(data_likelihood ~= 0)) / 1000); % actual zeros make problems in the log...
pois_vec = log_poisspdf(data_likelihood, lambda_vector);

%likelihood = sum(pois_vec(NFR_pos));
likelihood = sum(pois_vec);

% get the delta of the plus and minus 1, along with Peak Ratio and Width:
[plus1_dist, minus1_dist, peak_num_delta, plus_one_width_delta, minus_one_width_delta, peak_ratio_delta] ...
    = get_NFR_features(nuc_sum, data, NFR_pos);

end

