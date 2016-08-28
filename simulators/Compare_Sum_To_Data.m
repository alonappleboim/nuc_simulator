function [ feature ] = Compare_Sum_To_Data( nuc_sum , data )
%Compare_Sum_To_Data a function for comparison between a simulation and the
%experimental data - this is what we want to maximize when we are fitting
%parameters.
%   The parameters are the nucleosome state history of the simulation and
%   the experimental data (the reads) - the two vectors need to be of the
%   same region of the gene (the length doesn't have to be 2501 bp).
%   The function uses the s_hist as the probability distribution of the
%   nucleosomes, and then checks what is the chance to get the experimental
%   data, given the probability distribution.
%   Assumptions of the function:
%       - The MNase cuts exactly 147 bps (the exact nucleosomes)
%       - The distribution of the experimental data, given the probability
%         distribution from the simulation, is Poissonal.

% the normalized probability distribution from the simulation:
prob_dist = nuc_sum ./ sum(nuc_sum);
%prob_dist = ksdensity(1:length(prob_dist),1:length(prob_dist),'weights',double(prob_dist(1:end)),'width',5);

% the number of reads from the experimental data:
read_num = sum(data);

% find the average amount of reads we expect in every bp:
lambda_vector = prob_dist .* read_num;
%lambda_vector = ones(1,length(data)) .* read_num ./ length(data);

% get the log-chance that we get the data from the given distribution:
pois_vec = log_poisspdf(data, lambda_vector);

feature = sum(pois_vec);

end

