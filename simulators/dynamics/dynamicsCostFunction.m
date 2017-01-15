function [ cost ] = dynamicsCostFunction( timeExpansion, data_matrix, sim_matrix )
%dynamicsCostFunction calculate the sum of the likelihoods, given the time
%expansion coefficient.
%   we are assuming two vectors - data_matrix (containing the experimental
%   data) and sim_matrix (containing the simulation data).

genlen = 3500;
TSS = round(genlen/2);
NFR_pos = [TSS-299 : TSS+150];

exp_time = [0, 10, 20, 30, 45, 60, 120];
sim_indices = round(timeExpansion .* exp_time);

% deal with time expansions that are too big:
if (1+sim_indices(end) > length(sim_matrix(:,1)))
    cost = 100000;
else
    cost = 0;
    for i = 1 : length(sim_indices)
        likelihood = Compare_Sum_To_Data(sim_matrix(1+sim_indices(i),:), data_matrix(i,:), NFR_pos, true);
        cost = cost - likelihood;
    end
end

end

