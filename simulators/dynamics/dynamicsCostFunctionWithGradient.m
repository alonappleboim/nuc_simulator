function [ cost, gradient ] = dynamicsCostFunctionWithGradient( timeExpansion )
%dynamicsCostFunction calculate the sum of the likelihoods, given the time
%expansion coefficient.
%   we are assuming two vectors - data_matrix (containing the experimental
%   data) and sim_matrix (containing the simulation data).

cost = dynamicsCostFunction(timeExpansion);
gradient = dynamicsCostFunction(timeExpansion + 0.1) - dynamicsCostFunction(timeExpansion - 0.1);

end

