function [ timeExpansion, timeExpansions, costs ] = getTimeExpansion( data_matrix, sim_matrix )
%getTimeExpansion given the data and simulation matrices, the function
%returns the optimal time expansion constant that fits the simulation
%matrix onto the data matrix.

timeExpansions = [0.01 : 0.01 : 2];
costs = zeros(size(timeExpansions));

% check for every time expansion constant what the cost function is:
i=1;
for timeExp = timeExpansions
    costs(i) = dynamicsCostFunction(timeExp, data_matrix, sim_matrix);
    i = i+1;
end

% return the optimal time expansion constant:
[~, timeExpansionIndex] = min(costs);
timeExpansion = timeExpansions(timeExpansionIndex);

% plot the different time expansion coefficients costs:
%{
figure;
plot(timeExpansions(costs < 10000), costs(costs < 10000))
xlabel('Time Expansion Coefficient')
ylabel('Cost')
title('The Cost Of Different Time Coefficients')
%}

end


