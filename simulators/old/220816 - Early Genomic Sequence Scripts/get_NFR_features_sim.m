function [ mean_plus_one, mean_minus_one, mean_NFR_width ] = get_NFR_features_sim( s_hist )
%get_NFR_features_sim a function for getting the NFR features (the average distance
%between the +1 and -1 nucleosomes and the width) from a simulation.
%   given the state history and assuming the TSS is at position 1000, the
%   function calculates the mean distance between the +1 and -1 nucleosomes
%   and returns the features, using each time step in the simulation.

TSS = 1000;
initial_steps = 500;
max_distance = 400;
max_position = 1100;

distances = zeros(1,max_distance);
plus_one = zeros(1,max_position);
minus_one = zeros(1,max_position);

% extract the relevant locations for the two nucleosomes:
s_hist = s_hist(:, TSS-300 : TSS+100);

for i=initial_steps:size(s_hist,1)
    
    % find the positions of the nucs:
    nucs = find(s_hist(i,:));
    
    % calculate the distances:
    if (length(nucs) == 2)
        distances(nucs(2)-nucs(1)) = distances(nucs(2)-nucs(1)) + max(s_hist(i,:));
        plus_one((TSS-301) + nucs(2)) = plus_one((TSS-301) + nucs(2)) + max(s_hist(i,:));
        minus_one((TSS-301) + nucs(1)) = minus_one((TSS-301) + nucs(1)) + max(s_hist(i,:));

    end
    if (length(nucs) == 3) % if there are three, take the two closer to the TSS:
        distances(nucs(3)-nucs(2)) = distances(nucs(3)-nucs(2)) + max(s_hist(i,:));
        plus_one((TSS-301) + nucs(3)) = plus_one((TSS-301) + nucs(3)) + max(s_hist(i,:));
        minus_one((TSS-301) + nucs(2)) = minus_one((TSS-301) + nucs(2)) + max(s_hist(i,:));
    end 
end

% normalise the vectors:
distances = distances ./ sum(distances);
plus_one = plus_one ./ sum(plus_one);
minus_one = minus_one ./ sum(minus_one);

% calculate the mean distances and standard deviations:
mean_NFR_width = sum(distances .* [1:length(distances)]);
mean_plus_one = sum(plus_one .* [1:length(plus_one)]);
mean_minus_one = sum(minus_one .* [1:length(minus_one)]);

end

