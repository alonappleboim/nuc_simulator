function [ mean_dist, std_dist ] = get_NFR_width( s_hist, NFR_positions )
%get_NFR_width a function for getting the NFR width (the average distance
%between the +1 and -1 nucleosomes).
%   given the state history and the NFR positions of the simulation, the
%   function calculates the mean distance between the +1 and -1 nucleosomes
%   and returns the standard deviation as well.

initial_steps = 500;
max_distance = 400;

distances = zeros(1,max_distance);

% extract the relevant locations for the two nucleosomes:
s_hist = s_hist(:,NFR_positions(1)-100 : NFR_positions(end)+100);

for i=1:size(s_hist,1)
    
    % find the positions of the nucs:
    nucs = find(s_hist(i,:));
    
    % calculate the distances:
    if (length(nucs) == 2 && i > initial_steps)
        distances(nucs(2)-nucs(1)) = distances(nucs(2)-nucs(1)) + max(s_hist(i,:));
    end
        
end

% normalise the distance vector:
distances = distances ./ sum(distances);

% calculate the mean distances and standard deviations:
mean_dist = sum(distances .* [1:length(distances)]);
std_dist = std(distances);

end

