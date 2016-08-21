function [ mean_dist_1_2, mean_dist_2_3, std_1_2, std_2_3] = get_peak_distances( s_hist, NFR_positions )
%get_peak_distance a function for finding the average distance between the +1 and
%+2 nucleosomes, and between the +2 and +3 nucleosomes (and also the standard deviation).
%   The function accepts the state history of the simulation, along with
%   the NFR positions. It then looks in the +1+2 area for 2 nucs in every
%   step (and the same for +2+3), calculates the distance and then sums up
%   the whole thing.

initial_steps = 500;
max_distance = 300;

distances_1_2 = zeros(1,max_distance);
distances_2_3 = zeros(1,max_distance);

% extract the relevant locations for +1+2 and +2+3 nucleosomes:
s_hist_1_2 = s_hist(:,mean(NFR_positions):mean(NFR_positions)+max_distance);
s_hist_2_3 = s_hist(:,NFR_positions(end)+50:NFR_positions(end)+50+max_distance);

for i=1:size(s_hist,1)
    
    % find the positions of the nucs:
    nucs_1_2 = find(s_hist_1_2(i,:));
    nucs_2_3 = find(s_hist_2_3(i,:));
    
    % calculate the distances:
    if (length(nucs_1_2) == 2 && i > initial_steps)
        distances_1_2(nucs_1_2(2)-nucs_1_2(1)) = distances_1_2(nucs_1_2(2)-nucs_1_2(1)) + max(s_hist_1_2(i,:));
    end
    if (length(nucs_2_3) == 2 && i > initial_steps)
        distances_2_3(nucs_2_3(2)-nucs_2_3(1)) = distances_2_3(nucs_2_3(2)-nucs_2_3(1)) + max(s_hist_2_3(i,:));
    end
        
end

% normalise the distance vector:
distances_1_2 = distances_1_2 ./ sum(distances_1_2);
distances_2_3 = distances_2_3 ./ sum(distances_2_3);

% calculate the mean distances and standard deviations:
mean_dist_1_2 = sum(distances_1_2 .* [1:length(distances_1_2)]);
mean_dist_2_3 = sum(distances_2_3 .* [1:length(distances_2_3)]);
std_1_2 = std(distances_1_2);
std_2_3 = std(distances_2_3);

end

