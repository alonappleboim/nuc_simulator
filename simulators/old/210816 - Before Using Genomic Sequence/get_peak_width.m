function [ width, positions ] = get_peak_width( s_hist, NFR_positions )
%get_peak_width a function for finding the width (standard deviation of the +1 nucleosome.
%   the function accepts the state history and NFR position, and finds the
%   positions of the +1 nucleosome, calculating the standard deviation of
%   its position throughout the simulation. It returns the standard
%   deviation, along with a vector of the +1 positions (around the relevant
%   region - the position numbers are not like the state history).

initial_steps = 500;

% extract the relevant locations for the +1 nucleosome:
s_hist = s_hist(:,mean(NFR_positions):mean(NFR_positions)+200);
positions = zeros(1,size(s_hist,2));

for i=1:size(s_hist,1)
    
    % find the position nuc:
    nuc_pos = find(s_hist(i,:));
    
    % calculate the distances:
    if (length(nuc_pos) == 1 && i > initial_steps)
        positions(nuc_pos) = positions(nuc_pos) + max(s_hist(i,:));
    end
        
end

% normalise the distance vector:
positions = positions ./ sum(positions);

% calculate the mean distances and standard deviations:
width = std(positions);

end

