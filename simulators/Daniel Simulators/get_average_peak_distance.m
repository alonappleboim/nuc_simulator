function average_distance = get_average_peak_distance(smoothed)
% Given a smoothed vector of centers from a simulation (1500 bp, from TSS),
% the function returns the average distance between the 3 nucs that are downstream of the
% +1 nuc.

% make the peaks and positions vectors:
[peaks,positions] = findpeaks(smoothed, 'MinPeakHeight', 0.65*10^-3);

% find the average distance between the five nucleosomes:
%relevant_positions = positions(max_position:min((max_position+4),end));
relevant_positions = positions(1:min(3,end));
average_distance = (relevant_positions(end) - relevant_positions(1)) / (length(relevant_positions)-1);