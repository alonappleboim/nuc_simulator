function average_distance = get_average_peak_distance(centers, width)
% Given a vector of centers from a simulation and the width for the ksdensity smoothing,
% the function returns the average distance between the 5 nucs that are downstream of the
% +1 nuc. This is done by finding the first peak from the 1000 bp (the +1), and finding the
% differences accordingly (assuming each peak is correct with the findpeaks function - this
% means the given width has to be right).

% make the peaks and positions vectors:
centers = centers(1000:end);
smoothed = ksdensity(1:length(centers),1:length(centers),'weights',double(centers),'width',width);
[peaks,positions] = findpeaks(smoothed);

%{
% find the max peak (the +1 nucleosome):
peaks = transpose(peaks);
[max_peak,max_position] = max(peaks);
max_position = max_position(1);
%}

% find the average distance between the five nucleosomes:
%relevant_positions = positions(max_position:min((max_position+4),end));
relevant_positions = positions(1:min(4,end));
average_distance = (relevant_positions(end) - relevant_positions(1)) / (length(relevant_positions)-1);