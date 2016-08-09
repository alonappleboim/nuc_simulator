function [width_1, width_2] = get_peak_width(smoothed)
% Given a smoothed vector of centers from a simulation (1500 bp, from NFR start),
% the function returns the width of the +1 and +2 nucs (the distance between the two minima)

% make the peaks and positions vectors:
[peaks,peak_positions] = findpeaks(smoothed);
[minima,minima_positions] = findpeaks(smoothed.*(-1));

% find the widths:
width_1 = minima_positions(2) - minima_positions(1);
width_2 = minima_positions(3) - minima_positions(2);