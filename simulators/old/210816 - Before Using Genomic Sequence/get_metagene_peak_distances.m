function [distance_1_2, distance_2_3] = get_metagene_peak_distances(smoothed)
% Given a smoothed vector of centers (1500 bp, from TSS),
% the function returns the average distance between the 3 nucs that are downstream of the
% +1 nuc.

% make the peaks and positions vectors:
[peaks,positions] = findpeaks(smoothed);
[minima,minima_positions] = findpeaks(smoothed.*(-1));

% find the first peak:
for i=1:3
	if (positions(i) > minima_positions(1))
		peak1 = i;
        break
	end
end

% find the average distance between the five nucleosomes:
distance_1_2 = positions(peak1+1) - positions(peak1);
distance_2_3 = positions(peak1+2) - positions(peak1+1);