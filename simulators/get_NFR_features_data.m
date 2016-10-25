function [ plus_one_delta, minus_one_delta, peak_num_delta ] = get_NFR_features_data( sim_centers, data_centers, NFR_pos )
%get_NFR_features_data a function for extracting NFR features from the experimental
%data.
%   given the nuc sum vector from the data and from the simulation, along
%   with the NFR positions (a vector like [800:1100]), the function returns
%   the delta between the plus and minus one nucs of the simulation and the
%   data, and the delta of the number of peaks. We assume that the
%   right-most peak will be the +1 (which means the NFR position should be
%   appropriate...
%   In case there were less than two peaks in the nuc sum vector from the
%   data, the function will return NaN for the plus and minus one deltas.

% find the peaks of the data and the simulation:
sim_smooth = conv(sim_centers, gausswin(100), 'same');
sim_smooth = conv(sim_smooth, gausswin(20), 'same');
[sim_peaks,sim_positions] = findpeaks(sim_smooth, 'MinPeakHeight', sum(sim_smooth)./(1.5*length(sim_smooth)));
data_smooth = conv(data_centers, gausswin(100), 'same');
data_smooth = conv(data_smooth, gausswin(20), 'same');
[data_peaks,data_positions] = findpeaks(data_smooth, 'MinPeakHeight', sum(data_smooth)./(1.5*length(data_smooth)));

% keep just the NFR peaks:
sim_positions = sim_positions((sim_positions > NFR_pos(1)) & (sim_positions < NFR_pos(end)));
sim_peaks = sim_peaks((sim_positions > NFR_pos(1)) & (sim_positions < NFR_pos(end)));
data_positions = data_positions((data_positions > NFR_pos(1)) & (data_positions < NFR_pos(end)));
data_peaks = data_peaks((data_positions > NFR_pos(1)) & (data_positions < NFR_pos(end)));

% for every data peak, find the distance to the closest simulation peak:
deltas = zeros(size(data_positions));
for i = 1:length(data_positions)
    deltas(i) = min(abs(data_positions(i) - sim_positions));
end

% find the plus and minus one deltas (assumin the right-most peak is the +1)
if (length(deltas) > 1)
    plus_one_delta = deltas(end);
    minus_one_delta = deltas(end-1);
else
    % if there aren't at least two peaks:
    plus_one_delta = nan;
    minus_one_delta = nan;
end

peak_num_delta = abs(length(sim_peaks)-length(data_peaks));

end

