function [ plus_one_delta, minus_one_delta, peak_num_delta, ... 
    plus_one_width_delta, minus_one_width_delta, peak_ratio_delta ] ...
        = get_NFR_features_data( sim_centers, data_centers, NFR_pos )
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

%%

%%% find the PEAK LOCATIONS of the data and the simulation:
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

%%

%%% find the PEAK WIDTHS of the simulation and the data:

sim_plus1_pos = sim_positions(end);
sim_minus1_pos = sim_positions(end-1);
data_plus1_pos = data_positions(end);
data_minus1_pos = data_positions(end-1);

% find the closest new peak to the original peak (there is less smoothing now):
sim_smooth = conv(sim_centers, gausswin(50), 'same');
[~,sim_positions] = findpeaks(sim_smooth, 'MinPeakHeight', sum(sim_smooth)./(1.5*length(sim_smooth)));
[~, index] = min(abs(sim_plus1_pos - sim_positions));
sim_plus1_pos = sim_positions(index);
[~, index] = min(abs(sim_minus1_pos - sim_positions));
sim_minus1_pos = sim_positions(index);

data_smooth = conv(data_centers, gausswin(50), 'same');
[~,data_positions] = findpeaks(data_smooth, 'MinPeakHeight', sum(data_smooth)./(1.5*length(data_smooth)));
[~, index] = min(abs(data_plus1_pos - data_positions));
data_plus1_pos = data_positions(index);
[~, index] = min(abs(data_minus1_pos - data_positions));
data_minus1_pos = data_positions(index);

% define the distance to where we drop down by a quarter as the width:
sim_plus1_drop = sim_smooth(sim_plus1_pos) .* (3/4);
sim_minus1_drop = sim_smooth(sim_minus1_pos) .* (3/4);
data_plus1_drop = data_smooth(data_plus1_pos) .* (3/4);
data_minus1_drop = data_smooth(data_minus1_pos) .* (3/4);

% plus1 widths calculation:
temp = sim_smooth(1:sim_plus1_pos);
temp(temp > sim_plus1_drop) = 0;
temp = find(temp);
if (isempty(temp))
    plus_sim_left_index = NFR_pos(1);
else
    plus_sim_left_index = temp(end);
end

temp = data_smooth(1:data_plus1_pos);
temp(temp > data_plus1_drop) = 0;
temp = find(temp);
if (isempty(temp))
    plus_data_left_index = NFR_pos(1);
else
    plus_data_left_index = temp(end);
end 

temp = sim_smooth(sim_plus1_pos:end);
temp(temp > sim_plus1_drop) = 0;
temp = find(temp);
if (isempty(temp))
    plus_sim_right_index = NFR_pos(end);
else
    plus_sim_right_index = temp(1) + sim_plus1_pos;
end 

temp = data_smooth(data_plus1_pos:end);
temp(temp > data_plus1_drop) = 0;
temp = find(temp);
if (isempty(temp))
    plus_data_right_index = NFR_pos(end);
else
    plus_data_right_index = temp(1) + data_plus1_pos;
end 

plus_one_width_delta = abs(abs(plus_sim_right_index - plus_sim_left_index) - abs(plus_data_right_index - plus_data_left_index));

% minus1 widths calculation:
temp = sim_smooth(1:sim_minus1_pos);
temp(temp > sim_minus1_drop) = 0;
temp = find(temp);
if (isempty(temp))
    minus_sim_left_index = NFR_pos(1);
else
    minus_sim_left_index = temp(end);
end

temp = data_smooth(1:data_minus1_pos);
temp(temp > data_minus1_drop) = 0;
temp = find(temp);
if (isempty(temp))
    minus_data_left_index = NFR_pos(1);
else
    minus_data_left_index = temp(end);
end

temp = sim_smooth(sim_minus1_pos:end);
temp(temp > sim_minus1_drop) = 0;
temp = find(temp);
if (isempty(temp))
    minus_sim_right_index = NFR_pos(end);
else
    minus_sim_right_index = temp(1) + sim_minus1_pos;
end

temp = data_smooth(data_minus1_pos:end);
temp(temp > data_minus1_drop) = 0;
temp = find(temp);
if (isempty(temp))
    minus_data_right_index = NFR_pos(end);
else
    minus_data_right_index = temp(1) + data_minus1_pos;
end

minus_one_width_delta = abs(abs(minus_sim_right_index - minus_sim_left_index) - abs(minus_data_right_index - minus_data_left_index));

figure
plot(NFR_pos, sim_smooth(NFR_pos), 'r')
hold on
plot(NFR_pos, data_smooth(NFR_pos) .* (sum(sim_smooth(NFR_pos)) ./ sum(data_smooth(NFR_pos))), 'b')
plot(plus_sim_left_index, sim_plus1_drop, '^r')
plot(plus_sim_right_index, sim_plus1_drop, '^r')
plot(minus_sim_left_index, sim_minus1_drop, '^r')
plot(minus_sim_right_index, sim_minus1_drop, '^r')
plot(plus_data_left_index, data_plus1_drop.* (sum(sim_smooth(NFR_pos)) ./ sum(data_smooth(NFR_pos))), '^b')
plot(plus_data_right_index, data_plus1_drop.* (sum(sim_smooth(NFR_pos)) ./ sum(data_smooth(NFR_pos))), '^b')
plot(minus_data_left_index, data_minus1_drop.* (sum(sim_smooth(NFR_pos)) ./ sum(data_smooth(NFR_pos))), '^b')
plot(minus_data_right_index, data_minus1_drop.* (sum(sim_smooth(NFR_pos)) ./ sum(data_smooth(NFR_pos))), '^b')
legend('sim','wt')

%%

%%% find the PEAK HEIGHT RATIO feature:
peak_ratio_delta = 0;

end

