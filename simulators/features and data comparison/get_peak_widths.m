function [ plus_one_width_delta, minus_one_width_delta ] ...
    = get_peak_widths( sim_centers, data_centers, NFR_pos, sim_plus1_pos, sim_minus1_pos, data_plus1_pos, data_minus1_pos )
%GET_PEAK_WIDTHS get the plus1 and minus1 peak width deltas between the sim
%and the data
%   The function defines width as the distance we need to move from the
%   peak to drop down by a factor of 3/4. The function finds the peaks and
%   then the widths, and returns the deltas.

if (sum(isnan([sim_plus1_pos, sim_minus1_pos, data_plus1_pos, data_minus1_pos]) > 0))
    plus_one_width_delta = nan;
    minus_one_width_delta = nan;
    
else

    % find the closest new peak to the original peak after this smoothing:
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

    %{
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
    %}

end

end

