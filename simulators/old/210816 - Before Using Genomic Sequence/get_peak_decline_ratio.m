function decline_ratio = get_peak_decline_ratio(smoothed)
% Given a smoothed vector of centers from a simulation (1500 bp vector),
% the function returns the ratio between the second peak and the first one.

% make the peaks and positions vectors:
[peaks,positions] = findpeaks(smoothed, 'MinPeakHeight', 0.8*10^-3);

% find the ratio:
decline_ratio = 1 - peaks(2)/peaks(1);