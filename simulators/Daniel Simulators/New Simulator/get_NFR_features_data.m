function [ plus_one, minus_one, NFR_width ] = get_NFR_features_data( centers )
%get_NFR_features_data a function for extracting NFR features from the experimental
%data.
%   assuming the TSS is at 1000, the function finds the +1 and -1
%   nucleosome positions and the NFR width by smoothing the centers vector 
%   and finding the peaks.

smoothed = ksdensity([1:length(centers)],[1:length(centers)],'weights',double(centers),'width',10);
[peaks,positions] = findpeaks(smoothed,'MinPeakHeight',10^-3);

% find the +1 position:
temp = positions(positions > 950);
plus_one = temp(1);

% find the -1 position:
temp = positions(positions < 950);
minus_one = temp(end);

NFR_width = plus_one - minus_one;

end

