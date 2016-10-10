function [ plus_one, minus_one, NFR_width ] = get_NFR_features_data( centers )
%get_NFR_features_data a function for extracting NFR features from the experimental
%data.
%   assuming the TSS is at 1000, the function finds the +1 and -1
%   nucleosome positions and the NFR width by smoothing the centers vector 
%   and finding the peaks.

smoothed = ksdensity([1:length(centers)],[1:length(centers)],'weights',double(centers),'width',20);
[peaks,positions] = findpeaks(smoothed, 'MinPeakHeight', 2*10^(-4));

if (length(peaks) > 9)
    
    % get the top 10 peaks so we lose the small ones:
    [temp,indices] = sort(peaks,'descend');
    indices = sort(indices);
    indices = indices(1:10);
    positions = positions(indices);

    % find the +1 position:
    temp = positions(positions > 949);
    plus_one = temp(1);

    % find the -1 position:
    temp = positions(positions < 950);
    minus_one = temp(end);

    NFR_width = plus_one - minus_one;
   
else
    
    plus_one = 0;
    minus_one = 0;
    NFR_width = 0;
    
end

end

