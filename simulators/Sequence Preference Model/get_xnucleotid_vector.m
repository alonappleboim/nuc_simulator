function [ dinuc_vec ] = get_xnucleotid_vector( sequence, len )
%get_xnucleotid_vector given a gene sequence and length, the function returns a
%vector with 1s for the positions that have length A/Ts in a row, and 0 for the
%other positions.

% make vectors of A and T for the dinucleotid test:
AT_vec = zeros(size(sequence));
AT_vec(sequence == 'A') = 1;
AT_vec(sequence == 'T') = 1;

% find the places that have x A/T in a row:
dinuc_vec = conv(AT_vec, ones(1,len), 'same');
dinuc_vec(dinuc_vec < len) = 0;
dinuc_vec(dinuc_vec == len) = 1;

end

