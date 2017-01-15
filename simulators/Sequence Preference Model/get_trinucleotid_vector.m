function [ dinuc_vec ] = get_trinucleotid_vector( sequence )
%get_trinucleotid_vector given a gene sequence, the function returns a
%vector with 1s for the positions that have two A/T in a row, and 0 for the
%other positions.

% make vectors of A and T for the dinucleotid test:
AT_vec = zeros(size(sequence));
AT_vec(sequence == 'A') = 1;
AT_vec(sequence == 'T') = 1;

% find the places that have two A/T in a row:
dinuc_vec = conv(AT_vec, ones(1,3), 'same');
dinuc_vec(dinuc_vec < 3) = 0;
dinuc_vec(dinuc_vec == 3) = 1;

end

