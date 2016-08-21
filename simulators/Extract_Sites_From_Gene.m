function [ PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites ] = Extract_Sites_From_Gene( genome )
%Extract_Data_From_Gene The function that extracts Poly(dA:dT) and binding
%sites from the genome.
%   Given a 2501-bp-long genome, where the TSS is at position 1000, the
%   function goes over the known trans-factors binding sites and finds
%   their positions in the genome, along with the Poly(dA:dT) positions.
%   For the binding sites I used the data from
%   http://yetfasco.ccbr.utoronto.ca/ and used two signatures for each
%   trans factor - a strong one and a wea one (the vector that is returned
%   returns 0.5 for weak binding site centers and 1 for strong ones. For
%   the PolyA and PolyT, I just return 1 for positions that are the center
%   of a 5-bp-long PolyA or PolyT.

% define the REGEXs for the binding sites:
REB1_Strong_Bind = '(CCGGGTAA)|(GGCCCATT)';
REB1_Weak_Bind = '(GGGT)|(CCCA)';
ABF1_Strong_Bind = '(TCAC.....ACG)|(AGTG.....TGC)';
ABF1_Weak_Bind = '(TC.......ACG)|(AG.......TGC)';
RAP1_Strong_Bind = '(TGT.TGGGTG)|(ACA.ACCCAC)';
RAP1_Weak_Bind = '(G...GGGT)|(C...CCCA)';

% define the vectors that will be returned:
PolyA_Sites = zeros(1,2501);
PolyT_Sites = zeros(1,2501);
REB1_Sites = zeros(1,2501);
ABF1_Sites = zeros(1,2501);
RAP1_Sites = zeros(1,2501);

% make vectors of PolyA and PolyT centers of 5 or longer lengths:
PolyA_Sites(genome == 'A') = 1;
PolyA_Sites = conv(PolyA_Sites, ones(1,7), 'same');
PolyA_Sites(PolyA_Sites < 5) = 0;
PolyA_Sites(PolyA_Sites > 4) = 1;
PolyA_Sites = conv(PolyA_Sites,ones(1,5),'same');
PolyA_Sites(PolyA_Sites > 0) = 1;

PolyT_Sites(genome == 'T') = 1;
PolyT_Sites = conv(PolyT_Sites, ones(1,7), 'same');
PolyT_Sites(PolyT_Sites < 5) = 0;
PolyT_Sites(PolyT_Sites > 4) = 1;
PolyT_Sites = conv(PolyT_Sites,ones(1,5),'same');
PolyT_Sites(PolyT_Sites > 0) = 1;

% find the center of binding sites:
[start_i, end_i] = regexp(genome, REB1_Weak_Bind);
REB1_Sites(fix((start_i + end_i)/2)) = 0.5;
[start_i, end_i] = regexp(genome, REB1_Strong_Bind);
REB1_Sites(fix((start_i + end_i)/2)) = 1;

[start_i, end_i] = regexp(genome, ABF1_Weak_Bind);
ABF1_Sites(fix((start_i + end_i)/2)) = 0.5;
[start_i, end_i] = regexp(genome, ABF1_Strong_Bind);
ABF1_Sites(fix((start_i + end_i)/2)) = 1;

[start_i, end_i] = regexp(genome, RAP1_Weak_Bind);
RAP1_Sites(fix((start_i + end_i)/2)) = 0.5;
[start_i, end_i] = regexp(genome, RAP1_Strong_Bind);
RAP1_Sites(fix((start_i + end_i)/2)) = 1;

% MAKE THE VECTOR BE 3500-LONG:
PolyA_Sites = [PolyA_Sites , zeros(1,3500-2501)];
PolyT_Sites = [PolyT_Sites , zeros(1,3500-2501)];
REB1_Sites = [REB1_Sites , zeros(1,3500-2501)];
ABF1_Sites = [ABF1_Sites , zeros(1,3500-2501)];
RAP1_Sites = [RAP1_Sites , zeros(1,3500-2501)];

end

