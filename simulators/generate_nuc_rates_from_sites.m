function [ e_rate, r_rate, l_rate ] = ...
    generate_nuc_rates_from_sites( PolyA_sites, PolyT_sites, REB1_sites, ...
    ABF1_sites, RAP1_sites, nuc_width, REB1_width, ABF1_width, RAP1_width )
%generate_nuc_rates_from_sites generates the nucleosome eviction, right and
%left rate vectors.
%   given the PolyA, PolyT and trans sites, along with the nuc and trans widths,
%   the function returns the relevant rates for the nucleosomes.

XXX = 0;

PolyAT_sites = PolyA_sites + PolyT_sites;

% make the base rates:
e_rate = ones(size(REB1_sites)); 
r_rate = 0.1.*ones(size(REB1_sites)); 
l_rate = 0.1.*ones(size(REB1_sites)); 

% make the eviction rates from the sites:
RSC_evict = conv(PolyAT_sites, 0.2.*ones(1,20),'same');
REB1_evict = conv(REB1_sites, XXX.*ones(1,REB1_width + nuc_width),'same');
ABF1_evict = conv(ABF1_sites, XXX.*ones(1,ABF1_width + nuc_width),'same');
RAP1_evict = conv(RAP1_sites, XXX.*ones(1,RAP1_width + nuc_width),'same');
e_rate = e_rate + RSC_evict + REB1_evict + ABF1_evict + RAP1_evict;

% make the left-right sliding rates from the sites:
PolyAT_right = conv(PolyA_sites, 1.*ones(1,50),'same');
PolyAT_left  = conv(PolyT_sites, 1.*ones(1,50),'same');
r_rate = r_rate + PolyAT_right;
l_rate = l_rate + PolyAT_left;

end

