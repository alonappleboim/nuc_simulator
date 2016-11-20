load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\sequences_structure.mat')

NFR = [800:1100];
genes = find(~(sequences_structure(:,1) == ' '));
As = zeros(1,length(genes));
Ts = zeros(1,length(genes));
j=0;

for i = genes'
    j = j+1;
    genome = sequences_structure(i,NFR);
    PolyA_Sites = zeros(1,length(NFR));
    PolyT_Sites = zeros(1,length(NFR));
    
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
    
    As(j) = sum(PolyA_Sites);
    Ts(j) = sum(PolyT_Sites);
end

low_poly_genes = find(As+Ts < quantile(As+Ts,0.05));
high_poly_genes = find(As+Ts > quantile(As+Ts,0.95));
high_A_low_T_genes = find(As > quantile(As,0.8) & Ts < quantile(Ts,0.2));
low_A_high_T_genes = find(Ts > quantile(Ts,0.8) & As < quantile(As,0.2));
