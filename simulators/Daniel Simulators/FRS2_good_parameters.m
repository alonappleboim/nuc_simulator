
% Starting Variables:
FRS2_seq = sequences_structure(1001,:);
FRS2_PolyAT = zeros(1,2501);
FRS2_PolyA = zeros(1,2501);
FRS2_PolyT = zeros(1,2501);
FRS2_PolyT_5temp = zeros(1,2501);
FRS2_PolyT_7temp = zeros(1,2501);
FRS2_wt = all.wt_3h(1001,:);

% make vectors of PolyA and PolyT of 5 or larger lengths:
FRS2_PolyA(FRS2_seq == 'A') = 1;
FRS2_PolyT(FRS2_seq == 'T') = 1;
FRS2_PolyA = conv(FRS2_PolyA, ones(1,7), 'same');
FRS2_PolyT = conv(FRS2_PolyT, ones(1,7), 'same');

FRS2_PolyA(FRS2_PolyA < 5) = 0;
FRS2_PolyA(FRS2_PolyA > 4) = 1;
FRS2_PolyA = conv(FRS2_PolyA,ones(1,5),'same');
FRS2_PolyA(FRS2_PolyA > 0) = 1;
FRS2_PolyT(FRS2_PolyT < 5) = 0;
FRS2_PolyT(FRS2_PolyT > 4) = 1;
FRS2_PolyT = conv(FRS2_PolyT,ones(1,5),'same');
FRS2_PolyT(FRS2_PolyT > 0) = 1;
FRS2_PolyAT(FRS2_PolyA == 1) = 1;
FRS2_PolyAT(FRS2_PolyT == 1) = 1;

% run a simulation with the PolyATs as parameters:
[s_hist, s_hist_coverage, peak_distance_1_2, peak_distance_2_3] = run_simulation_FRS2(FRS2_PolyAT, FRS2_PolyA, FRS2_PolyT, 'n_steps', 20000, 'linker_len', 40);

% plot wild type gene with PolyAT and the simulation (wt smoothed):
smoothed_wt = ksdensity([1:length(FRS2_wt)],[1:length(FRS2_wt)],'weights',double(FRS2_wt),'width',5);
figure;
plot(smoothed_wt(1:end-1),'g')
hold on
plot(s_hist_coverage(1:2500),'b')
plot(FRS2_PolyAT(1:end-1) .* mean(smoothed_wt(1:end-1)), 'r')
legend('wild-type','simulation','Poly(dA:dT)')
