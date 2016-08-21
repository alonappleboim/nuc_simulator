
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

% run 5 simulation with the PolyATs as parameters:
[s_hist, s_hist_coverage1, peak_distance_1_2, peak_distance_2_3] = run_simulation_FRS2(FRS2_PolyAT, FRS2_PolyA, FRS2_PolyT, 'n_steps', 30000, 'linker_len', 40);
[s_hist, s_hist_coverage2, peak_distance_1_2, peak_distance_2_3] = run_simulation_FRS2(FRS2_PolyAT, FRS2_PolyA, FRS2_PolyT, 'n_steps', 30000, 'linker_len', 40);
[s_hist, s_hist_coverage3, peak_distance_1_2, peak_distance_2_3] = run_simulation_FRS2(FRS2_PolyAT, FRS2_PolyA, FRS2_PolyT, 'n_steps', 30000, 'linker_len', 40);
[s_hist, s_hist_coverage4, peak_distance_1_2, peak_distance_2_3] = run_simulation_FRS2(FRS2_PolyAT, FRS2_PolyA, FRS2_PolyT, 'n_steps', 30000, 'linker_len', 40);
[s_hist, s_hist_coverage5, peak_distance_1_2, peak_distance_2_3] = run_simulation_FRS2(FRS2_PolyAT, FRS2_PolyA, FRS2_PolyT, 'n_steps', 30000, 'linker_len', 40);

% plot the 5 simulations vs the real data:
smoothed_wt = ksdensity([1:length(FRS2_wt)],[1:length(FRS2_wt)],'weights',double(FRS2_wt),'width',5);
figure;
plot(smoothed_wt(1:end-1),'b')
hold on
plot(s_hist_coverage1(1:2500),'r')
plot(s_hist_coverage2(1:2500),'g')
plot(s_hist_coverage3(1:2500),'y')
plot(s_hist_coverage4(1:2500),'k')
plot(s_hist_coverage5(1:2500),'m')
