
% Starting Variables:
FRS2_seq = sequences_structure(1001,:);
FRS2_wt = all.wt_3h(1001,:);

% run the simulation:
[full_s_hist, nuc_s_hist, peak_distance_1_2, peak_distance_2_3, NFR_width] = ...
    run_simulation_from_genome(FRS2_seq, 'linker_len', 40, 'n_steps', 25000);

centers_vector = sum(nuc_s_hist(:,:));
s_hist_coverage = ksdensity(1:length(centers_vector),1:length(centers_vector),'weights',double(centers_vector(1:end)),'width',10);

% plot wild type gene with PolyAT and the simulation (wt smoothed):
smoothed_wt = ksdensity([1:length(FRS2_wt)],[1:length(FRS2_wt)],'weights',double(FRS2_wt),'width',5);
figure;
plot(smoothed_wt(1:end-1),'b')
hold on
plot(s_hist_coverage(1:2500),'r')
legend('wild-type','simulation')