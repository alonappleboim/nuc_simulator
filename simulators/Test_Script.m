
% Starting Variables:
FRS2_seq = sequences_structure(1001,:);
FRS2_wt = all.wt_3h(1001,:);

[ PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites ] = ...
    Extract_Sites_From_Gene(FRS2_seq, 3500);

% run the simulation:
[full_s_hist, nuc_s_hist] = ...
    run_simulation_from_genome(FRS2_seq, 'linker_len', 10, 'n_steps', 25000);

centers_vector = sum(nuc_s_hist(5000:end,:));
s_hist_coverage = ksdensity(1:length(centers_vector),1:length(centers_vector),'weights',double(centers_vector(1:end)),'width',5);

% plot wild type gene with PolyAT and the simulation (wt smoothed):
smoothed_wt = ksdensity([1:length(FRS2_wt)],[1:length(FRS2_wt)],'weights',double(FRS2_wt),'width',5);
figure;
plot(smoothed_wt(1:end-1),'b')
%plot(FRS2_wt(1:2500) ./ sum(FRS2_wt),'b')
hold on
plot(s_hist_coverage(1:2500),'r')
%plot(centers_vector(1:2500) ./ sum(centers_vector(1:2500)),'r')
plot(PolyA_Sites(1:2500) .* mean(smoothed_wt),'k')
plot(PolyT_Sites(1:2500) .* mean(smoothed_wt),'m')
legend('wild-type','simulation','PolyA (right)','PolyT (left)')