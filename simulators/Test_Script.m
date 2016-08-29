
% Starting Variables:
FRS2_seq = sequences_structure(1001,:);
FRS2_wt =  wt_3h(1001,:); %all.wt_3h(1001,:);

[ PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites ] = ...
    Extract_Sites_From_Gene(FRS2_seq, 3500);

% run the simulation:
[nuc_sum, time, nuc_s_hist, REB1_s_hist, ABF1_s_hist, RAP1_s_hist] = ...
    run_simulation_from_genome(FRS2_seq, 'linker_len', 10, 'n_steps', 25000, ...
    'poly_rate', 0, 'REB1_a_rate', 0, 'REB1_e_rate', 0, 'ABF1_a_rate', 0, ...
                    'ABF1_e_rate', 0, 'RAP1_a_rate', 0, 'RAP1_e_rate', 0,...
                    'TF_evic_intensity', 0.001, 'RSC_evic_intensity', 0.001, ...
                    'RSC_evic_length', 20, 'RSC_slide_length', 40, ...
                    'RSC_slide_intensity', 6);
                
centers_vector = nuc_sum;
s_hist_coverage = ksdensity(1:length(centers_vector),1:length(centers_vector),'weights',double(centers_vector(1:end)),'width',5);

% plot wild type gene with PolyAT and the simulation (wt smoothed):
smoothed_wt = ksdensity([1:length(FRS2_wt)],[1:length(FRS2_wt)],'weights',double(FRS2_wt),'width',5);

figure;
plot(smoothed_wt(1:end-1),'c')
%plot(FRS2_wt(1:2500) ./ sum(FRS2_wt),'b')
hold on
plot(s_hist_coverage(1:2500),'r')
%plot(centers_vector(1:2500) ./ sum(centers_vector(1:2500)),'r')
plot(PolyA_Sites(1:2500) .* mean(smoothed_wt),'k')
plot(PolyT_Sites(1:2500) .* mean(smoothed_wt),'m')
plot(REB1_Sites(1:2500) .* 4 .* mean(smoothed_wt), 'g')
plot(ABF1_Sites(1:2500) .* 4 .* mean(smoothed_wt), 'b')
plot(RAP1_Sites(1:2500) .* 4 .* mean(smoothed_wt), 'y')
legend('wild-type','simulation','PolyA (right)','PolyT (left)', 'REB1', 'ABF1', 'RAP1')
xlabel('Position')
ylabel('Intensity')

nuc_sum = nuc_sum(701:1100) ./ sum(nuc_sum(701:1100));
nuc_sum = nuc_sum .* sum(FRS2_wt(701:1100));
figure;
plot(FRS2_wt(701:1100),'g')
hold on
plot(nuc_sum, 'r')