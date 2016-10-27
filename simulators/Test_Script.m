
Gene_id = 23;
genlen = 3500;
TSS = round(genlen/2);
NFR_pos = [TSS-299 : TSS+150];

% Starting Variables:
FRS2_seq = sequences_structure(Gene_id,:);
FRS2_wt =  wt_3h(Gene_id,:);

% make the wt data the right length:
buffer = genlen - 2501;
right_buffer = fix((buffer-500)/2);
left_buffer = right_buffer + 500;
if (right_buffer + left_buffer < buffer)
    left_buffer = left_buffer + 1;
end
FRS2_wt = [zeros(1,left_buffer), FRS2_wt, zeros(1,right_buffer)];

[ PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites ] = ...
    Extract_Sites_From_Gene(FRS2_seq, genlen);

centers_vector = zeros(1,genlen);

% run the simulation:
for i = 1:1
    [nuc_sum, time, nuc_s_hist, REB1_s_hist, ABF1_s_hist, RAP1_s_hist] = ...
        run_simulation_from_genome(FRS2_seq, 'linker_len', 10, 'n_steps', 25000, ...
        'poly_rate', 0, 'REB1_a_rate', 0.0001, 'REB1_e_rate', 0.0001, 'ABF1_a_rate', 0.0001, ...
                    'ABF1_e_rate', 0.0001, 'RAP1_a_rate', 0.0001, 'RAP1_e_rate', 0.0001,...
                    'TF_evic_intensity', 0.0001, 'RSC_evic_intensity', 0.05, ...
                    'RSC_evic_length', 60, 'RSC_slide_length', 120, ...
                    'RSC_slide_intensity', 0.5, 'gen_len', genlen, 'slide_len', 3);
                
    centers_vector = centers_vector + nuc_sum;
end

s_hist_coverage = ksdensity(1:length(centers_vector),1:length(centers_vector),'weights',double(centers_vector(1:end)),'width',5);

% plot wild type gene with PolyAT and the simulation (wt smoothed):
smoothed_wt = ksdensity([1:length(FRS2_wt)],[1:length(FRS2_wt)],'weights',double(FRS2_wt),'width',5);

%{
%%%%%%%%%%%%%%%%%%%%%%%%
% trying out the MNase digestion:
window = Get_MNase_Window(20);
centers_vector = conv(centers_vector, window, 'same');
%%%%%%%%%%%%%%%%%%%%%%%%
%}

[likelihood, plus1, minus1] = Compare_Sum_To_Data(centers_vector, FRS2_wt, NFR_pos, true)
%plus_one_turnover(nuc_s_hist, time, [700,900], 5)

figure;
plot(smoothed_wt(1:end-1),'b')
%plot(FRS2_wt(1:2500) ./ sum(FRS2_wt),'b')
hold on
plot(s_hist_coverage,'r')
%plot(centers_vector(1:2500) ./ sum(centers_vector(1:2500)),'r')
plot(PolyA_Sites .* mean(smoothed_wt),'k')
plot(PolyT_Sites .* mean(smoothed_wt),'m')
plot(REB1_Sites .* 4 .* mean(smoothed_wt), 'g')
plot(ABF1_Sites .* 4 .* mean(smoothed_wt), 'c')
plot(RAP1_Sites .* 4 .* mean(smoothed_wt), 'y')
legend('wild-type','simulation','PolyA (right)','PolyT (left)', 'REB1', 'ABF1', 'RAP1')
xlabel(['Position (TSS at ' num2str(fix(genlen/2)) ')'])
ylabel('Intensity')

%{
% plot the NFR smoothed:
s_hist_coverage = ksdensity(1:length(FRS2_wt),1:length(FRS2_wt),'weights',double(centers_vector(1:length(FRS2_wt))),'width',5);
smoothed_wt = ksdensity([1:length(FRS2_wt)],[1:length(FRS2_wt)],'weights',double(FRS2_wt),'width',5);
figure;
plot(smoothed_wt(NFR_pos) .* sum(FRS2_wt(NFR_pos)),'g')
hold on
plot(s_hist_coverage(NFR_pos) .* sum(FRS2_wt(NFR_pos)), 'r')
%}

% plot the NFR normal:
centers_vector = conv(centers_vector,gausswin(5)./sum(gausswin(5)),'same');
centers_vector = centers_vector(NFR_pos) ./ sum(centers_vector(NFR_pos));
centers_vector = centers_vector .* sum(FRS2_wt(NFR_pos));
FRS2_wt = conv(FRS2_wt,gausswin(5)./sum(gausswin(5)),'same');
figure;
plot(FRS2_wt(NFR_pos),'b')
hold on
plot(centers_vector, 'r')
legend('wild-type', 'simulation')
