%%

Gene_id = 9;
genlen = 3500;
TSS = round(genlen/2);
NFR_pos = [TSS-399 : TSS+200];

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

centers_vector = zeros(20,genlen);

% run the simulation:
for i = 1:20
    [nuc_sum, time, nuc_s_hist, REB1_s_hist, ABF1_s_hist, RAP1_s_hist] = ...
        run_simulation_from_genome(FRS2_seq, 'linker_len', 10, 'n_steps', 25000, ...
        'poly_rate', 0, 'REB1_a_rate', 0.0001, 'REB1_e_rate', 0.0001, 'ABF1_a_rate', 0.0001, ...
                    'ABF1_e_rate', 0.0001, 'RAP1_a_rate', 0.0001, 'RAP1_e_rate', 0.0001,...
                    'TF_evic_intensity', 0.0001, 'RSC_evic_intensity', 0.05, ...
                    'RSC_evic_length', 60, 'RSC_slide_length', 120, ...
                    'RSC_slide_intensity', 0.5, 'gen_len', genlen, 'slide_len', 3);
                
    centers_vector(i,:) = nuc_sum;
end

% plot wild type gene with PolyAT and the simulation (wt smoothed):
smoothed_wt = zeros(20,genlen);

%%

FRS_smooth = 0;
FRS_smooth = conv(FRS2_wt, gausswin(100), 'same');
FRS_smooth = conv(FRS_smooth, gausswin(20), 'same');
[FRS_peaks,FRS_positions] = findpeaks(FRS_smooth, 'MinPeakHeight', sum(FRS_smooth)./length(FRS_smooth));

%%
for i = 1:20
    smoothed_wt(i,:) = conv(centers_vector(i,:), gausswin(100), 'same');
    smoothed_wt(i,:) = conv(smoothed_wt(i,:), gausswin(20), 'same');
end

%%

for i = 1:20
    figure;
    plot(NFR_pos,smoothed_wt(i,NFR_pos),'b')
    hold on
    FRS_smooth = FRS_smooth .* (sum(smoothed_wt(i,:)) / sum(FRS_smooth));
    plot(NFR_pos,FRS_smooth(NFR_pos),'g')

    [peaks,positions] = findpeaks(smoothed_wt(i,:), 'MinPeakHeight', (sum(smoothed_wt(i,:)))./(1.5*length(smoothed_wt(i,:))));
    [FRS_peaks,FRS_positions] = findpeaks(FRS_smooth, 'MinPeakHeight', (sum(FRS_smooth))./(1.5*length(FRS_smooth)));

    plot(positions((positions > NFR_pos(1)) & (positions < NFR_pos(end))), peaks((positions > NFR_pos(1)) & (positions < NFR_pos(end))), '^r')
    plot(FRS_positions((FRS_positions > NFR_pos(1)) & (FRS_positions < NFR_pos(end))), FRS_peaks((FRS_positions > NFR_pos(1)) & (FRS_positions < NFR_pos(end))), '^k')
end

%%