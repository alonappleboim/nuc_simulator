features = zeros(1,6);

% Starting Variables:
FRS2_seq = sequences_structure(1001,:);
FRS2_wt = round(wt_3h(1001,:));

%RSC_evic_intensity = [0.05, 0.1, 0.15];
RSC_evic_intensity = 0.1.*ones(1,3);
%RSC_evic_length = [10, 20, 30];
RSC_evic_length = 20.*ones(1,3); % 20
%RSC_slide_intensity = [2 , 4 , 6];
RSC_slide_intensity = 4.*ones(1,3);
%RSC_slide_length = [20, 40, 60];
RSC_slide_length = 40.*ones(1,3); %40

figure
for i=1:6
    % run the simulation:
    [nuc_sum1] = ...
        run_simulation_from_genome(FRS2_seq, 'linker_len', 10, 'n_steps', 25000, ...
        'RSC_evic_intensity', RSC_evic_intensity(1), 'RSC_evic_length', RSC_evic_length(1), ... 
        'RSC_slide_intensity', RSC_slide_intensity(1), 'RSC_slide_length', RSC_slide_length(1));
    [nuc_sum2] = ...
        run_simulation_from_genome(FRS2_seq, 'linker_len', 10, 'n_steps', 25000, ...
        'RSC_evic_intensity', RSC_evic_intensity(1), 'RSC_evic_length', RSC_evic_length(1), ... 
        'RSC_slide_intensity', RSC_slide_intensity(1), 'RSC_slide_length', RSC_slide_length(1));
    [nuc_sum3] = ...
        run_simulation_from_genome(FRS2_seq, 'linker_len', 10, 'n_steps', 25000, ...
        'RSC_evic_intensity', RSC_evic_intensity(1), 'RSC_evic_length', RSC_evic_length(1), ... 
        'RSC_slide_intensity', RSC_slide_intensity(1), 'RSC_slide_length', RSC_slide_length(1));
    [nuc_sum4] = ...
        run_simulation_from_genome(FRS2_seq, 'linker_len', 10, 'n_steps', 25000, ...
        'RSC_evic_intensity', RSC_evic_intensity(1), 'RSC_evic_length', RSC_evic_length(1), ... 
        'RSC_slide_intensity', RSC_slide_intensity(1), 'RSC_slide_length', RSC_slide_length(1));
    [nuc_sum5] = ...
        run_simulation_from_genome(FRS2_seq, 'linker_len', 10, 'n_steps', 25000, ...
        'RSC_evic_intensity', RSC_evic_intensity(1), 'RSC_evic_length', RSC_evic_length(1), ... 
        'RSC_slide_intensity', RSC_slide_intensity(1), 'RSC_slide_length', RSC_slide_length(1));
    
    nuc_sum = nuc_sum1 + nuc_sum2 + nuc_sum3 + nuc_sum4 + nuc_sum5;
    features(i) = Compare_Sum_To_Data(nuc_sum(:, 701:1100), FRS2_wt(701:1100));
    
    sim_data = nuc_sum(701:1100);
    sim_data = (sim_data./sum(sim_data)) .* sum(FRS2_wt(701:1100));
    
    subplot(2,3,i)
    plot(FRS2_wt(701:1100),'g')
    hold on
    plot(sim_data, 'r')
end

for i=1:6
    subplot(2,3,i)
    title(['feature=' num2str(features(i))])
    xlabel(num2str(mean(features)))
end