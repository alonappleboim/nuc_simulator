features = zeros(1,6);

% Starting Variables:
FRS2_seq = sequences_structure(1001,:);
FRS2_wt = round(all.wt_3h(1001,:) ./ (10^9));

RSC_evic_intensity = [0.05, 0.1, 0.15];
%RSC_evic_intensity = 0.1.*ones(1,3);
%RSC_evic_length = [10, 20, 30];
RSC_evic_length = 20.*ones(1,3);
%RSC_slide_intensity = [2 , 4 , 6];
RSC_slide_intensity = 4.*ones(1,3);
%RSC_slide_length = [20, 40, 60];
RSC_slide_length = 40.*ones(1,3);

for j = 1:3

figure

for i=1:6
    % run the simulation:
    [full_s_hist, nuc_s_hist] = ...
        run_simulation_from_genome(FRS2_seq, 'linker_len', 10, 'n_steps', 25000, ...
        'RSC_evic_intensity', RSC_evic_intensity(j), 'RSC_evic_length', RSC_evic_length(j), ... 
        'RSC_slide_intensity', RSC_slide_intensity(j), 'RSC_slide_length', RSC_slide_length(j));
    
    features(i) = Compare_Sim_To_Data(nuc_s_hist(:, 701:1100), FRS2_wt(701:1100));
    sim_data = sum(nuc_s_hist(5000:end, :));
    sim_data = sim_data(701:1100);
    sim_data = (sim_data./sum(sim_data)) .* sum(FRS2_wt(701:1100));
    
    subplot(2,3,i)
    plot(FRS2_wt(701:1100),'g')
    hold on
    plot(sim_data, 'r')
end

for i=1:6
    subplot(2,3,i)
    title(['feature=' num2str(features(i)) ' RSC_evic_intensity=' num2str(RSC_evic_intensity(j))])
    xlabel(num2str(mean(features)))
end
end