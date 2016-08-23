features = zeros(1,6);

% Starting Variables:
FRS2_seq = sequences_structure(1001,:);
FRS2_wt = round(all.wt_3h(1001,:) ./ (10^9));

figure

for i=1:6
    % run the simulation:
    [full_s_hist, nuc_s_hist] = ...
        run_simulation_from_genome(FRS2_seq, 'linker_len', 10, 'n_steps', 25000);
    
    features(i) = Compare_Sim_To_Data(nuc_s_hist(:, 701:1100), FRS2_wt(701:1100));
    sim_data = sum(nuc_s_hist(5000:end, :));
    sim_data = sim_data(701:1100);
    sim_data = (sim_data./sum(sim_data)) .* sum(FRS2_wt(701:1100));
    
    if (i<7)
        subplot(2,3,i)
        plot(FRS2_wt(701:1100),'g')
        hold on
        plot(sim_data, 'r')
    end  
end

for i=1:6
    subplot(2,3,i)
    title(num2str(features(i)))
end

figure
plot(features,'.')