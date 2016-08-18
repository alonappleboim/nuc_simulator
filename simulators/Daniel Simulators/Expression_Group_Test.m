
metagene_split = zeros(5,2501);
metagene_split_smooth = zeros(5,2501);
distances_1_2 = zeros(1,5);
distances_2_3 = zeros(1,5);

% split the metagene into expression levels:
metagene_split(1,:) = nansum(all.wt_3h((midlog_expression > 1000),:));
metagene_split(2,:) = nansum(all.wt_3h((midlog_expression < 1000 & midlog_expression > 500),:));
metagene_split(3,:) = nansum(all.wt_3h((midlog_expression < 500 & midlog_expression > 300),:));
metagene_split(4,:) = nansum(all.wt_3h((midlog_expression < 300 & midlog_expression > 200),:));
metagene_split(5,:) = nansum(all.wt_3h((midlog_expression < 200),:));

% normalise and smooth the metagene:
for i=1:5
    metagene_split(i,:) = metagene_split(i,:) ./ sum(metagene_split(i,:));
    metagene_split_smooth(i,:) = ksdensity(1:length(metagene_split(i,:)),1:length(metagene_split(i,:)),'weights',double(metagene_split(i,:)),'width',20);
    metagene_split_smooth(i,:) = metagene_split_smooth(i,:) ./ sum(metagene_split_smooth(i,:));
end

% find the peak distances of the metagene:
for i=1:5
    [distances_1_2(i),distances_2_3(i)] = get_metagene_peak_distances(metagene_split_smooth(i,1000:end-1));
end

% run simulations with different poly rates:
poly_rates = 8.5:-2:0.5;
sim_distances_1_2 = zeros(1,5);
sim_distances_2_3 = zeros(1,5);
for i=1:5
    [tmp1,tmp2,sim_distances_1_2(i),sim_distances_2_3(i)] = run_simulation('n_steps', 25000, 'poly_rate', poly_rates(i),'linker_len',1);
end


%{
figure
hold on
plot(metagene_split_smooth(1,:), 'r')
plot(metagene_split_smooth(2,:), 'g')
plot(metagene_split_smooth(3,:), 'b')
plot(metagene_split_smooth(4,:), 'y')
plot(metagene_split_smooth(5,:), 'k')
%}

figure
hold on
plot(distances_1_2,'*g')
plot(sim_distances_1_2,'sr')
errorbarxy([1:5], sim_distances_1_2, 0.001.*ones(1,5), 2.*ones(1,5),{'.','r','r'});
title('Distance between +1 and +2 nucs - Poly rates and gene expression rates')
xlabel(['Metagene - 5 expression groups from large to small' char(10) 'Simulation - poly rates from 8.5 to 0.5 in jumps of 2'])
ylabel('Distance Between Peaks')
legend('metagene','simulation','errorbar')

figure
hold on
plot(distances_2_3,'*g')
plot(sim_distances_2_3,'sr')
errorbarxy([1:5], sim_distances_2_3, 0.001.*ones(1,5), 2.*ones(1,5),{'.','r','r'});
title('Distance between +2 and +3 nucs - Poly rates and gene expression rates')
xlabel(['Metagene - 5 expression groups from large to small' char(10) 'Simulation - poly rates from 8.5 to 0.5 in jumps of 2'])
ylabel('Distance Between Peaks')
legend('metagene','simulation')

