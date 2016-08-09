linker_lengths = 1:10;
distances_1_2 = zeros(10,10);
distances_2_3 = zeros(10,10);
ratios = zeros(10,10);

for i=1:10
    for j=1:10
    clc;
    fprintf('ITERATION %i-%i',i,j);
    [tmp1, tmp2, distances_1_2(j,i), distances_2_3(j,i)] = run_simulation('poly_rate',linker_lengths(i));
    end
end

[meta_distance_1_2, meta_distance_2_3] = get_peak_distances(smooth_metagene_5);

variance_distance_1_2 = var(distances_1_2);
variance_distance_2_3 = var(distances_2_3);
mean_distances_1_2 = mean(distances_1_2);
mean_distances_2_3 = mean(distances_2_3);

figure;
subplot(1,2,1)
bar(linker_lengths,mean_distances_1_2)
hold on
errorbarxy(linker_lengths,mean_distances_1_2,0.001.*ones(size(linker_lengths)), variance_distance_1_2.^0.5,{'.','r','r'})
hold on
bar(11, meta_distance_1_2, 'g')
hold on
plot([1 11],[meta_distance_1_2 meta_distance_1_2], 'g')
hold on
title('Distance between +1 and +2 Peaks VS Polymerase Parameter')
xlabel('Polymerase Parameter (left slide rate)')
ylabel('Peak Distance')

subplot(1,2,2)
bar(linker_lengths,mean_distances_2_3)
hold on
errorbarxy(linker_lengths,mean_distances_2_3,0.001.*ones(size(linker_lengths)), variance_distance_2_3.^0.5,{'.','r','r'})
hold on
bar(11, meta_distance_2_3, 'g')
hold on
plot([1 11],[meta_distance_2_3 meta_distance_2_3], 'g')
hold on
title('Distance between +2 and +3 Peaks VS Polymerase Parameter')
xlabel('Polymerase Parameter (left slide rate)')
ylabel('Peak Distance')

annotation('textbox',[0 0.8 0.2 0.2], ...
			'String',['genome length = 3500' char(10) ...
						'number of steps per simulation = 20000' char(10) ...
						'polymerase parameter = 3' char(10) ...
						'NFR location = 1000:1200' char(10) ...
						'Polymerase location = 1200:2800'])

hold off
