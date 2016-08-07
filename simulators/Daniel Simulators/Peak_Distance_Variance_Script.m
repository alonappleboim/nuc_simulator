poly_params = 0:7;
distances = zeros(10,8);
ratios = zeros(10,8);
for i=1:8
    for j=1:10
    clc;
    fprintf('ITERATION %i-%i',i,j);
    [distances(j,i) , ratios(j,i), simulation] = simulate_polymerase(poly_params(i));
    end
end

meta_distance = get_average_peak_distance(smooth_metagene_5);
meta_ratio = get_peak_decline_ratio(smooth_metagene_5);

variance_distance = var(distances);
variance_ratio = var(ratios);

mean_distances = mean(distances);
mean_ratios = mean(ratios);

figure;
bar(poly_params,mean_distances)
hold on
errorbarxy(poly_params,mean_distances,0.001.*ones(size(poly_params)), variance_distance.^0.5,{'.','r','r'})
hold on
bar(8, meta_distance, 'g')
hold on
plot([-1 9],[meta_distance meta_distance], 'g')
hold on
title('distances')
hold off

figure;
bar(poly_params,mean_ratios)
hold on
errorbarxy(poly_params,mean_ratios,0.001.*ones(size(poly_params)), variance_ratio.^0.5,{'.','r','r'})
hold on
bar(8, meta_ratio, 'g')
hold on
plot([-1 9],[meta_ratio meta_ratio], 'g')
hold on
title('ratios')
hold off