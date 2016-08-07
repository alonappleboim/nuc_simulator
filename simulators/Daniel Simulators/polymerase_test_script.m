poly_params = 0:0.5:7;
distances = zeros(1,15);
ratios = zeros(1,15);
simulations = zeros(3,1500,15);
for i=1:15
    clc;
    fprintf('At iteration %i',i);
    [dist_a , ratio_a, simulations(1,:,i)] = simulate_polymerase(poly_params(i));
    clc;
    fprintf('At iteration %i',i);
    [dist_b , ratio_b , simulations(2,:,i)] = simulate_polymerase(poly_params(i));
    clc;
    fprintf('At iteration %i',i);
    [dist_c , ratio_c , simulations(3,:,i)] = simulate_polymerase(poly_params(i));
    distances(i) = (dist_a+dist_b+dist_c)/3;
    ratios(i) = (ratio_a + ratio_b + ratio_c)/3;
end

meta_distance = get_average_peak_distance(smooth_metagene_5);
meta_ratio = get_peak_decline_ratio(smooth_metagene_5);

bar(poly_params,ratios)
hold on
bar(8, meta_ratio, 'g')
plot([-1 9],[meta_ratio meta_ratio], 'g')