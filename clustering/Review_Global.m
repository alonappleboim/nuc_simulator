create_params_global_search;

over_2_count = zeros(1,length(params));
over_3_count = zeros(1,length(params));
sum_count = zeros(1,length(params));

for i = 1:length(params)
    load(['C:\Users\Daniel\Documents\MATLAB\Friedman Lab\results\wt_global\results_' num2str(i) '.mat'])
    over_2_count(i) = length(over_2);
    over_3_count(i) = length(over_3);
    sum_count(i) = sum_ratio;
end