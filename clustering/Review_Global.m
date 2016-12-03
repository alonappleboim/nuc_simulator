%%

create_params_global_search;

over_2_count = zeros(1,length(params));
over_3_count = zeros(1,length(params));
sum_count = zeros(1,length(params));
genes_over_3_per_set = zeros(length(params),54);

for i = 1:length(params)
    load(['C:\Users\Daniel\Documents\MATLAB\Friedman Lab\results\wt_global\results_' num2str(i) '.mat'])
    over_2_count(i) = length(over_2);
    over_3_count(i) = length(over_3);
    sum_count(i) = sum_ratio;
    genes_over_3_per_set(i, 1:length(over_3)) = over_3;
end

[~,over_3_best_indices] = sort(over_3_count,'descend');

% find the 10 parameter sets that add the most new genes:
parameter_sets = zeros(1,10);
current_set = [0];
added_per_step = zeros(1,10);
for i = 1:10
    
    % find the best current one (has the biggest count using setdiff)
    optimal_count = 0;
    optimal_index = 0;
    for j = 1:length(params)
        if (length(setdiff(genes_over_3_per_set(over_3_best_indices(j),:),current_set)) > optimal_count)
            optimal_count = length(setdiff(genes_over_3_per_set(over_3_best_indices(j),:),current_set));
            optimal_index = j;
        end
    end
    
    current_set = union(current_set,genes_over_3_per_set(over_3_best_indices(optimal_index),:));
    parameter_sets(i) = optimal_index;
    added_per_step(i) = optimal_count;
end

params(:,over_3_best_indices(parameter_sets(parameter_sets~=0)))
added_per_step