n_steps = 5000:5000:25000;
n_steps(6) = 50000;
ks_widths = 5:5:25;
centers = zeros(6,5,3500); % n_steps, sim_num, s_hist_summation...
smoothed = zeros(6,5,3500,5); % the final 5 is for the smooth widths...
distance_1_2 = zeros(6,5,5,1); % n_steps, sim_num, width_used, distance...
distance_2_3 = zeros(6,5,5,1); % n_steps, sim_num, width_used, distance...

% run the simulations:
for i=1:6 % n_steps
    for j=1:5 % sim_num
        clc;
        fprintf('n_steps %i, sim_num %i',i,j);
        [s_hist_temp, tmp, dist_1_2, dist2_3] = run_simulation('n_steps',n_steps(i));
        centers(i,j,:) = sum(s_hist_temp);
        for k=1:5
            smoothed(i,j,:,k) = ksdensity(1:length(centers(i,j,:)),1:length(centers(i,j,:)),'weights',double(centers(i,j,:)),'width',ks_widths(k));
            [distance_1_2(i,j,k,1) , distance_2_3(i,j,k,1)] = get_peak_distances(squeeze(smoothed(i,j,1100:end,k)));
        end
    end
end

% see how each smoothing effects the convergence, for 20000 steps:
%{
figure
subplot(2,3,1)
hold on
plot(squeeze(centers(4,1,1100:1698)),'r')
plot(squeeze(centers(4,2,1100:1698)),'g')
plot(squeeze(centers(4,3,1100:1698)),'b')
plot(squeeze(centers(4,4,1100:1698)),'k')
plot(squeeze(centers(4,5,1100:1698)),'y')
title('Without Smoothing')
hold off

for i=1:5
    subplot(2,3,i+1)
    hold on
    plot(squeeze(smoothed(4,1,1100:1698,i)),'r')
    plot(squeeze(smoothed(4,2,1100:1698,i)),'g')
    plot(squeeze(smoothed(4,3,1100:1698,i)),'b')
    plot(squeeze(smoothed(4,4,1100:1698,i)),'k')
    plot(squeeze(smoothed(4,5,1100:1698,i)),'y')
    title(['Width = ' num2str(ks_widths(i))])
    hold off
end
%}

% see how the number of steps effect the convergence, for width 10:
figure
for i=1:6
    subplot(2,3,i)
    hold on
    plot(squeeze(smoothed(i,1,1100:1398,4)),'r')
    plot(squeeze(smoothed(i,2,1100:1398,4)),'g')
    plot(squeeze(smoothed(i,3,1100:1398,4)),'b')
    plot(squeeze(smoothed(i,4,1100:1398,4)),'k')
    plot(squeeze(smoothed(i,5,1100:1398,4)),'y')
    title(['steps = ' num2str(n_steps(i))])
    hold off
end
