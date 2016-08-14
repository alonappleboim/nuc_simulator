n_steps = 5000:5000:25000;
n_steps(6) = 50000;
distance_1_2 = zeros(6,10); % n_steps, sim_num
distance_2_3 = zeros(6,10); % n_steps, sim_num
width = zeros(6,10); % n_steps, sim_num
NFR_dist = zeros(6,10); % n_steps, sim_num

% run the simulations:
for i=1:6 % n_steps
    for j=1:10 % sim_num
        clc;
        fprintf('n_steps %i, sim_num %i',i,j);
        [s_hist_temp, tmp, distance_1_2(i,j), distance_2_3(i,j)] = run_simulation('n_steps',n_steps(i));
        [width(i,j), tmp] = get_peak_width(s_hist_temp, [1000:1200]);
        [NFR_dist(i,j), tmp] = get_NFR_width(s_hist_temp, [1000:1200]);
    end
end

% see how the number of steps effect the convergence:
figure
for i=1:6
    subplot(2,3,i)
    hold on
    plot(distance_1_2(i,:),'.')
    title(['+1+2 dist, steps = ' num2str(n_steps(i))])
    hold off
end
figure
for i=1:6
    subplot(2,3,i)
    hold on
    plot(distance_2_3(i,:),'.')
    title(['+2+3 dist, steps = ' num2str(n_steps(i))])
    hold off
end
figure
for i=1:6
    subplot(2,3,i)
    hold on
    plot(width(i,:),'.')
    title(['+1 width, steps = ' num2str(n_steps(i))])
    hold off
end
figure
for i=1:6
    subplot(2,3,i)
    hold on
    plot(NFR_dist(i,:),'.')
    title(['NFR width, steps = ' num2str(n_steps(i))])
    hold off
end
