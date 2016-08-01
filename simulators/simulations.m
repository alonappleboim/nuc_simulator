
%% %%% PARAMETERS %%%%%

params = struct();
gen_len = 8000;
params.a_rate = 10.*ones(1,gen_len); params.a_rate(3600:4400) = 0;
params.e_rate = ones(1,gen_len); params.e_rate(3600:4400) = 1;
params.r_rate = 10.*ones(1,gen_len); params.r_rate(3600:4400) = 0;
params.l_rate = 10.*ones(1,gen_len); params.l_rate(3600:4400) = 0;
params.nuc_width = 147;
params.slide_len = 20;
params.nuc_footprint = ones(1,(params.nuc_width.*2) - 1);
params.linker_len = 200;

%%

[time, s_hist] = gillespie(params, 'n_steps', 2e4, 's0', zeros(1,gen_len));

%% %%% FINAL CALCULATIONS %%%%%

% get the number of times each bp had a nuceosome center on it
centers_vector = sum(s_hist(:,:));
smooth_vector =  conv(centers_vector, ones(1,params.nuc_width), 'same');

% rescale graphs
smooth_vector = smooth_vector .* (max(centers_vector)/mean(smooth_vector));

%% %%% OUTPUT GRAPHS %%%%%

%{
figure;
plot(centers_vector)
title('Number of Nucleosome Centers VS Base Pair')
figure;
plot(smooth_vector)
title('Number of Nucleosome Coverage VS Base Pair')
%}
figure;
plot(centers_vector,'r')
hold on
plot(smooth_vector,'g')
legend('centers','coverage')
title('Combined')
xlabel('Position')
ylabel('Time')
hold off
