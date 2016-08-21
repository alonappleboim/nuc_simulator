
%% %%% PARAMETERS %%%%%

n_steps = 40000;
params = struct();
gen_len = 4000;
params.nuc_width = 147;
params.slide_len = 20;
params.a_rate = ones(1,gen_len); 
params.e_rate = ones(1,gen_len);
params.r_rate = 5.*ones(1,gen_len); 
params.l_rate = 5.*ones(1,gen_len); 
params.nuc_footprint = ones(1,(params.nuc_width.*2) - 1);
params.linker_len = 1;

% NFRs
%{
params.a_rate(750:1050) = 0.3;
params.a_rate(850:950) = 0.01;
params.e_rate(800:1000) = 20;
params.e_rate(850:950) = 500;
params.r_rate(750:900) = 1;
params.l_rate(750:800) = 20;
params.l_rate(800:900) = 100;
params.r_rate(901:951) = 100;
params.r_rate(951:1050) = 20;
params.l_rate(901:1050) = 1;
%}
params.a_rate(750:1050) = 0;
params.r_rate(750:1050) = 0;
params.l_rate(750:1050) = 0;


%%

[time, s_hist] = gillespie(params, 'n_steps', n_steps, 's0', zeros(1,gen_len));

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
