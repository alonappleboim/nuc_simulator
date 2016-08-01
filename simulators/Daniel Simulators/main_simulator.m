
%% %%% PARAMETERS %%%%%

n_steps = 40000;
params = struct();
gen_len = 4000;
params.nuc_width = 147;
params.slide_len = 20;
params.a_rate = ones(1,gen_len); 
params.e_rate = ones(1,gen_len);
params.r_rate = 0.1*ones(1,gen_len); 
params.l_rate = 0.1*ones(1,gen_len); 
params.nuc_footprint = ones(1,(params.nuc_width.*2) - 1);
params.linker_len = 1;

% NFRs
params.a_rate(750:1050) = 0.1;
params.e_rate(750:1050) = 5;
params.e_rate(850:950) = 50;
params.r_rate(900:1050) = 20;
params.l_rate(750:900) = 20;

%%

[time, s_hist] = gillespie_sim(params, 'n_steps', n_steps, 's0', zeros(1,gen_len));

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
