
%% %%% PARAMETERS %%%%%

n_steps = 15000;
params = struct();
gen_len = 6000;
params.nuc_width = 147;
params.slide_len = 10;
params.a_rate = ones(1,gen_len); 
params.e_rate = ones(1,gen_len);
params.r_rate = 0.1*ones(1,gen_len); 
params.l_rate = 0.1*ones(1,gen_len); 
params.nuc_footprint = ones(1,(params.nuc_width.*2) - 1);
params.linker_len = 5;

% NFRs
params.r_rate(1050:1200) = 30;
params.l_rate(1000:1050) = 30;

% Polymerase Parameters
params.l_rate(1200:3200) = 5;

%%

[time, s_hist] = gillespie(params, 'n_steps', n_steps, 's0', zeros(1,gen_len));

%% %%% FINAL CALCULATIONS %%%%%

% get the number of times each bp had a nuceosome center on it
centers_vector = sum(s_hist(:,161:2660));
%smooth_vector =  conv(centers_vector, ones(1,51), 'same');
smooth_vector = ksdensity(1000:length(centers_vector)-1,1000:length(centers_vector)-1,'weights',double(centers_vector(1000:end-1)),'width',20);

% rescale graphs
%smooth_vector = smooth_vector .* (max(centers_vector)/mean(smooth_vector));

%% %%% OUTPUT GRAPHS %%%%%

%{
figure;
plot(centers_vector)
title('Number of Nucleosome Centers VS Base Pair')
%}

figure;
plot(smooth_vector,'g')
title('Number of Nucleosome Coverage VS Base Pair')

%{
figure;
plot(centers_vector,'r')
hold on
plot(smooth_vector,'g')
legend('centers','coverage')
title('Combined')
xlabel('Position')
ylabel('Time')
hold off
%}