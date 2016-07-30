% SIMULATOR - Now With GILLESPIE For Extra Strength %
% this file is split to cells with this "%%" - each can be executed
% separately, which is sometimes very convinient...

%% clear all?
clear;
clc;
close all;


%% %%% PARAMETERS %%%%%

REP_EVERY = 1000;
gen_len = 4000;                % the length of the genome
nuc_len = 147;                 % the length of a single nucleosome (odd number for convenience)
jump_len = 20; 				   % the length of a single jump (left or right) - even for convenience
nuc_vec = ones(1,(nuc_len*2)-1); % the vector for the in_vec convulution
jump_vec = ones(1,nuc_len+jump_len); % the vector for the jumping convulutions
a_rate = ones(1,gen_len);     % a rate vector for the assembly of a nucleosome into each bp
e_rate = ones(1,gen_len);     % a rate vector for the evicition of a nucleosome from each bp
r_rate = ones(1,gen_len);     % a rate vector for the right slide from each bp
l_rate = zeros(1,gen_len);     % a rate vector for the left slide from each bp
num_steps = 5000;              % the number of steps to simulate
rand_nums = rand(1,num_steps);  % random numbers for calculating the dt of each change (GILLESPIE)
state = zeros(1,gen_len);     % the initial state

state_hist = zeros(num_steps+1,gen_len); % the state history
time = zeros(num_steps,1); % keeping track of simulation time
state_hist(1,:) = state; 

% change the input probabilities so that we have an NFR in the middle:
nfr_hw = 200;
a_rate((gen_len/2-nfr_hw):(gen_len/2)+nfr_hw) = 0;

%% %%% MAIN LOOP %%%%%

for i=1:num_steps
    
    % report progress:
    if(mod(i,REP_EVERY)==0)
        clc;
        fprintf('At iteration %i (%d%%)...\n', i, 100*i./num_steps);
    end;
    
	% generate the assembly, eviction and slide vectors:
    evic = state; % which bp can have a nuc evicted from them
	
	assem = conv(state,nuc_vec,'same'); % which bp can have a nuc inserted to them
	assem(assem~=0) = 1;
	assem = ~assem;
	
	temp_jump_vec = conv(state,jump_vec,'same');
	temp_jump_vec(temp_jump_vec ~= 0) = 1;
    temp_jump_vec = ~temp_jump_vec; % ones are the positions that indicate that there is enough space for a jump
	
	temp_right_vec = find(temp_jump_vec) - (jump_len/2) - (2*fix(nuc_len/2)); % find the bps that can jump right
	temp_left_vec = find(temp_jump_vec) + (jump_len/2) + (2*fix(nuc_len/2)); % find the bps that can jump left
	
	right_vec = zeros(1,gen_len); % which bp have a nuc that can jump right
    temp = temp_right_vec > 0;
	right_vec(temp_right_vec(temp)) = 1; % remove bps that are left of the edge (can't jump right)
	left_vec = zeros(1,gen_len); % which bp have a nuc that can jump left
    temp = temp_left_vec <= gen_len;
	left_vec(temp_left_vec(temp)) = 1; % remove bps that are right of the edge (can't jump left)
    
    % make sure the edges are NFRs:
    %{
    for j=1:(fix(nuc_len/2))
        in_vec(j)=0;
        in_vec(end-j+1)=0;
    end
    %}
    
	% generate the change_matrix:
	ep = evic .* e_rate; % probability of eviction
	ap = assem .* a_rate;  % .. of assembly
	lp = left_vec .* state .* l_rate; % .. sliding left 
	rp = right_vec .* state .* r_rate; % .. and sliding right
	p = ep + ap + lp + rp;
    sp = sum(p);
    p = p./sp; % normalize p
	
	% find how much time passed until the next change (simulating exponential distribution):
	dt = -log(rand_nums(i))./sp;
	
	% sample a bp for a change, and sample a state change:
	bp = find(mnrnd(1,p));
	change = find(mnrnd(1,[ep(bp),ap(bp),lp(bp),rp(bp)]));
	
	% make the change in the state:
	if (change == 2)
		state(bp) = 1;
    else
        state(bp) = 0;
        if (change == 3)
            state(bp-jump_len) = 1;
        end
        if (change == 4)
            state(bp+jump_len) = 1;
        end
	end
	state_hist(i+1,:) = state; % we multiply by dt to get the duration of the state
    time(i+1) = time(i) + dt;
end
clc;

%% %%% FINAL CALCULATIONS %%%%%

% get the number of times each bp had a nuceosome center on it
centers_vector = sum(state_hist(:,:));
smooth_vector =  conv(centers_vector, ones(1,nuc_len), 'same');

% make the two graphs be of the same order of magnitude:
smooth_vector_smaller = smooth_vector .* (max(centers_vector)/mean(smooth_vector));

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
plot(smooth_vector_smaller,'g')
legend('centers','smooth')
title('Combined')
xlabel('Position')
ylabel('Time')
hold off
