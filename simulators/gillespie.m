function [t,s_hist] = gillespie(model_params, varargin)
%
% gillespie simulation for nucleosome positioning 
%
% Arguments:
%  modelparams - a struct with the following fields:
%    nuc_w - width of a nucleosome
%    slide_s - sliding step per unit time
%    a - assembly rate per position along genome
%    e - eviction rate per position along genome
%    l - left slide rate per position along genome
%    r - right slide rate per position along genome
%    nuc_fp - nucleosome footprint for assembly between 0..1. 0
%             means completely free, 1 means zero chance of assembly.
%
% Name/Value arguments:
%    n_steps - number of simulation steps. default = genome size * 10.
%    n_inst - number of instances to simulate. default = 10.
%    s0 - state at t=0. default = random state.
%    report - update every <report> itrations. default = n_steps/10. If
%             'no', then no report is generated
%  

<<<<<<< HEAD
% parse arguments
mp = model_params;
mp.genlen = length(mp.a);
defaults = struct('n_steps', round(mp.genlen/mp.nuc_w*10),...
                  'n_inst', 10,...
                  's0', rand_state(mp.genlen, mp.nuc_w), ...
                  'report', nan);
sp = parse_namevalue_pairs(defaults, varargin); %simulation parameters
if isnan(sp.report)
    sp.to_report = true;
    sp.report = sp.n_steps./10;
else
    if strcmp(sp.to_report, 'no'), sp.to_report = false;
    else sp.to_report = true; end
end

% prep simulation
slide_k = ones(1, mp.nuc_w+mp.slide_s); % the kernel for the sliding test
rand_nums = rand(1, sp.n_steps);  % random numbers for calculating the dt of each change (GILLESPIE)
=======
%% clear all?
clear;
clc;
close all;

%sdjkg
%% %%% PARAMETERS %%%%%

REP_EVERY = 1000;
gen_len = 8000;                % the length of the genome
nuc_len = 147;                 % the length of a single nucleosome (odd number for convenience)
jump_len = 20; 				   % the length of a single jump (left or right) - even for convenience
nuc_vec = ones(1,(nuc_len*2)-1); % the vector for the in_vec convulution
jump_vec = ones(1,nuc_len+jump_len); % the vector for the jumping convulutions
a_rate = 10*ones(1,gen_len);     % a rate vector for the assembly of a nucleosome into each bp
e_rate = ones(1,gen_len);     % a rate vector for the evicition of a nucleosome from each bp
r_rate = 10*ones(1,gen_len);     % a rate vector for the right slide from each bp
l_rate = 10*ones(1,gen_len);     % a rate vector for the left slide from each bp
num_steps = 20000;              % the number of steps to simulate
rand_nums = rand(1,num_steps);  % random numbers for calculating the dt of each change (GILLESPIE)
state = zeros(1,gen_len);     % the initial state

state_hist = zeros(num_steps+1,gen_len); % the state history
time = zeros(num_steps,1); % keeping track of simulation time
state_hist(1,:) = state; 

% change the input probabilities so that we have an NFR in the middle:
nfr_hw = 400;
a_rate((gen_len/2-nfr_hw):(gen_len/2)+nfr_hw) = 0;
l_rate((gen_len/2-nfr_hw):(gen_len/2)+nfr_hw) = 0;
r_rate((gen_len/2-nfr_hw):(gen_len/2)+nfr_hw) = 0;
>>>>>>> 912e43d72aa988d7a3ae92bb758b453b9e830560

s_hist = [sp.s0; zeros(sp.n_steps, mp.genlen)]; % the state history
t = zeros(sp.n_steps, 1); % keeping track of simulation time

% time loop
state = sp.s0;
for i = 1:sp.n_steps
    % report progress:
    if sp.to_report && mod(i,sp.report) == 0
        clc;
        fprintf('At iteration %i (%d%%)...\n', i, 100*i./sp.n_steps);
    end;
    
	% generate the assembly, eviction and sliding rate vectors
    er = state .* mp.e; % eviction rate
    free_dna = (1-min(1,conv(state, mp.nuc_fp, 'same')));
	ar = free_dna.*mp.a; %assembly rate
	lr = zeros(1,mp.genlen);
    rr = zeros(1,mp.genlen);
% 	can_slide = conv(state, slide_k, 'same');
% 	temp_jump_vec(temp_jump_vec ~= 0) = 1;
%     temp_jump_vec = ~temp_jump_vec; % ones are the positions that indicate that there is enough space for a jump
	
% 	temp_right_vec = find(temp_jump_vec) - (jump_len/2) - (2*fix(nuc_len/2)); % find the bps that can jump right
% 	temp_left_vec = find(temp_jump_vec) + (jump_len/2) + (2*fix(nuc_len/2)); % find the bps that can jump left
% 	
% 	right_vec = zeros(1,gen_len); % which bp have a nuc that can jump right
%     temp = temp_right_vec > 0;
% 	right_vec(temp_right_vec(temp)) = 1; % remove bps that are left of the edge (can't jump right)
% 	left_vec = zeros(1,gen_len); % which bp have a nuc that can jump left
%     temp = temp_left_vec <= gen_len;
% 	left_vec(temp_left_vec(temp)) = 1; % remove bps that are right of the edge (can't jump left)

    % sample next time point
    all_rates = er + ar + lr + rr;
    dt = -log(rand_nums(i))./sum(all_rates);
    
    % sample a position, and a state change:
	p = [er,ar,lr,rr];
    nzidx = find(p~=0);
    p = p(p~=0)./sum(p);
	ind = nzidx(mnrnd(1,p)==1);
    pos = mod(ind, mp.genlen);
    step = (ind-pos)./mp.genlen + 1;
	if pos == 0
        pos = mp.genlen;
        step = step - 1;
    end
    
<<<<<<< HEAD
=======
	% generate the change_matrix:
	ep = evic .* e_rate; % probability of eviction
	ap = assem .* a_rate;  % .. of assembly
	lp = left_vec .* state .* l_rate; % .. sliding left 
	rp = right_vec .* state .* r_rate; % .. and sliding right
	p = ep + ap + lp + rp;
    sp = sum(p);
    
    ep = ep./p; % normalize the vectors for mnrnd
    ap = ap./p;
    lp = lp./p;
    rp = rp./p;
    p = p./sp;
	
	% find how much time passed until the next change (simulating exponential distribution):
	dt = -log(rand_nums(i))./sp;
	
	% sample a bp for a change, and sample a state change:
	bp = find(mnrnd(1,p));
	change = find(mnrnd(1,[ep(bp),ap(bp),lp(bp),rp(bp)]));
	
>>>>>>> 912e43d72aa988d7a3ae92bb758b453b9e830560
	% make the change in the state:
	if (step == 2) %assembly
		state(pos) = 1;
    else %eviction
        state(pos) = 0;
        if (step == 3) %left slide
            state(pos-jump_len) = 1;
        end
        if (step == 4) %right slide
            state(pos+jump_len) = 1;
        end
	end
	s_hist(i+1,:) = state;
    t(i+1) = t(i) + dt;
end

if sp.to_report, clc; end;
end

function s = rand_state(N,w)
    s = rand(N,1) < 1./(2.*w+1);
    while true
        idx = find(s);
        r = find(diff(idx) < 2*w+1);
        if isempty(r), break; end;
        s(idx(r)) = 0;
    end
    s = 1.*s.';
end