function [nuc_sum, time, nuc_s_hist, nuc_evics, REB1_s_hist, ABF1_s_hist, RAP1_s_hist] = ... 
        Simulate(model_params, varargin)
%
% gillespie simulation for nucleosome positioning, including trans factors
% and using genomic data.
%
% Arguments:
%  modelparams - a struct with the following fields:
%    nuc_width - width of a nucleosome
%    REB1_width - width of REB1 trans factor
%    ABF1_width - width of ABF1 trans factor
%    RAP1_width - width of RAP1 trans factor
%    slide_len - sliding step per unit time
%    linker_len - the average length of the linkers between the nucleosomes
%    nuc_a_rate - assembly rate per position along genome
%    nuc_e_rate - eviction rate per position along genome
%    nuc_l_rate - left slide rate per position along genome
%    nuc_r_rate - right slide rate per position along genome
%    REB1_a_rate - REB1 assembly rate
%    REB1_e_rate - REB1 eviction rate
%    ABF1_a_rate - ABF1 assembly rate
%    ABF1_e_rate - ABF1 eviction rate
%    RAP1_a_rate - RAP1 assembly rate
%    RAP1_e_rate - RAP1 eviction rate
%
% Name/Value arguments:
%    n_steps - number of simulation steps. default = genome size * 10.
%    s0 - state at t=0. default = random state.
%    report - update every <report> itrations. default = n_steps/10. If
%             'no', then no report is generated
%  
% the returned state history has 1s for nucs, 2s for REB1s, 3s for ABF1s 
% and 4s for RAP1s. In addidtion, we return a nuc_state_history, which only
% includes the nucleosomes, multiplyed by the time of each state.

% parse arguments
params = model_params;
params.genlen = length(params.nuc_a_rate);
defaults = struct('n_steps', round((params.genlen/params.nuc_width)*100),...
                  's0', rand_state(params.genlen, params.nuc_width), ...
                  'report', nan);
extra_inputs = parse_namevalue_pairs(defaults, varargin); %simulation parameters

% if "report" was added to the parameters - don't report:
if isnan(extra_inputs.report)
    extra_inputs.to_report = true;
    extra_inputs.report = extra_inputs.n_steps./20;
else
    extra_inputs.to_report = false;
end

% allocate memory:
nuc_s_hist = false(1 + extra_inputs.n_steps, params.genlen);
%%% REB1_s_hist = false(1 + extra_inputs.n_steps, params.genlen);
%%% ABF1_s_hist = false(1 + extra_inputs.n_steps, params.genlen);
%%% RAP1_s_hist = false(1 + extra_inputs.n_steps, params.genlen);
nuc_sum = zeros(1, params.genlen);
nuc_evics = zeros(1, params.genlen);

% prep simulation
rand_nums = rand(1, extra_inputs.n_steps);  % random numbers for calculating the dt of each change (GILLESPIE)
nuc_s_hist(1, extra_inputs.s0 == 1) = true;
time = zeros(extra_inputs.n_steps+1, 1); % keeping track of simulation time

% build the linker vectors of the nucleosomes:
x = -10 : 20/params.linker_len : 10;
y = fliplr(x);
linker_left_vector = 1 ./ (1 + exp(-x)); % a sigmoid to be added to the nuc_vector
linker_right_vector = 1 ./ (1 + exp(-y)); % a sigmoid to be added to the nuc_vector

% make the vectors for the nucleosome (sliding and assembly):
nuc_vector = ones(1,(params.nuc_width.*2) - 1); % a kernel that shows where a nucleosome cant assemble
assem_vector = [linker_left_vector,nuc_vector,linker_right_vector];
slide_vec = ones(1, params.nuc_width + params.slide_len); % the kernel for the sliding test
slide_right_vec = [linker_left_vector,slide_vec,zeros(1,params.linker_len*2 + 1)];
slide_left_vec = [zeros(1,params.linker_len*2 + 1),slide_vec,linker_right_vector];

% time loop
state = extra_inputs.s0;
nuc_state = zeros(1,params.genlen);
nuc_state(state == 1) = 1;
%%% REB1_state = zeros(1,params.genlen);
%%% ABF1_state = zeros(1,params.genlen);
%%% RAP1_state = zeros(1,params.genlen);
for i = 1:extra_inputs.n_steps
    % report progress:
    if ((extra_inputs.to_report == true) && (mod(i,extra_inputs.report) == 0))
        clc;
        fprintf('At iteration %i (%d%%)...\n', i, 100*i./extra_inputs.n_steps);
    end;
                
	% generate the nucleosome eviction rate vector
    evic_rate = nuc_state .* params.nuc_e_rate;
	
	% generate the nucleosome assembly rate vector
    free_nuc_dna = (1-min(1,conv(nuc_state, assem_vector, 'same')));
    %%% free_REB1_dna = (1-min(1,conv(REB1_state, ones(1,params.nuc_width+params.REB1_width), 'same')));
    %%% free_ABF1_dna = (1-min(1,conv(ABF1_state, ones(1,params.nuc_width+params.ABF1_width), 'same')));
    %%% free_RAP1_dna = (1-min(1,conv(RAP1_state, ones(1,params.nuc_width+params.RAP1_width), 'same')));
    free_dna = free_nuc_dna; %%% .* free_REB1_dna .* free_ABF1_dna .* free_RAP1_dna;
	assem_rate = free_dna.*params.nuc_a_rate; % assembly rate
	
	% generate the left and right sliding rate vectors
 	temp_right_slide_vector = 1-min(1,conv(nuc_state, slide_right_vec, 'same'));
    %%% temp_right_slide_vector = temp_right_slide_vector - min(1,conv(REB1_state, ones(1,params.REB1_width + params.slide_len + fix(params.nuc_width/2)), 'same'));
    %%% temp_right_slide_vector = temp_right_slide_vector - min(1,conv(ABF1_state, ones(1,params.ABF1_width + params.slide_len + fix(params.nuc_width/2)), 'same'));
    %%% temp_right_slide_vector = temp_right_slide_vector - min(1,conv(RAP1_state, ones(1,params.RAP1_width + params.slide_len + fix(params.nuc_width/2)), 'same'));
    temp_right_slide_vector(temp_right_slide_vector < 0) = 0;
    temp_left_slide_vector = 1-min(1,conv(nuc_state, slide_left_vec, 'same'));
    %%% temp_left_slide_vector = temp_left_slide_vector - min(1,conv(REB1_state, ones(1,params.REB1_width + params.slide_len + fix(params.nuc_width/2)), 'same'));
    %%% temp_left_slide_vector = temp_left_slide_vector - min(1,conv(ABF1_state, ones(1,params.ABF1_width + params.slide_len + fix(params.nuc_width/2)), 'same'));
    %%% temp_left_slide_vector = temp_left_slide_vector - min(1,conv(RAP1_state, ones(1,params.RAP1_width + params.slide_len + fix(params.nuc_width/2)), 'same'));
    temp_left_slide_vector(temp_left_slide_vector < 0) = 0;

	right_vec = circshift(temp_right_slide_vector, -(fix(params.slide_len/2))-(fix(params.nuc_width/2))-(2*params.linker_len), 2);
	right_vec(end - ((fix(params.slide_len/2))+(fix(params.nuc_width/2))+(2*params.linker_len)):end) = 0;
	left_vec = circshift(temp_left_slide_vector, (fix(params.slide_len/2))+(fix(params.nuc_width/2))+(2*params.linker_len), 2);
	left_vec(1:(fix(params.slide_len/2))+(fix(params.nuc_width/2))+(2*params.linker_len)) = 0;
    	
	right_rate = params.nuc_r_rate.*right_vec.*nuc_state;
	left_rate = params.nuc_l_rate.*left_vec.*nuc_state;
    
    % generate the trans factors eviction rates:
    %%% REB1_evic_rate = REB1_state .* params.REB1_e_rate;
    %%% ABF1_evic_rate = ABF1_state .* params.ABF1_e_rate;
    %%% RAP1_evic_rate = RAP1_state .* params.RAP1_e_rate;
    
    % generate the trans factors assembly rates:
    %%% free_REB1_dna = (1-min(1,conv(REB1_state, ones(1,2*params.REB1_width), 'same')));
    %%% free_ABF1_dna = (1-min(1,conv(ABF1_state, ones(1,2.*params.ABF1_width), 'same')));
    %%% free_RAP1_dna = (1-min(1,conv(RAP1_state, ones(1,2.*params.RAP1_width), 'same')));
    %%% free_dna = free_nuc_dna .* free_REB1_dna .* free_ABF1_dna .* free_RAP1_dna;
    %%% REB1_assem_rate = free_dna .* (1 - REB1_state) .* params.REB1_a_rate;
    %%% ABF1_assem_rate = free_dna .* (1 - ABF1_state) .* params.ABF1_a_rate;
    %%% RAP1_assem_rate = free_dna .* (1 - RAP1_state) .* params.RAP1_a_rate;
    		
    % sample next time point
    all_rates = evic_rate + assem_rate + left_rate + right_rate; %%% + ...
        %%% REB1_assem_rate + ABF1_assem_rate + RAP1_assem_rate + ...
        %%% REB1_evic_rate + ABF1_evic_rate + RAP1_evic_rate;
    dt = -log(rand_nums(i))./sum(all_rates);
    
    % sample a position, and a state change:
	all_rates = [evic_rate,assem_rate,left_rate,right_rate]; %%% , ...
        %%% REB1_evic_rate, REB1_assem_rate, ABF1_evic_rate, ...
        %%% ABF1_assem_rate, RAP1_evic_rate, RAP1_assem_rate];
    possible_rates = find(all_rates~=0);
    all_rates = all_rates(all_rates~=0)./sum(all_rates); % normalize all_rates
	ind = possible_rates(mnrnd(1,all_rates)==1);
    bp = mod(ind, params.genlen);
    change = (ind-bp)./params.genlen + 1;
    if (bp == 0)
        bp = params.genlen;
        change = change - 1;
    end
  
	% make the change in the state:
    if (change == 1) % nuc eviction
        nuc_state(bp) = 0;
        nuc_evics(bp) = nuc_evics(bp) + 1;
    end
    if (change == 2) % nuc assembly
		nuc_state(bp) = 1;
    end
    if (change == 3) % nuc left slide
        nuc_state(bp) = 0;
        nuc_state(bp - params.slide_len) = 1;
    end
    if (change == 4) % nuc right slide
        nuc_state(bp) = 0;
        nuc_state(bp + params.slide_len) = 1;
    end
    
    %{
    if (change == 5) % REB1 eviction
        REB1_state(bp) = 0;
    end
    if (change == 6) % REB1 assembly
        REB1_state(bp) = 1;
    end
    if (change == 7) % ABF1 eviction
        ABF1_state(bp) = 0;
    end
    if (change == 8) % ABF1 assembly
        ABF1_state(bp) = 1;
    end
    if (change == 9) % RAP1 eviction
        RAP1_state(bp) = 0;
    end
    if (change == 10) % RAP1 assembly
        RAP1_state(bp) = 1;
    end
    %}
    
    nuc_s_hist(i+1, nuc_state > 0) = true;
    %%% REB1_s_hist(i+1, REB1_state > 0) = true;
    %%% ABF1_s_hist(i+1, ABF1_state > 0) = true;
    %%% RAP1_s_hist(i+1, RAP1_state > 0) = true;
    time(i+1) = time(i) + dt;

    if (i > 5000) % start taking states into account only after a while
        nuc_sum = nuc_sum + (nuc_state .* dt);
    end
end

%%% for now, until we actually use the TFs in the simulations:
REB1_s_hist = 0; 
ABF1_s_hist = 0;
RAP1_s_hist = 0;
%%%

if extra_inputs.to_report, clc; end;
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