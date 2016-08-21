function [time,state_history] = gillespie(model_params, varargin)
%
% gillespie simulation for nucleosome positioning 
%
% Arguments:
%  modelparams - a struct with the following fields:
%    nuc_width - width of a nucleosome
%    slide_len - sliding step per unit time
%    a_rate - assembly rate per position along genome
%    e_rate - eviction rate per position along genome
%    l_rate - left slide rate per position along genome
%    r_rate - right slide rate per position along genome
%    linker_len - the average length of the linkers between the nucleosomes
%
% Name/Value arguments:
%    n_steps - number of simulation steps. default = genome size * 10.
%    n_inst - number of instances to simulate. default = 10.
%    s0 - state at t=0. default = random state.
%    report - update every <report> itrations. default = n_steps/10. If
%             'no', then no report is generated
%  

% parse arguments
params = model_params;
params.genlen = length(params.a_rate);
defaults = struct('n_steps', round((params.genlen/params.nuc_width)*100),...
                  'n_inst', 10,...
                  's0', rand_state(params.genlen, params.nuc_width), ...
                  'report', nan);
extra_inputs = parse_namevalue_pairs(defaults, varargin); %simulation parameters
if isnan(extra_inputs.report)
    extra_inputs.to_report = true;
    extra_inputs.report = extra_inputs.n_steps./20;
else
    if strcmp(extra_inputs.to_report, 'no'), extra_inputs.to_report = false;
    else extra_inputs.to_report = true; end
end

% prep simulation
rand_nums = rand(1, extra_inputs.n_steps);  % random numbers for calculating the dt of each change (GILLESPIE)
state_history = [extra_inputs.s0; zeros(extra_inputs.n_steps, params.genlen)]; % the state history
time = zeros(extra_inputs.n_steps, 1); % keeping track of simulation time

x = -10 : 20/params.linker_len : 10;
y = fliplr(x);
linker_left_vector = 1 ./ (1 + exp(-x)); % a sigmoid to be added to the nuc_vector
linker_right_vector = 1 ./ (1 + exp(-y)); % a sigmoid to be added to the nuc_vector

nuc_vector = ones(1,(params.nuc_width.*2) - 1); % a kernel that shows where a nucleosome cant assemble
assem_vector = [linker_left_vector,nuc_vector,linker_right_vector];
slide_vec = ones(1, params.nuc_width + params.slide_len); % the kernel for the sliding test
%slide_vec = [linker_left_vector,slide_vec,linker_right_vector]; (ADD ZEROS INSTEAD OF THE IRRELIVENT)
slide_right_vec = [linker_left_vector,slide_vec,zeros(1,params.linker_len*2 + 1)];
slide_left_vec = [zeros(1,params.linker_len*2 + 1),slide_vec,linker_right_vector];

% time loop
state = extra_inputs.s0;
for i = 1:extra_inputs.n_steps
    % report progress:
    if extra_inputs.to_report && mod(i,extra_inputs.report) == 0
        clc;
        fprintf('At iteration %i (%d%%)...\n', i, 100*i./extra_inputs.n_steps);
    end;
    
	% generate the eviction rate vector
    evic_rate = state .* params.e_rate; % eviction rate
	
	% generate the assembly rate vector
    free_dna = (1-min(1,conv(state, assem_vector, 'same')));
	assem_rate = free_dna.*params.a_rate; % assembly rate
	
	% generate the left and right sliding rate vectors
 	temp_right_slide_vector = 1-min(1,conv(state, slide_right_vec, 'same'));
   	temp_left_slide_vector = 1-min(1,conv(state, slide_left_vec, 'same'));

	right_vec = circshift(temp_right_slide_vector, -(fix(params.slide_len/2))-(fix(params.nuc_width/2))-(2*params.linker_len), 2);
	right_vec(end - ((fix(params.slide_len/2))+(fix(params.nuc_width/2))+(2*params.linker_len)):end) = 0;
	left_vec = circshift(temp_left_slide_vector, (fix(params.slide_len/2))+(fix(params.nuc_width/2))+(2*params.linker_len), 2);
	left_vec(1:(fix(params.slide_len/2))+(fix(params.nuc_width/2))+(2*params.linker_len)) = 0;
	
	right_rate = params.r_rate.*right_vec.*state;
	left_rate = params.l_rate.*left_vec.*state;
    		
    % sample next time point
    all_rates = evic_rate + assem_rate + left_rate + right_rate;
    dt = -log(rand_nums(i))./sum(all_rates);
    
    % sample a position, and a state change:
	all_rates = [evic_rate,assem_rate,left_rate,right_rate];
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
    if (change == 2) %assembly
		state(bp) = 1;
    else %eviction
        state(bp) = 0;
        if (change == 3) %left slide
            state(bp - params.slide_len) = 1;
        end
        if (change == 4) %right slide
            state(bp + params.slide_len) = 1;
        end
    end
	state_history(i+1,:) = state.*dt;
    time(i+1) = time(i) + dt;
end

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