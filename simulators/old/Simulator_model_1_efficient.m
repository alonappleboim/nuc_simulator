% SIMULATOR - only input/output %

%%%%% PARAMETERS %%%%%

gen_len = 4000;                % the length of the genome
nuc_len = 147;                 % the length of a single nucleosome (odd number for convenience)
nuc_vec = ones(1,(nuc_len*2)-1); % the vector for the in_vec convulution
i_prob = ones(1,gen_len);     % a probability vector for the input of a nucleosome into each bp
o_prob = ones(1,gen_len);     % a probability vector for the output of a nucleosome from each bp
num_steps = 20000;              % the number of steps to simulate
state = zeros(1,gen_len);     % the initial state

state_matrix = zeros(num_steps+1,gen_len); % the matrix keeping track of all the states
state_matrix(1,:) = state; 

% change the input probabilities so that we have an NFR in the middle:
for i=1800:2200
    i_prob(i) = 0;
end

%%%%% BIG LOOP %%%%%

for i=1:num_steps
    
    % check progress:
    if(mod(i,1000)==0)
        i
    end
    
	% generate the input and output vectors:
    out_vec = state; % which bp can have a nuc extracted from them
	in_vec = conv(state,nuc_vec,'same'); % which bp can have a nuc inserted to them
	in_vec(in_vec~=0) = 1;
	in_vec = ~in_vec;
    
    % make sure the edges are NFRs:
    %for j=1:(fix(nuc_len/2))
    %    in_vec(j)=0;
    %    in_vec(end-j+1)=0;
    %end
        
	% generate the change_matrix:
	change_matrix = zeros(gen_len,3);
	change_matrix(:,2) = out_vec .* o_prob; % the chances of output
	change_matrix(:,3) = in_vec .* i_prob;  % the chances of input
	change_matrix(:,1) = change_matrix(:,2) + change_matrix(:,3);
	change_matrix(:,1) = change_matrix(:,1) ./ sum(change_matrix(:,1)); % total chances, normalized
	
	% choose a bp that will be changed and decide what changes:
	bp = find(mnrnd(1,change_matrix(:,1)));
	bp_line = change_matrix(bp,2:end);
	bp_line = bp_line ./ sum(bp_line);
	change = find(mnrnd(1,bp_line)); % 1 is out, 2 is in
	
	% make the change in the state:
	if (change == 1)
		state(bp) = 0;
	end
	if (change == 2)
		state(bp) = 1;
	end
	state_matrix(i+1,:) = state;
	
end

%%%%% FINAL CALCULATIONS %%%%%

% get the number of times each bp had a nuceosome center on it
centers_vector = sum(state_matrix(:,:));
smooth_vector =  conv(centers_vector, ones(1,nuc_len), 'same');

% make the two graphs be of the same order of magnitude:
smooth_vector_smaller = smooth_vector .* (max(centers_vector)/mean(smooth_vector));

%%%%% OUTPUT GRAPHS %%%%%

figure;
plot(centers_vector)
title('Number of Nucleosome Centers VS Base Pair')
figure;
plot(smooth_vector)
title('Number of Nucleosome Coverage VS Base Pair')
figure;
plot(centers_vector,'r')
hold on
plot(smooth_vector_smaller,'g')
legend('centers','smooth')
title('Combined')
hold off
