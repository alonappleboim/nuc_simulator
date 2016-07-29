% SIMULATOR - allowing input, output and jumping of nucleosomes %

%%%%% PARAMETERS %%%%%

gen_len = 4000;                % the length of the genome
nuc_len = 147;                 % the length of a single nucleosome (odd number for convenience)
jump_len = 20; 				   % the length of a single jump (left or right) - even for convenience
nuc_vec = ones(1,(nuc_len*2)-1); % the vector for the in_vec convulution
jump_vec = ones(1,nuc_len+jump_len); % the vector for the jumping convulutions
i_prob = ones(1,gen_len);     % a probability vector for the input of a nucleosome into each bp
o_prob = ones(1,gen_len);     % a probability vector for the output of a nucleosome from each bp
r_prob = ones(1,gen_len);     % a probability vector for the right jumps from each bp
l_prob = ones(1,gen_len);     % a probability vector for the left jumps from each bp
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
    
	% generate the input, output and jump vectors:
    out_vec = state; % which bp can have a nuc extracted from them
	
	in_vec = conv(state,nuc_vec,'same'); % which bp can have a nuc inserted to them
	in_vec(in_vec~=0) = 1;
	in_vec = ~in_vec;
	
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
    %for j=1:(fix(nuc_len/2))
    %    in_vec(j)=0;
    %    in_vec(end-j+1)=0;
    %end
        
	% generate the change_matrix:
	change_matrix = zeros(gen_len,5);
	change_matrix(:,2) = out_vec .* o_prob; % the chances of output
	change_matrix(:,3) = in_vec .* i_prob;  % the chances of input
	change_matrix(:,4) = left_vec .* state .* l_prob; % the chances of jumping left
	change_matrix(:,5) = right_vec .* state .* r_prob; % the chances of jumping right
	change_matrix(:,1) = change_matrix(:,2) + change_matrix(:,3) + change_matrix(:,4) + change_matrix(:,5);
	change_matrix(:,1) = change_matrix(:,1) ./ sum(change_matrix(:,1)); % total chances, normalized
	
	% choose a bp that will be changed and decide what changes:
	bp = find(mnrnd(1,change_matrix(:,1)));
	bp_line = change_matrix(bp,2:end);
	bp_line = bp_line ./ sum(bp_line);
	change = find(mnrnd(1,bp_line)); % 1 is out, 2 is in, 3 is left, 4 is right
	
	% make the change in the state:
	if (change == 1)
		state(bp) = 0;
	end
	if (change == 2)
		state(bp) = 1;
	end
	if (change == 3)
		state(bp) = 0;
		state(bp-jump_len) = 1;
	end
	if (change == 4)
		state(bp) = 0;
		state(bp+jump_len) = 1;
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
