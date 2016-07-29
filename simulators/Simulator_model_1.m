% SIMULATOR - only input/output %

%%%%% PARAMETERS %%%%%

gen_len = 2000;                % the length of the genome
nuc_len = 147;                 % the length of a single nucleosome (odd number for convenience)
i_prob = ones(1,gen_len);     % a probability vector for the input of a nucleosome into each bp
o_prob = ones(1,gen_len);     % a probability vector for the output of a nucleosome from each bp
num_steps = 20000;              % the number of steps to simulate
state = zeros(1,gen_len);     % the initial state
state_matrix = zeros(num_steps+1,gen_len); % the matrix keeping track of all the states
state_matrix(1,:) = state; 

% change the input probabilities so that we have an NFR in the middle:
for i=800:1200
    i_prob(i) = 0;
end

%%%%% BIG LOOP %%%%%

for i=1:num_steps
    
    % check progress:
    if(mod(i,1000)==0)
        i
    end
    
    in_vec = ones(1,gen_len); % which bp can have a nuc inserted to them
    out_vec = zeros(1,gen_len); % which bp can have a nuc extracted from them
    
    % generate the input and output vectors:
    %for j=1:(fix(nuc_len/2)) %% this part decides if the edges are NFRs %%
    %    in_vec(j)=0;
    %    in_vec(end-j+1)=0;
    %end
    for j=1:gen_len
        if (state(j) == 1)
            out_vec(j) = 1;
            for l=max(1,j-(2*fix(nuc_len/2))):min(gen_len,j+(2*fix(nuc_len/2))) % cancel possible insertions
                in_vec(l) = 0;
            end
        end
    end
        
    % decide which step is taken, using a change_matrix, where each row is
    % a possible position change, that is represented by three values: 
    % [the chance of the change, the index of the nucleosome, in or out]
    k = 2;
    change_matrix = [0,0,0]; 
    for j=1:gen_len % make the change_matrix in this loop
        if (in_vec(j) == 1)
            change_matrix(k,:) = [change_matrix(k-1,1)+i_prob(j),j,1];
            k=k+1;
        end
        if (out_vec(j) == 1)
            change_matrix(k,:) = [change_matrix(k-1,1)+o_prob(j),j,0];
            k=k+1;
        end
    end
    
    % randomly choose a position to switch to - position_change:
    total = change_matrix(end,1);
    random_number = rand * total;
    for j=size(change_matrix,1):-1:1
        if (random_number > change_matrix(j,1))
            position_change = change_matrix(j+1,:);
            break
        end
    end
    
    % change the state!
    state(position_change(2)) = position_change(3);
    state_matrix(i,:) = state;
end

%%%%% FINAL CALCULATIONS %%%%%

% get the number of times each bp had a nuceosome center on it
centers_vector = zeros(1,gen_len);
smooth_vector =  zeros(1,gen_len);
for i=1:num_steps
    for j=1:gen_len
        if(state_matrix(i,j) == 1)
            centers_vector(j) = centers_vector(j)+1;
            for l=max(1,j-fix(nuc_len/2)):min(gen_len,j+fix(nuc_len/2))
                smooth_vector(l) = smooth_vector(l)+1;
            end
        end
    end
end

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
