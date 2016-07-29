%%%%%%%%%%%
%
% In this example, I will show a simple implementation of the Gillespie
% algorithm. The example we will work on is a birth/death, or creation
% destruction, process. Two things can happen:
% 1 - the creation of a molecule/particle
% 2 - the destruction of a molecule/particle
% We will only consider one species of particles, and they will be
% independent of each other. We will monitor the overall count of
% particles, which we will call nn. (I am using two n, because otherwise
% when you want to find the variable in the code using the 'find' function,
% you will have a lot of trouble.)
%
% The general idea is that we first define the two 'reaction channels'.
% This means we define the events that will take place in the two cases of
% creation or destruction, which I will call 'reactions'. Each of these
% reactions has an associated 'reaction channel'. When the reaction channel
% is fired, a specific change happens in the system. The Gillespie
% simulation consists of stochastically determined series of individual
% firing events, where always exactly one reaction channel is executed
% exactly once. After each of these reactions, the 'propensities' of each
% channel to fire are calculated, to tell us how likely it is that for each
% reaction channel to be fired in the next reaction, and how long we will
% have to wait till then.
%
%%%%%%%%%%%
 
%%%% Definition of model parameters
% In the model we are using, two events can take place - the creation and
% the destruction of a molecule. The first one takes place at a constant
% rate:
 
kk = 1; % Constant rate of molecule creation
 
% The second one, destruction, happens at a rate that is proportional to
% the number of molecules that currently exist. This rate will later be
% multiplied with the number nn.
 
gamma = 0.1; % Destruction rate
 
%%%% Definition of general simulation parameters
% As a first thing, we have to set how many individual reactions should
% maximally be simulated
 
R_max = 1000; % maximal simulation steps
 
%%%% Preallocation of containers for the simulation time course
% To make the code more efficient, we 'pre-allocate' the variabels we will
% use to store the molecule counts and time points of our Gillespie
% simulation. This will make the code run faster - when MatLab has to
% change the size (number of elements) of a variable, what it will do is
% to destroy the old variable, make a new one in the right size and work
% with it. This is slow, and when it is done for every reaction step, we
% will make a very slow program. The way around it is to reserve enough
% space in a variable beforehand, and then save all the molecule counts and
% reaction times into there.
 
nn_cont = zeros(1,R_max+1); % Container for nn time course
tt_cont = zeros(1,R_max+1); % Container for time points of reactions
 
%%%% Pre-drawing of random numbers for simulation
% Very much like the pre-allocation of variables, we draw all the random
% numbers we will need beforehand. MatLab is fast at making large sets of
% numbers at once, while calling for these activities is very slow. So,
% what we want to do is call for the activity once, and have all the random
% numbers ready to work with later.
 
rand_nums = rand(2,R_max); % TWO random numbers for each reaction step
 
%%%% Definition of state changes for each reaction channel
% The next thing to define is the change in molecule number that occurs
% when a specific reaction channel is fired. We have two reaction channels,
% so we need two state changes, one for each reaction channel.
 
state_changes = [+1,-1];
 
% The first one is for the 'creation' reaction channel, where one molecule
% is added. The second one is for the 'destruction' channel, where one
% molecule is removed.
 
%%%% The actual simluation
% To simulate the time course, we now go step-by-step from reaction to
% reaction, until we have reached the maximal number of reactions. As a
% first thing, we initialize our simulation at the time point 0
tt = 0; % current time
nn = 0; % current molecule count
 
% We also need a counter that tells us how many reactions have occured
% until now.
R_counter = 0; % reaction counter
 
% Now we can run a while loop, which is repeated until we reach the maximum
% number of reactions
while R_counter < R_max
 
R_counter = R_counter + 1; % Tell the counter one more reaction happened
 
% As a first thing, we need to calculate the propensities for both of
% the reaction channels to fire. Let us first allocate them, and then
% assign them, for clarity's sake
aa = zeros(1,2); % Pre-allocated propensity vector, no values yet
 
% Now, we assign the propensity for a molecule to be created. This
% happens at a constant rate, so that the propensity is also constant.
aa(1) = kk;
 
% Next, we assign the propensity for a molecule to be destroyed. This
% happens at a rate that is proportional to the number of currently
% existent molecules - two molecules are double as likely to lead to
% the destruction of a molecule than one molecule is.
aa(2) = gamma.*nn;
 
% Now we draw the waiting time till the next reaction channel is
% fired. For this we use one of the pre-drawn random numbers
deltat_t = -log(rand_nums(1,R_counter))./sum(aa);
 
% We also need to find out which of the reaction channels is fired at
% that time. For this we have to choose from the possible reaction
% channels, according to their relative propensity. The idea here is to
% construct an interval that ranges from 0 to 1, [0,1], which contains
% different segments that correspond to the different reaction
% channels. Then we draw a random number between 0 and 1, and dependent
% on which segments it falls into, we choose the reaction channel that
% is fired. It would look a little bit like this:
%
%   rand_num=0.5            x
%             0[-----------------------,---]1 <- interval with segments
%   channel               1              2
%
% Here the random number for choosing the reaction channel is roughly
% 0.5, which means we choose the weighted segment that covers the
% point 0.5 in the interval. Here, that would refer to reaction channel
% 1, so the next channel that would be fired would be channel 1.
%
% You might realize that channels that have a length of 0 could mess up
% this procedure, so we have to remove them from the following
% computation. It looks a little complicated, but it will avoid a lot
% of trouble in simulations that are more complicated than this simple
% example.
 
% So let us first get only the info from the valid reaction channels,
% that is the ones with >0 propensity.
 
valid_inds = aa > 0; % Find the indices that refer to valid reactions
valid_aa = aa(valid_inds); % Use only valid reactions
valid_changes = state_changes(valid_inds); % Use only valid changes
 
% Next, we can construct the intervals described above
 
selection_intvl = cumsum(valid_aa); % Cumulative sum
% Then, we have to normalize the interval to a [0,1]. We do this by
% division by the last element.
selection_intvl = selection_intvl./selection_intvl(end);
 
% What is left to do, is to pick the segment as described above, and
% execute the according reaction channel's change. The values in
% selection_intvl refer to the upper ends of the segments of the
% interval. So, what we really need to find is the segment, which has
% the smallest upper end that is greater than the random number.
 
% The following operation is simply an implementation of what I
% described above. It will probably not be immediately clear. Go to the
% the MatLab command line, and type 'help find', that should explain to
% you how this works in detail.
selected_ind = find(selection_intvl>rand_nums(2,R_counter),1,'first');
selected_change = valid_changes(selected_ind);
 
% To finish the loop, we apply the time and state change we
% stochastically determined above, and store the current state of the
% simulation into the pre-allocated storage variables.
 
tt = tt + deltat_t;
nn = nn + selected_change;
 
tt_cont(R_counter+1) = tt;
nn_cont(R_counter+1) = nn;
 
end
 
% Let us see what the outcome of this is
plot(tt_cont,nn_cont,'k-')
xlabel('t')
ylabel('n')