% assume batch index (bi) and size (BS) are set in environment
addpath(genpath('~/Dropbox/MATLAB/'))
addpath(genpath('~/Dropbox/ChIP_Analysis/Scripts/'))
disp('added paths to matlab instance')

%% data prep
Nucs = normalize('~/Dropbox/ChIP_Analysis/Files/Nucs_19jan2014.mat',...
                 {'log', 'input_clust','zscore','time'},...
                 {1, -1, 1, 1},...
                 'array');
[N,M,T] = size(Nucs.Data);
disp('data loaded')

%% do work
output = struct();
output.verbosity = 1;
fopts = struct();
fopts.n_trials = 20;
%output.plot = true;
%output.prompt = false;
starti = (bi - 1) * BS + 1
endi = min(bi * BS, N)
params = struct();

disp('starting loop...')

fit = interp_pulsefit(Nucs.Data(starti:endi,:,:), Nucs.timepoints, output);

for i = 1:length(fit.param_names)
    fit.(fit.param_names{i}) = squeeze(fit.params(:,:,i));
end
fit = rmfield(fit, 'params');
fit = rmfield(fit, 'param_names');
fit.nuc_inds = starti:endi;

fname = sprintf('~/Dropbox/MATLAB/chip_analysis/dynamics/fit_results/nuc_fit_batch_%i', bi);
disp(['writing to: ',fname]);
save(fname,'fit')
