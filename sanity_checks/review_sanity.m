load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\wt_centers.mat')
load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\sequences_structure.mat')

for i = 50:75

    try
        load(['C:\Users\Daniel\Documents\MATLAB\nuc_simulator\sanity_checks\results\sanity_results_' num2str(i)])
    catch a
        continue
    end

    figure
    contour(features)
    xlabel('RSC Length Parameter')
    ylabel('RSC Eviction Intensity Parameter')
    legend('Likelihood Test Result')
    title(['Likelihood Result as a function of RSC Length and Eviction Intensity - gene ' num2str(i)])
end