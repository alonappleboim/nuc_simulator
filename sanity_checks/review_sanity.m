load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\wt_centers.mat')
load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\sequences_structure.mat')

feat = zeros(1,80);

for i = 1:80

    try
        load(['C:\Users\Daniel\Documents\MATLAB\nuc_simulator\sanity_checks\results_slide\sanity_results_' num2str(i)])
    catch a
        continue
    end

    feat(i) = max(max(features));
    
    figure
    [C, h] = contour(features,7, 'Fill', 'on');
    clabel(C, h);
    xlabel('RSC Length Parameter')
    ylabel('RSC Slide Intensity Parameter')
    legend('Likelihood Test Result')
    title(['Likelihood Result as a function of RSC Length and Eviction Intensity - gene ' num2str(i)])
    
end


%plot(feat, '^r')