create_full_params;

wt_length = zeros(1,100);
wt_evic = zeros(1,100);
wt_slide = zeros(1,100);
sth1_length = zeros(1,100);
sth1_evic = zeros(1,100);
sth1_slide = zeros(1,100);

for i = 1:100
    try
        load(['C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\results\wt_101016\results_' num2str(i) '.mat'])
    catch a
        wt_length(i) = nan;
        wt_evic(i) = nan;
        wt_slide(i) = nan;
        continue
    end
    parameters = params(: , best_sim_index);
    wt_length(i) = parameters(4);
    wt_evic(i) = parameters(5);
    wt_slide(i) = parameters(6);
    
    try
        load(['C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\results\sth1_6h_201016\results_' num2str(i) '.mat'])
    catch a
        sth1_length(i) = nan;
        sth1_evic(i) = nan;
        sth1_slide(i) = nan;
        continue
    end
    parameters = params(: , best_sim_index);
    sth1_length(i) = parameters(4);
    sth1_evic(i) = parameters(5);
    sth1_slide(i) = parameters(6);
end

figure
subplot(3,1,1)
plot(wt_evic,'^r')
hold on
plot(sth1_evic,'ob')
title('Wild Type VS STH1 6h')
ylabel('Eviction Intensity')
xlabel('Gene ID')
legend('wild type','sth1 6h')
subplot(3,1,2)
plot(wt_slide,'^r')
hold on
plot(sth1_slide,'ob')
ylabel('Slide Intensity')
xlabel('Gene ID')
legend('wild type','sth1 6h')
subplot(3,1,3)
plot(wt_length,'^r')
hold on
plot(sth1_length,'ob')
ylabel('RSC Length')
xlabel('Gene ID')
legend('wild type','sth1 6h')

