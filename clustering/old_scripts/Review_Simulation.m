load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\wt_centers.mat')
load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\sequences_structure.mat')

genlen = 3500;
TSS = fix(genlen/2);
NFR_pos = [TSS-299 : TSS+150];

gene_id = 200;

create_full_params_301016;

seq = sequences_structure(gene_id,:);
wt =  wt_3h(gene_id,:);

% make the wt data the right length:
buffer = genlen - 2501;
right_buffer = fix((buffer-500)/2);
left_buffer = right_buffer + 500;
if (right_buffer + left_buffer < buffer)
    left_buffer = left_buffer + 1;
end
wt = [zeros(1,left_buffer), wt, zeros(1,right_buffer)];

[ PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites ] = ...
    Extract_Sites_From_Gene(seq, genlen);

nuc_sum = nuc_sums(best_sim_index, :);
coverage = ksdensity(1:length(nuc_sum),1:length(nuc_sum),'weights',double(nuc_sum(1:end)),'width',5);
smoothed_wt = ksdensity([1:length(wt)],[1:length(wt)],'weights',double(wt),'width',5);

% full gene plot:
figure;
plot(smoothed_wt,'b')
hold on
plot(coverage,'r')
plot(PolyA_Sites(:) .* mean(smoothed_wt),'k')
plot(PolyT_Sites(:) .* mean(smoothed_wt),'m')
plot(REB1_Sites(:) .* 4 .* mean(smoothed_wt), 'g')
plot(ABF1_Sites(:) .* 4 .* mean(smoothed_wt), 'c')
plot(RAP1_Sites(:) .* 4 .* mean(smoothed_wt), 'y')
legend('wild-type','simulation','PolyA (right)','PolyT (left)', 'REB1', 'ABF1', 'RAP1')
xlabel(['Position (TSS at ' num2str(fix(genlen/2)) ')'])
ylabel('Intensity')
title(['Likelihood = ' num2str(features(best_sim_index)) char(10) 'Gene Index = ' num2str(gene_id)])

% make the textbox string:
textbox_string = [...
    'Poly Rate = ' num2str(params(1, best_sim_index)) char(10) ...
    'TF Eviction Effect = ' num2str(params(3, best_sim_index)) char(10) ...
    'RSC Eviction Effect = ' num2str(params(5, best_sim_index)) char(10) ...
    'RSC Eviction Length = ' num2str(params(4, best_sim_index)) char(10) ...
    'RSC Slide Effect = ' num2str(params(6, best_sim_index)) char(10) ...
    'RSC Slide Length = ' num2str(params(4, best_sim_index) * 2) ...
    ];

data = conv(wt(NFR_pos), gausswin(10)./sum(gausswin(10)), 'same');

% feature plot (only NFR):
nuc_sum = conv(nuc_sum, gausswin(10)./sum(gausswin(10)), 'same');
nuc_sum = nuc_sum(NFR_pos) .* sum(data) ./ sum(nuc_sum(NFR_pos));
figure;
plot(data, 'b')
hold on
plot(nuc_sum, 'r')
legend('wild-type','simulation')
xlabel('Position (TSS at 400)')
ylabel('Intensity')
title(['Likelihood = ' num2str(features(best_sim_index)) char(10) 'Gene Index = ' num2str(gene_id)])
annotation('textbox','String',textbox_string,'FitBoxToText','on')