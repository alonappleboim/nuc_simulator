%%

Gene_id = 23;
genlen = 3500;
TSS = round(genlen/2);
NFR_pos = [TSS-299 : TSS+150];

% Starting Variables:
FRS2_seq = sequences_structure(Gene_id,:);
FRS2_wt =  wt_3h(Gene_id,:);

% make the wt data the right length:
buffer = genlen - 2501;
right_buffer = fix((buffer-500)/2);
left_buffer = right_buffer + 500;
if (right_buffer + left_buffer < buffer)
    left_buffer = left_buffer + 1;
end
FRS2_wt = [zeros(1,left_buffer), FRS2_wt, zeros(1,right_buffer)];

[ PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites ] = ...
    Extract_Sites_From_Gene(FRS2_seq, genlen);

% run the simulation:
centers_vector = zeros(15,genlen);

for i=1:20
    [nuc_sum, time, nuc_s_hist, REB1_s_hist, ABF1_s_hist, RAP1_s_hist] = ...
                run_simulation_from_genome(FRS2_seq, 'linker_len', 10, 'n_steps', 10000, ...
                'poly_rate', 0, 'REB1_a_rate', 0.0001, 'REB1_e_rate', 0.0001, 'ABF1_a_rate', 0.0001, ...
                            'ABF1_e_rate', 0.0001, 'RAP1_a_rate', 0.0001, 'RAP1_e_rate', 0.0001,...
                            'TF_evic_intensity', 0.0001, 'RSC_evic_intensity', 0.0005, ...
                            'RSC_evic_length', 10, 'RSC_slide_length', 20, ...
                            'RSC_slide_intensity', 0.001, 'gen_len', genlen, 'slide_len', 3);
    centers_vector(15,:) = centers_vector(1,:) + nuc_sum;
end

for j = 1:14
    for i = 1:20
        [nuc_sum, time, nuc_s_hist, REB1_s_hist, ABF1_s_hist, RAP1_s_hist] = ...
            run_simulation_from_genome(FRS2_seq, 'linker_len', 10, 'n_steps', 10000, ...
            'poly_rate', 0, 'REB1_a_rate', 0.0001, 'REB1_e_rate', 0.0001, 'ABF1_a_rate', 0.0001, ...
                        'ABF1_e_rate', 0.0001, 'RAP1_a_rate', 0.0001, 'RAP1_e_rate', 0.0001,...
                        'TF_evic_intensity', 0.0001, 'RSC_evic_intensity', 0.05, ...
                        'RSC_evic_length', 10*j, 'RSC_slide_length', 20*j, ...
                        'RSC_slide_intensity', 0.6, 'gen_len', genlen, 'slide_len', 3);

        centers_vector(j,:) = centers_vector(j,:) + nuc_sum;
    end
end

%%

likelihood = zeros(1,15);
plus1 = zeros(1,15);
minus1 = zeros(1,15);
for j=1:15
    [likelihood(j), plus1(j), minus1(j)] = Compare_Sum_To_Data(centers_vector(j,:), FRS2_wt, NFR_pos, true);
end

std_plus1 = std(plus1)
std_minus1 = std(minus1)
std_likelihood = std(likelihood)

%%

figure
plot(plus1,'^b')
hold on
plot(minus1,'or')
plot((-likelihood ./ 5),'xk')
legend('+1','-1','Likelihood / (-5)')
xlabel('RSC Length Parameter (a.u.)')
ylabel('Feature Score')
title('Features as a Function of RSC Length Parameter (Gene 23)')

%%

smooth_wt = conv(FRS2_wt, gausswin(100), 'same');
smooth_wt = conv(smooth_wt, gausswin(20), 'same');

smooth_vector = zeros(15,genlen);
for j=1:15
    smooth_vector(j,:) = conv(centers_vector(j,:), gausswin(100), 'same');
    smooth_vector(j,:) = conv(smooth_vector(j,:), gausswin(20), 'same');
end

%%

figure
plot(FRS2_wt(NFR_pos), 'b')
hold on
plot(centers_vector(15,NFR_pos) .* (sum(FRS2_wt(NFR_pos)) ./ sum(centers_vector(15,NFR_pos))),'r')
legend('wild type', 'simulation')
xlabel('Position in NFR (bp)')
ylabel('Reads')
title('The Simulation Results When the Distance Between Peaks is ~20 (not smoothed)')

figure
plot(smooth_wt(NFR_pos), 'b')
hold on
plot(smooth_vector(15,NFR_pos) .* (sum(smooth_wt(NFR_pos)) ./ sum(smooth_vector(15,NFR_pos))),'r')
legend('wild type', 'simulation')
xlabel('Position in NFR (bp)')
ylabel('Reads')
title('The Simulation Results When the Distance Between Peaks is ~20 (smoothed)')

