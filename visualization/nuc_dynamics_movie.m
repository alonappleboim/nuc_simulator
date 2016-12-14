function nuc_dynamics_movie(s_hist, t, path, gene_id, start_param, end_param, varargin)
% record a movie of the simulation for easy viewing.
% 
% Arguments:
%  s_hist - state_history
%  t - absoulte simulation time
%  path - to which mp4 video is saved.
%
% Name/Value Arguments:
%  sample_frame - amount of samples shown in each frame. default = 200.
%  frame_overlap - overlap between frames in %. default = .2 (20%).
%  frame_rate - frames per second. default = 12 fps.
%  frame_size - video frame size. default = [16, 8] (cm).
%  cmap - colormap. default = white->yellow->red->black
%

load('C:\Users\Daniel\Documents\MATLAB\Friedman Lab\Experiment Data\sequences_structure.mat')
create_params_sth1;
genlen = 3500;
TSS = 1750;
NFR_pos = [TSS-299:TSS+150];

defs = struct('sample_frame', 200, 'frame_overlap', .2, ...
              'frame_rate', 12, 'cmap', AdvancedColormap('w wy wyr wrk rrk'), ...
              'frame_size', [16,8]);
args = parse_namevalue_pairs(defs, varargin);

% make the vectors of the rates for before and after the sth1 damage:
[ PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites ] = ...
    Extract_Sites_From_Gene(sequences_structure(gene_id,:), genlen);
sim_params = params(: , start_param);
[ start_nuc_a_rate, start_nuc_e_rate, start_nuc_r_rate, start_nuc_l_rate ] = ...
    generate_rates_from_sites(PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites, ...
        'poly_rate',0, ...
        'REB1_a_rate', 0.0001, 'REB1_e_rate', 0.0001, ...
        'ABF1_a_rate', 0.0001, 'ABF1_e_rate', 0.0001, ...
        'RAP1_a_rate', 0.0001, 'RAP1_e_rate', 0.0001, ...
        'TF_evic_intensity', sim_params(1), ...
        'RSC_evic_length', sim_params(2), 'RSC_slide_length', sim_params(2).*sim_params(3), ...
        'RSC_evic_intensity', sim_params(4), ...
        'RSC_slide_intensity', sim_params(4)*sim_params(5));
sim_params = params(: , end_param);
[ end_nuc_a_rate, end_nuc_e_rate, end_nuc_r_rate, end_nuc_l_rate ] = ...
    generate_rates_from_sites(PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites, ...
        'poly_rate',0, ...
        'REB1_a_rate', 0.0001, 'REB1_e_rate', 0.0001, ...
        'ABF1_a_rate', 0.0001, 'ABF1_e_rate', 0.0001, ...
        'RAP1_a_rate', 0.0001, 'RAP1_e_rate', 0.0001, ...
        'TF_evic_intensity', sim_params(1), ...
        'RSC_evic_length', sim_params(2), 'RSC_slide_length', sim_params(2).*sim_params(3), ...
        'RSC_evic_intensity', sim_params(4), ...
        'RSC_slide_intensity', sim_params(4)*sim_params(5));

% get the data of the 0m and the 6h:
load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_0m_centers.mat')
data_0m = create_gene_buffer(data(gene_id, :),genlen);
data_0m = data_0m(NFR_pos);
load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_6h_centers.mat')
data_6h = create_gene_buffer(data(gene_id, :),genlen);
data_6h = data_6h(NFR_pos);

data_0m = conv(data_0m, gausswin(10)./sum(gausswin(10)), 'same');
data_6h = conv(data_6h, gausswin(10)./sum(gausswin(10)), 'same');
data_6h = data_6h .* (sum(data_0m) / sum(data_6h));
 
v = VideoWriter(path, 'MPEG-4');
v.FrameRate = args.frame_rate;
open(v);
[T, ~] = size(s_hist);
T = T/2;
s_hist = 1.*(s_hist>0);
F = figure('units','normalized','outerposition',[0 0 1 1]);
colormap(args.cmap);

%TODO: use interp2 to interpolate to simulation time
ker = fspecial('gaussian',[2,150],50);

% the loop that makes the movie:
for i = 500:round(args.frame_overlap*args.sample_frame):T
    clf;
    if i+args.sample_frame > T, break; end;
    
    % getting the current frame data:
    fr = s_hist(i:T+i , NFR_pos);
    centers = sum(fr);
    centers = conv(centers, gausswin(10)./sum(gausswin(10)), 'same');
    centers = centers ./ sum(centers);
    
    area(data_0m,'FaceColor',[0 0.9 0.9]);
    hold all;
    area(data_6h,'FaceColor',[0 0.6 0.6]);
    plot(centers .* sum((data_0m)), 'LineWidth',2);
    axis([1 length(NFR_pos) 0 max([max(data_0m) max(data_6h) max(centers .* sum((data_0m+data_6h)/2))])]);
    legend('0m','6h','sim')
    title(['time = ' num2str(t(i)) char(10) 'simulation step = ' num2str(((i-500)/100))+1])
    
    %{
    % making the top screen:
    subplot(10,1,1:6)
    imagesc(fr)
    title(['Gene ' num2str(gene_id) ' Simulation'])
    set(gca,'xtick',[],'ytick',1:50:args.sample_frame,...
        'yticklabel',i+1:50:i+args.sample_frame, 'ticklength', [0 0]);
    ylabel('time');
    
    % making the middle screen:
    subplot(10,1,7:8)
    area(nuc_sum./max(nuc_sum), 'FaceColor', [0.8 0.8 0.8])
    hold all;
    plot(coverage./max(coverage), 'b', 'linewidth', 2);
    plot(centers./max(centers), 'r', 'linewidth', 2)
    set(gca,'ytick', [], 'ticklength', [0 0]);
    legend('total sim centers', 'Location', 'NorthWest')
    set(gcf, 'Color',[1 1 1], 'PaperPosition',[0 0 args.frame_size], ...
        'PaperSize', args.frame_size);
    
    % making the lower screen:
    subplot(10,1,9:10)
    rectangle('Position',[1350 0 600 1],'FaceColor','y')
    hold all
    plot(polyA./max(polyA),'b', 'linewidth', 2)
    plot(polyT./max(polyT),'r', 'linewidth', 2)
    legend('PolyA (pushes right)','PolyT (pushes left)', 'Location', 'NorthWest')
    xlabel(['Position (TSS at ' num2str(fix(length(s_hist(1,:))/2)) ')'])
    %}
    
    fr = getframe(F);
    writeVideo(v, fr);
    pause(.1);
end

close(v)