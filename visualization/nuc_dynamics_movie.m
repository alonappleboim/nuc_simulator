function nuc_dynamics_movie(s_hist, t, path, gene_id, first_n_steps, window_size, timeExpansion, jump_len, varargin)
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
create_params_genome;
genlen = 3500;
TSS = 1750;
NFR_pos = [TSS-299:TSS+150];

defs = struct('sample_frame', 200, 'frame_overlap', .2, ...
              'frame_rate', 12, 'cmap', AdvancedColormap('w wy wyr wrk rrk'), ...
              'frame_size', [16,8]);
args = parse_namevalue_pairs(defs, varargin);

% get the data of the 0m and the 2h:
load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_0m_centers.mat')
data_0m = create_gene_buffer(data(gene_id, :),genlen);
data_0m = data_0m(NFR_pos);
load('C:\Users\Daniel\Documents\MATLAB\nuc_simulator\clustering\experiment_data\sth1_2h_centers.mat')
data_6h = create_gene_buffer(data(gene_id, :),genlen);
data_6h = data_6h(NFR_pos);

data_0m = conv(data_0m, gausswin(10)./sum(gausswin(10)), 'same');
data_6h = conv(data_6h, gausswin(10)./sum(gausswin(10)), 'same');
data_6h = data_6h .* (sum(data_0m) / sum(data_6h));
 
v = VideoWriter(path, 'MPEG-4');
v.FrameRate = args.frame_rate;
open(v);
T = 120 * timeExpansion;
%[T, ~] = size(s_hist);
%s_hist = 1.*(s_hist>0);
F = figure('units','normalized','outerposition',[0 0 1 1]);
colormap(args.cmap);

% the loop that makes the movie:
%for i = first_n_steps-(window_size):round(args.frame_overlap*args.sample_frame) / 2:first_n_steps+T*jump_len-(window_size)
for i = first_n_steps-(window_size) : 10 : first_n_steps+T*jump_len-(window_size)
    clf;
    if i+args.sample_frame > first_n_steps+T*jump_len, break; end;
    
    % getting the current frame data:
    fr = s_hist(i:(window_size)+i , NFR_pos);
    centers = sum(fr);
    centers = conv(centers, gausswin(10)./sum(gausswin(10)), 'same');
    centers = centers ./ sum(centers);
    
    area(data_0m,'FaceColor',[0.9 0.9 0.9]);
    hold all;
    area(data_6h,'FaceColor',[0.6 0.6 0.6]);
    plot(centers .* sum((data_0m)), 'LineWidth',2);
    axis([1 length(NFR_pos) 0 max([max(data_0m) max(data_6h) max(centers .* sum((data_0m+data_6h)/2))])]);
    legend('0 minutes','2 hours','Simulation')
    title(['Gene ' num2str(gene_id) char(10) 'Time Expansion = ' num2str(timeExpansion.*jump_len) char(10) 'Experiment Time In Minutes = ' num2str(round((i-(first_n_steps-(window_size)))/(jump_len*timeExpansion)))])
    xlabel('Position (TSS at 300)')
    ylabel('Nucleosome Intensity')
    
    fr = getframe(F);
    writeVideo(v, fr);
    pause(.1);
end

close(v)