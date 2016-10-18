function nucsim_movie(s_hist, t, path, seq, varargin)
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

defs = struct('sample_frame', 200, 'frame_overlap', .2, ...
              'frame_rate', 12, 'cmap', AdvancedColormap('w wy wyr wrk rrk'), ...
              'frame_size', [16,8]);
args = parse_namevalue_pairs(defs, varargin);

[ PolyA_Sites, PolyT_Sites, REB1_Sites, ABF1_Sites, RAP1_Sites ] = Extract_Sites_From_Gene(seq, 2500);
PolyA_Sites = [PolyA_Sites, zeros(1,999)];
PolyT_Sites = [PolyT_Sites, zeros(1,999)];
polyA = conv(PolyA_Sites, gausswin(40), 'same');
polyT = conv(PolyT_Sites, gausswin(40), 'same');

v = VideoWriter(path, 'MPEG-4');
v.FrameRate = args.frame_rate;
open(v);

[T, ~] = size(s_hist);
s_hist = 1.*(s_hist>0);
F = figure;
colormap(args.cmap);

%TODO: use interp2 to interpolate to simulation time
ker = fspecial('gaussian',[2,150],50);
for i = 0:round(args.frame_overlap*args.sample_frame):T
    clf;
    if i+args.sample_frame > T, break; end;
    fr = s_hist(i+1:i+args.sample_frame,:); % get the current s_hist data
    centers = sum(fr);
    fr = conv2(fr, ker, 'same');
    coverage = sum(fr);
    subplot(10,1,1:6)
    imagesc(fr)
    set(gca,'xtick',[],'ytick',1:50:args.sample_frame,...
        'yticklabel',i+1:50:i+args.sample_frame, 'ticklength', [0 0]);
    ylabel('time');
    
    subplot(10,1,7:8)
    rectangle('Position',[600 0 400 1],'FaceColor','y')
    hold all
    plot(polyA./max(polyA),'b', 'linewidth', 2)
    plot(polyT./max(polyT),'r', 'linewidth', 2)
    legend('PolyA','PolyT')
    
    subplot(10,1,9:10)
    plot(coverage./max(coverage), 'linewidth', 2);
    hold all;
    plot(centers./max(centers), 'r', 'linewidth', 2)
    set(gca,'ytick', [], 'ticklength', [0 0]);
    xlabel('position')
    set(gcf, 'Color',[1 1 1], 'PaperPosition',[0 0 args.frame_size], ...
        'PaperSize', args.frame_size);
    fr = getframe(F);
    writeVideo(v, fr);
    pause(.1);
end

close(v)