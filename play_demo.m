function play_demo(recordVideo)
% Play the two comparison clips with their listener path and source direction.
% play_demo(true) saves silent AVI previews for video editing instead.
if nargin == 0, recordVideo = false; end
folder = fileparts(mfilename('fullpath'));
data = load(fullfile(folder, 'results', 'results.mat'));
files = {'world_comparison.wav', 'head_comparison.wav'};
conditions = [2, 4];
labels = {'World-fixed landmark', 'Head-fixed landmark'};

for clip = 1:2
    [audio, fs] = audioread(fullfile(folder, 'audio', files{clip}));
    path = data.geometry{conditions(clip)};
    duration = size(audio,1)/fs;
    f = figure('Color','w','Position',[100 100 1000 550]);
    subplot(1,2,1); plot(path(:,2),path(:,3),'Color',[0.7 0.7 0.7]); hold on;
    listener = plot(path(1,2),path(1,3),'bo','MarkerFaceColor','b','MarkerSize',8);
    facing = quiver(path(1,2),path(1,3),2,0,0,'k','LineWidth',1.5);
    axis equal; xlim([-2 25]); ylim([-3 16]); grid on;
    xlabel('World X (m)'); ylabel('World Y (m)'); title('Listener position and heading');
    subplot(1,2,2); plot(path(:,1),path(:,5)*180/pi,'LineWidth',1.3); hold on;
    marker = plot(0,path(1,5)*180/pi,'ro','MarkerFaceColor','r');
    xline(8,'--'); xline(12,'--'); grid on; xlim([0 20]); ylim([-160 10]);
    xlabel('Time (s)'); ylabel('Head-relative landmark angle (deg)');
    heading = sgtitle(labels{clip});

    if recordVideo
        filename = fullfile(folder,'results',sprintf('preview_%d.avi',clip));
        movie = VideoWriter(filename,'Motion JPEG AVI');
        movie.FrameRate = 10;
        open(movie);
        times = 0:0.1:20;
    else
        player = audioplayer(audio,fs);
        play(player);
        timer = tic;
        times = 0:0.04:duration;
    end

    for frame = 1:numel(times)
        if ~isgraphics(f), break; end
        if recordVideo, t = times(frame); else, t = toc(timer); end
        [~, index] = min(abs(path(:,1)-t));
        x = path(index,2); y = path(index,3); psi = path(index,4);
        set(listener,'XData',x,'YData',y);
        set(facing,'XData',x,'YData',y,'UData',2*cos(psi),'VData',2*sin(psi));
        set(marker,'XData',t,'YData',path(index,5)*180/pi);
        heading.String = sprintf('%s - %.1f s',labels{clip},t);
        drawnow;
        if recordVideo
            imageFile = fullfile(tempdir,'elec5305_demo_frame.png');
            print(f,imageFile,'-dpng','-r80');
            writeVideo(movie,imread(imageFile));
        else
            if t >= duration, break; end
            pause(0.04);
        end
    end
    if recordVideo, close(movie); else, stop(player); end
    if isgraphics(f), close(f); end
end
end
