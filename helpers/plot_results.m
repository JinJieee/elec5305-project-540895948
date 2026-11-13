function plot_results(folder, names, geometry, cue, references, footstepLevels, hopStudy, fs, N, landmarks)
resultFolder = fullfile(folder, 'results');
labels = {'Straight', 'Left world', 'Right world', 'Left head', 'Right head', 'Left inside'};
set(groot, 'defaultAxesFontSize', 11, 'defaultFigureColor', 'w');

f = figure('Position', [100 100 1000 420]);
subplot(1,2,1); hold on;
for k = 1:3
    plot(geometry{k}(:,2), geometry{k}(:,3), 'LineWidth', 1.4, 'DisplayName', labels{k});
end
axis equal; grid on; xlabel('World X (m)'); ylabel('World Y (m)'); legend('Location','best');
title('Listener paths');
subplot(1,2,2); hold on;
plot([0 9 12.82], [0 0 3.82], 'k-', 'LineWidth', 1.3);
scatter([100 100], [50 -50], 45, [0.4 0.4 0.4], 'filled');
scatter(100, -60, 55, 'r', 'filled'); scatter(100, 60, 55, 'b', 'filled');
text(107,-60,'Outside');
text(107,60,'Inside');
xlim([-5 140]); ylim([-80 80]);
xlabel('World X (m)'); ylabel('World Y (m)'); grid on; title('Left-turn source positions');
save_plot(f, fullfile(resultFolder, 'paths.png'));

f = figure('Position', [100 100 950 600]);
subplot(2,1,1); hold on;
for k = 2:6
    plot(geometry{k}(:,1), geometry{k}(:,5)*180/pi, 'LineWidth', 1.2, 'DisplayName', labels{k});
end
xline(8,'--','HandleVisibility','off'); xline(12,'--','HandleVisibility','off'); grid on; ylabel('Head-relative azimuth (deg)');
legend('Location','best','NumColumns',3); title('Landmark direction');
subplot(2,1,2); hold on;
for k = 2:3
    plot(geometry{k}(:,1), geometry{k}(:,4)*180/pi, 'LineWidth', 1.2, 'DisplayName', labels{k});
end
xline(8,'--','HandleVisibility','off'); xline(12,'--','HandleVisibility','off'); grid on; xlabel('Time (s)'); ylabel('Heading (deg)');
save_plot(f, fullfile(resultFolder, 'azimuth.png'));

f = figure('Position', [100 100 950 600]);
for panel = 1:2
    subplot(2,1,panel); hold on;
    for k = [2 3 4 6]
        if panel == 1, values = cue{k}.ILD_dB; else, values = cue{k}.ITD_us; end
        plot(cue{k}.Time_s, values, 'LineWidth', 1, 'DisplayName', labels{k});
    end
    xline(8,'--','HandleVisibility','off'); xline(12,'--','HandleVisibility','off'); grid on;
    if panel == 1
        ylabel('ILD, left/right (dB)'); legend('Location','best','NumColumns',4);
        title('Isolated landmark: 50-ms analysis window');
    else
        ylabel('ITD, left minus right (us)'); xlabel('Time (s)');
    end
end
save_plot(f, fullfile(resultFolder, 'binaural_cues.png'));

f = figure('Position', [100 100 950 600]);
subplot(2,1,1); hold on;
for k = 1:3
    balance = measure_cues(references{k}, fs, N);
    average = conv(balance.ILD_dB, ones(25,1)/25, 'valid');
    plot(balance.Time_s(13:end-12), average, 'LineWidth', 1.2, 'DisplayName', labels{k});
end
xline(8,'--','HandleVisibility','off'); xline(12,'--','HandleVisibility','off'); grid on; legend('Location','best');
xlabel('Time (s)'); ylabel('Left/right balance (dB)'); title('Reference pair, 0.5-s average');
subplot(2,1,2); hold on;
plot(footstepLevels.Distance_m, 20*log10(footstepLevels.PlaybackRMS_L), 'LineWidth', 1.2, 'DisplayName','Left ear');
plot(footstepLevels.Distance_m, 20*log10(footstepLevels.PlaybackRMS_R), 'LineWidth', 1.2, 'DisplayName','Right ear');
grid on; xlabel('Distance travelled (m)'); ylabel('250-ms RMS (dBFS)'); legend('Location','best');
title('Continuous binaural recording, before the pooled output gain');
save_plot(f, fullfile(resultFolder, 'balance_footsteps.png'));

f = figure('Position', [100 100 900 380]);
subplot(1,2,1); bar(hopStudy.Hop_samples(1:2), hopStudy.Error_dB(1:2));
xlabel('Hop (samples)'); ylabel('Error relative to 128-sample hop (dB)'); grid on;
subplot(1,2,2); bar(hopStudy.Hop_samples, hopStudy.Render_s);
xlabel('Hop (samples)'); ylabel('Rendering time (s)'); grid on;
save_plot(f, fullfile(resultFolder, 'hop_study.png'));

f = figure('Position', [100 100 1050 620]);
windowLength = 1024; hop = 256;
window = 0.5 - 0.5*cos(2*pi*(0:windowLength-1)'/windowLength);
starts = 1:hop:N-windowLength+1;
time = (starts-1+windowLength/2)/fs;
frequency = (0:windowLength/2)'*fs/windowLength;
conditions = [2 4];
for row = 1:2
    for ear = 1:2
        spectrum = zeros(windowLength/2+1, numel(starts));
        for frame = 1:numel(starts)
            samples = starts(frame)+(0:windowLength-1);
            z = fft(landmarks{conditions(row)}(samples,ear).*window);
            spectrum(:,frame) = abs(z(1:windowLength/2+1))/sum(window);
        end
        subplot(2,2,(row-1)*2+ear);
        imagesc(time, frequency, 20*log10(spectrum+1e-12)); axis xy;
        ylim([0 10000]); clim([-90 -35]);
        scale = colorbar; scale.Label.String = 'Magnitude (dB)';
        xlabel('Time (s)'); ylabel('Frequency (Hz)');
        title(sprintf('%s, ear %d', labels{conditions(row)}, ear));
    end
end
save_plot(f, fullfile(resultFolder, 'spectrograms.png'));
end

function save_plot(f, filename)
print(f, filename, '-dpng', '-r160');
end
