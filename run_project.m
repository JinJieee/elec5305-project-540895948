%% ELEC5305 project - auditory cues for walking and turning
clear; close all; clc;
folder = fileparts(mfilename('fullpath'));
addpath(fullfile(folder, 'helpers'));
resultFolder = fullfile(folder, 'results');
audioFolder = fullfile(folder, 'audio');
if ~isfolder(resultFolder), mkdir(resultFolder); end
if ~isfolder(audioFolder), mkdir(audioFolder); end

%% Trial settings
fs = 48000;
duration = 20;
N = fs * duration;
frameLength = 1024;
hop = 256;
route.speed = 1.5;
route.walkStart = 2;
route.turnStart = 8;
route.turnDuration = 4;
route.walkStop = 18;
starts = (-frameLength + hop : hop : N - 1)';
frameTime = min(duration, max(0, (starts + frameLength/2) / fs));

%% Load the course HRIRs and keep the horizontal directions
data = load(fullfile(folder, 'data', 'HRIR_ANL+_cj.mat'));
hrir = data.hrir;
assert(hrir.sampFreq == fs);
horizontal = abs(hrir.sourceSphCoord(:,2)) < 1e-8;
[azimuth, order] = sort(hrir.sourceSphCoord(horizontal,1) * 180/pi);
responses = hrir.impulseResponses(:, :, horizontal);
responses = responses(:, :, order);

% One gain for the whole database preserves differences between directions.
energy = reshape(sum(sum(responses.^2, 1), 2), [], 1);
responses = responses / mean(sqrt(energy / 2));
azimuth = [azimuth; azimuth(1) + 360];
responses = cat(3, responses, responses(:,:,1));
filterLength = size(responses, 1);

%% Make the same source signals and footstep stream for every condition
[sources, footstep, footstepLevels, footstepInfo] = make_signals(fullfile(folder, 'data'), fs, duration, route);
feet = [footstep; zeros(filterLength-1, 2)];
writetable(footstepLevels, fullfile(resultFolder, 'footstep_levels.csv'));

names = {'straight', 'left_world', 'right_world', 'left_head', 'right_head', 'left_inside'};
turns = [0, 1, -1, 1, -1, 1];
mixtures = cell(6,1);
landmarks = cell(6,1);
references = cell(6,1);
geometry = cell(6,1);
peak = 0;

%% Render the six conditions
for k = 1:numel(names)
    fprintf('Rendering %s...\n', names{k});
    [position, heading] = listener_path(frameTime, turns(k), route);
    reference = zeros(size(feet));
    referenceY = [50, -50];
    for source = 1:2
        angle = atan2(referenceY(source)-position(:,2), 100-position(:,1)) - heading;
        reference = reference + render_binaural(sources(:,source), angle, responses, azimuth, hop);
    end

    angle = zeros(size(heading));
    landmark = zeros(size(feet));
    if k > 1
        sourceY = -60 * turns(k);
        if k == 6, sourceY = 60; end
        angle = unwrap(atan2(sourceY-position(:,2), 100-position(:,1)) - heading);
        if k == 4 || k == 5
            % Match the outside landmark at turn onset, then hold its direction.
            entryLength = route.speed * (route.turnStart-route.walkStart);
            angle(:) = atan2(sourceY, 100-entryLength);
        end
        landmark = render_binaural(sources(:,3), angle, responses, azimuth, hop);
    end
    mixtures{k} = feet + reference + landmark;
    landmarks{k} = landmark;
    references{k} = reference;
    geometry{k} = [frameTime, position, heading, angle];
    peak = max([peak; abs(mixtures{k}(:)); abs(landmark(:)); abs(reference(:)); abs(feet(:))]);
end

%% Export audio and measure the isolated landmark
outputGain = 0.89 / peak;
cue = cell(6,1);
summary = zeros(6,5);
[~, turnStart] = min(abs(frameTime-route.turnStart));
[~, turnEnd] = min(abs(frameTime-route.turnStart-route.turnDuration));

for k = 1:numel(names)
    audiowrite(fullfile(audioFolder, [names{k} '.wav']), mixtures{k}*outputGain, fs, 'BitsPerSample', 24);
    if k == 1
        cue{k} = measure_cues(references{k}, fs, N);
    else
        cue{k} = measure_cues(landmarks{k}, fs, N);
    end
    writetable(cue{k}, fullfile(resultFolder, [names{k} '_cues.csv']));
    pathTable = array2table(geometry{k}, 'VariableNames', ...
        {'Time_s', 'X_m', 'Y_m', 'Heading_rad', 'Azimuth_rad'});
    writetable(pathTable, fullfile(resultFolder, [names{k} '_path.csv']));

    change = geometry{k}(turnEnd,5) - geometry{k}(turnStart,5);
    gain = NaN;
    if k > 1
        headingChange = geometry{k}(turnEnd,4) - geometry{k}(turnStart,4);
        gain = -change / headingChange;
    end
    selected = cue{k}.Time_s >= 8.1 & cue{k}.Time_s <= 11.9;
    ildRange = max(cue{k}.ILD_dB(selected)) - min(cue{k}.ILD_dB(selected));
    itdRange = max(cue{k}.ITD_us(selected)) - min(cue{k}.ITD_us(selected));
    if k == 1, itdRange = NaN; end  % Two sources do not have a single landmark ITD.
    summary(k,:) = [gain, change*180/pi, ildRange, itdRange, max(abs(mixtures{k}(:)*outputGain))];
end
summary = array2table(summary, 'VariableNames', ...
    {'TurnGain', 'AzimuthChange_deg', 'ILDRange_dB', 'ITDRange_us', 'Peak'});
summary = addvars(summary, names', 'Before', 1, 'NewVariableNames', 'Condition');
writetable(summary, fullfile(resultFolder, 'summary.csv'));

% These two clips omit the reference pair, so only the landmark frame changes.
audiowrite(fullfile(audioFolder, 'world_comparison.wav'), (feet+landmarks{2})*outputGain, fs, 'BitsPerSample', 24);
audiowrite(fullfile(audioFolder, 'head_comparison.wav'), (feet+landmarks{4})*outputGain, fs, 'BitsPerSample', 24);
audiowrite(fullfile(audioFolder, 'footsteps_only.wav'), feet*outputGain, fs, 'BitsPerSample', 24);

%% Check the fixed filter case and compare update rates
fixedAngle = 0.3 * ones(size(frameTime));
fixed = render_binaural(sources(:,3), fixedAngle, responses, azimuth, hop);
h = interpolate_hrir(0.3, azimuth, responses);
ordinary = [conv(sources(:,3),h(:,1)), conv(sources(:,3),h(:,2))];
convolutionError = max(abs(fixed(:)-ordinary(:)));
assert(convolutionError < 1e-11);
assert(footstepInfo.gainError_dB < 1e-10);
assert(max(summary.Peak) <= 0.890001);

hops = [512, 256, 128];
rendered = cell(3,1);
seconds = zeros(3,1);
for k = 1:3
    s = (-frameLength+hops(k):hops(k):N-1)';
    t = min(duration, max(0, (s+frameLength/2)/fs));
    [p, heading] = listener_path(t, 1, route);
    angle = atan2(-60-p(:,2), 100-p(:,1)) - heading;
    timer = tic;
    rendered{k} = render_binaural(sources(:,3), angle, responses, azimuth, hops(k));
    seconds(k) = toc(timer);
end
errorDB = zeros(3,1);
for k = 1:2
    errorDB(k) = 20*log10(norm(rendered{k}-rendered{3},'fro') / norm(rendered{3},'fro'));
end
errorDB(3) = NaN; % Reference output, rather than an independent error measurement.
hopStudy = table(hops', seconds, errorDB, 'VariableNames', {'Hop_samples','Render_s','Error_dB'});
writetable(hopStudy, fullfile(resultFolder, 'hop_study.csv'));

%% Compare landmark positions at matched turn-onset distances
placementStudy = placement_study(folder, sources(:,3), responses, azimuth, ...
    frameTime, route, fs, hop, N);

%% Save results and draw the figures
save(fullfile(resultFolder, 'results.mat'), 'summary', 'cue', 'geometry', 'footstepLevels', 'footstepInfo', ...
    'route', 'fs', 'duration', 'hop', 'hopStudy', 'placementStudy', 'convolutionError', 'outputGain');
plot_results(folder, names, geometry, cue, references, footstepLevels, hopStudy, fs, N, landmarks);
fid = fopen(fullfile(resultFolder, 'run_log.txt'), 'w');
fprintf(fid, 'MATLAB %s\nStatic convolution error: %.12g\n', version, convolutionError);
fprintf(fid, 'Source: fountain_mono.wav; landmark 0-20 s, references 20-40 and 40-60 s\n');
fprintf(fid, 'Source RMS: references 0.012, landmark 0.018; roll-off 3.5-6.5 kHz\n');
fprintf(fid, 'Binaural source interval: %.2f to %.2f s\n', footstepInfo.sourceStart_s, footstepInfo.sourceDuration_s);
fprintf(fid, 'Common footstep gain: %.12g dB\n', footstepInfo.gain_dB);
fprintf(fid, 'Footstep gain error: %.12g dB\n', footstepInfo.gainError_dB);
fprintf(fid, 'Original segment L/R level difference: %.12g dB\n', footstepInfo.inputILD_dB);
fprintf(fid, 'Common output gain: %.12g\n', outputGain);
fclose(fid);
disp(summary);
disp(hopStudy);
fprintf('Finished. Audio and results are in the project folder.\n');
