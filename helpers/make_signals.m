function [sources, footsteps, levels, info] = make_signals(folder, fs, duration, route)
N = round(duration * fs);
t = (0:N-1)' / fs;
fade = min(1, min(t/0.08, (duration-t)/0.08));
[fountain, sourceRate] = audioread(fullfile(folder, 'fountain_mono.wav'));
assert(sourceRate == fs && size(fountain,2) == 1);
assert(numel(fountain) >= 3*N);

% Use different portions for the references, but the same landmark in every trial.
sourceStarts = [20, 40, 0];
sourceRMS = [0.012, 0.012, 0.018];
frequency = (0:N-1)' * fs/N;
frequency = min(frequency, fs-frequency);
lowpass = ones(N,1);
lowpass(frequency >= 6500) = 0;
rolloff = frequency > 3500 & frequency < 6500;
lowpass(rolloff) = 0.5 + 0.5*cos(pi*(frequency(rolloff)-3500)/3000);
sources = zeros(N,3);
for k = 1:3
    first = round(sourceStarts(k)*fs) + 1;
    recordingPart = fountain(first:first+N-1);
    recordingPart = recordingPart - mean(recordingPart);
    recordingPart = real(ifft(fft(recordingPart).*lowpass));
    recordingPart = recordingPart * (sourceRMS(k)/sqrt(mean(recordingPart.^2)));
    sources(:,k) = recordingPart .* fade;
end

% Read the binaural recording without mixing, aligning or normalizing its ears.
[recording, sampleRate] = audioread(fullfile(folder, 'balanced_denoised_clip.wav'));
assert(sampleRate == fs && size(recording,2) == 2);
walkLength = round((route.walkStop-route.walkStart) * fs);
assert(size(recording,1) >= walkLength);
segment = recording(1:walkLength, :);
info.sourceStart_s = 0;
info.sourceDuration_s = walkLength / fs;
info.sourceRecordingDuration_s = size(recording,1) / fs;
info.inputRMS = sqrt(mean(segment.^2, 1));
info.gain = 0.045 / sqrt(mean(segment.^2, 'all'));

% Only a common gain and common 40-ms boundary fades are applied.
edgeLength = round(0.04 * fs);
envelope = ones(walkLength,1);
envelope(1:edgeLength) = linspace(0,1,edgeLength)';
envelope(end-edgeLength+1:end) = linspace(1,0,edgeLength)';
playback = segment .* envelope * info.gain;
footsteps = zeros(N,2);
first = round(route.walkStart * fs) + 1;
footsteps(first:first+walkLength-1, :) = playback;

info.inputILD_dB = 20*log10(info.inputRMS(1)/info.inputRMS(2));
info.gain_dB = 20*log10(info.gain);

% Use interior windows so that the boundary fades do not enter the gain check.
windowLength = round(0.25 * fs);
starts = windowLength:windowLength:walkLength-2*windowLength;
values = zeros(numel(starts),9);
for k = 1:numel(starts)
    samples = starts(k)+(1:windowLength);
    inputRMS = sqrt(mean(segment(samples,:).^2,1));
    outputRMS = sqrt(mean(playback(samples,:).^2,1));
    gain = 20*log10(outputRMS ./ inputRMS);
    sourceTime = (starts(k)+windowLength/2)/fs;
    values(k,:) = [sourceTime, route.walkStart+sourceTime, route.speed*sourceTime, ...
        inputRMS, outputRMS, gain];
end
levels = array2table(values, 'VariableNames', ...
    {'SourceTime_s','TrialTime_s','Distance_m','InputRMS_L','InputRMS_R', ...
    'PlaybackRMS_L','PlaybackRMS_R','Gain_L_dB','Gain_R_dB'});
info.gainError_dB = max(abs([levels.Gain_L_dB;levels.Gain_R_dB]-info.gain_dB));
end
