function cue = measure_cues(audio, fs, signalLength)
% Positive ILD means more left-ear energy; positive ITD means left-ear delay.
windowLength = round(0.05 * fs);
hop = round(0.02 * fs);
starts = 0:hop:signalLength-windowLength;
lags = (-round(0.001 * fs):round(0.001 * fs))';
fftLength = 2^nextpow2(2 * windowLength - 1);
values = zeros(numel(starts), 4);

for frame = 1:numel(starts)
    samples = starts(frame) + (1:windowLength);
    block = audio(samples, :);
    block = block - mean(block, 1);
    energy = mean(block.^2, 1);
    ild = 10 * log10((energy(1) + 1e-20) / (energy(2) + 1e-20));

    correlation = real(ifft(fft(block(:,1), fftLength) .* ...
        conj(fft(block(:,2), fftLength))));
    correlation = correlation(mod(lags, fftLength) + 1);
    [peak, index] = max(correlation);
    offset = 0;
    if index > 1 && index < numel(correlation)
        denominator = correlation(index-1) - 2*correlation(index) + correlation(index+1);
        if abs(denominator) > 1e-20
            offset = 0.5 * (correlation(index-1) - correlation(index+1)) / denominator;
        end
    end
    itd = (lags(index) + offset) / fs * 1e6;
    coefficient = peak / (sqrt(sum(block(:,1).^2) * sum(block(:,2).^2)) + 1e-20);
    values(frame, :) = [(starts(frame) + windowLength/2)/fs, ild, itd, coefficient];
end
cue = array2table(values, 'VariableNames', ...
    {'Time_s', 'ILD_dB', 'ITD_us', 'Correlation'});
end
