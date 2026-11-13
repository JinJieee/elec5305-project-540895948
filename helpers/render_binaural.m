function output = render_binaural(source, angle, responses, azimuth, hop)
% Window the input, filter each frame, then keep the complete filter tails.
source = source(:);
N = numel(source);
frameLength = 1024;
filterLength = size(responses, 1);
starts = (-frameLength + hop : hop : N - 1)';
window = 0.5 - 0.5 * cos(2*pi*(0:frameLength-1)'/frameLength);
weight = zeros(N, 1);

for frame = 1:numel(starts)
    first = max(0, starts(frame));
    last = min(N, starts(frame) + frameLength);
    samples = (first:last-1)';
    weight(samples + 1) = weight(samples + 1) + window(samples - starts(frame) + 1);
end

output = zeros(N + filterLength - 1, 2);
fftLength = 2^nextpow2(frameLength + filterLength - 1);
for frame = 1:numel(starts)
    start = starts(frame);
    first = max(0, start);
    last = min(N, start + frameLength);
    samples = (first:last-1)';
    block = zeros(frameLength, 1);
    block(samples - start + 1) = source(samples + 1) .* ...
        window(samples - start + 1) ./ weight(samples + 1);

    h = interpolate_hrir(angle(frame), azimuth, responses);
    filtered = real(ifft(fft(block, fftLength) .* fft(h, fftLength), [], 1));
    filtered = filtered(1:frameLength + filterLength - 1, :);
    first = max(0, start);
    last = min(size(output, 1), start + size(filtered, 1));
    output(first + 1:last, :) = output(first + 1:last, :) + ...
        filtered(first - start + 1:last - start, :);
end
end
