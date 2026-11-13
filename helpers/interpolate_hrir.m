function h = interpolate_hrir(angle, azimuth, responses)
% Blend neighbouring measured responses. The last direction repeats the first.
angle = mod(angle * 180/pi + 180, 360) - 180;
left = find(azimuth <= angle, 1, 'last');
left = min(left, numel(azimuth) - 1);
fraction = (angle - azimuth(left)) / (azimuth(left + 1) - azimuth(left));
h = (1 - fraction) * responses(:, :, left) + fraction * responses(:, :, left + 1);
end
