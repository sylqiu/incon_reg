function normalized_intensity = intensity_normalization(intensity)

v_min = min(intensity(:)); v_max = max(intensity(:));

normalized_intensity = (intensity - v_min) ./ (v_max - v_min);


end