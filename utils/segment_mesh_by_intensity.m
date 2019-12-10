function seg_mask = segment_mesh_by_intensity(intensity)
% assume intensity normalized
lower = quantile(intensity,0.05);
higher = quantile(intensity, 0.95);
seg_mask = ones(size(intensity, 1),1) * 0;
seg_mask(intensity <= lower) = 1;
seg_mask(intensity >= higher) = 1;

end