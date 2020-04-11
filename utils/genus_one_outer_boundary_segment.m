function [source_outer_boundary_index, source_corner_offset] = genus_one_outer_boundary_segment(source_outer_boundary_index, source_corner_index)

source_corner_offset = ismember(source_outer_boundary_index,source_corner_index);
tmp = 1:size(source_outer_boundary_index,1);
source_corner_offset = tmp(source_corner_offset);
source_outer_boundary_index = circshift(source_outer_boundary_index, -source_corner_offset(1)+1, 1);
source_corner_offset = source_corner_offset - source_corner_offset(1) + 1;
end