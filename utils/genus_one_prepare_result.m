function [source_vertex_reg_intersect_index, correspondence_mask, target_correspondence_mask...
    source_vertex_reg_3D, source_face_reg, displace, dist,...
    target_intensity_reg, intensity_diff] = ...
    genus_one_prepare_result(source_face, source_vertex, source_vertex_reg,...
                    target_face, target_vertex, flat_target_vertex,...
                    target_outer_boundary_index, target_inner_boundary_index, target_corner_index,...
                    source_intensity, target_intensity,...
                    source_inner_boundary_index)
               
%% find the common domain submesh from inconsistent shape registration

% assume only one boundary component

target_boundary_pos = flat_target_vertex(target_inner_boundary_index,:);

% find the vertices in the common domain
source_vertex_reg_intersect_index = ~in_or_out(target_boundary_pos, source_vertex_reg);   
source_vertex_reg_intersect_pos = source_vertex_reg(source_vertex_reg_intersect_index,:);

correspondence_mask = zeros(size(source_vertex_reg,1), 1);
correspondence_mask(source_vertex_reg_intersect_index, 1) = 1;

source_reg_boundary_pos = source_vertex_reg(source_inner_boundary_index,:);
target_vertex_intersect_index = ~in_or_out(source_reg_boundary_pos, flat_target_vertex);   
target_correspondence_mask = zeros(size(target_vertex,1), 1);
target_correspondence_mask(target_vertex_intersect_index, 1) = 1;

tmp = 1:size(source_vertex, 1);
source_vertex_reg_intersect_index = tmp(source_vertex_reg_intersect_index);

% extend target vertex for interpolation
[target_outer_boundary_index, target_corner_offset] = genus_one_outer_boundary_segment(target_outer_boundary_index, target_corner_index);
% [target_outer_boundary_index, target_corner_offset] = genus_one_outer_boundary_segment(target_outer_boundary_index, target_corner_index);
target_corner_id = target_outer_boundary_index(target_corner_offset); % ordered
% target_corner_id = target_outer_boundary_index(target_corner_offset); % ordered
basis1 = flat_target_vertex(target_corner_id(3),:) - flat_target_vertex(target_corner_id(2),:);
basis2 = flat_target_vertex(target_corner_id(1),:) - flat_target_vertex(target_corner_id(2),:);
nv_target = size(flat_target_vertex,1);
nf_target = size(target_face,1);
flat_target_vertex_ext = zeros(nv_target*9,2);
target_face_ext = zeros(nf_target*9,3);
target_vertex_ext_coord = repmat(target_vertex, 9, 1);
target_intensity_ext = repmat(target_intensity, 9, 1);
% figure; hold on;
for j = 0:2
    for i = 0:2
        flat_target_vertex_ext((i + 3*j)*nv_target+(1:nv_target),:) = flat_target_vertex + (-1 + j)*repmat(basis1, nv_target, 1) + (-1 + i)*repmat(basis2, nv_target, 1);
        target_face_ext((i + 3*j)*nf_target+(1:nf_target),:) = target_face + (i + 3*j)*nv_target;
%         plot(source_vertex_ext(:,1), source_vertex_ext(:,2), '.'); drawnow; pause;
    end
end

% get barycentric coordinate of registration with respect to the target mesh 
Tri_target = triangulation(target_face_ext, flat_target_vertex_ext);
tid = pointLocation(Tri_target, source_vertex_reg_intersect_pos);
source_vertex_reg_bary = cartesianToBarycentric(Tri_target, tid, source_vertex_reg_intersect_pos);
% use the barycentric coordinate to get 3D registration
Tri_target_3D = triangulation(target_face_ext, target_vertex_ext_coord);
source_vertex_reg_3D_sub = barycentricToCartesian(Tri_target_3D, tid, source_vertex_reg_bary);
source_vertex_reg_3D = arap(source_vertex, source_face, source_vertex_reg_intersect_index', source_vertex_reg_3D_sub);
% find the faces in the common domain
f_list = ismember(source_face(:,1),source_vertex_reg_intersect_index') | ismember(source_face(:,2),source_vertex_reg_intersect_index') |...
    ismember(source_face(:,3), source_vertex_reg_intersect_index');
source_face_reg = source_face(f_list,:);
% interpolate target intensity to registered
target_intensity_reg = griddata(flat_target_vertex_ext(:,1), flat_target_vertex_ext(:,2), target_intensity_ext, source_vertex_reg(:,1), source_vertex_reg(:,2));
target_intensity_reg(isnan(target_intensity_reg)) = 0;
% distance from initial alignment
displace = source_vertex_reg_3D - source_vertex;
% displace_3D = displace_3D - mean(displace_3D, 2);
dist = sqrt(sum((displace).^2, 2));
intensity_diff = abs(source_intensity - target_intensity_reg);
end