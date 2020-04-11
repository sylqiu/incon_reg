function [g3,g4,Mg,Sg,mask,scale] = genus_one_combine_to_same_grid(source_vertex, source_intensity, ...
                                target_vertex, target_intensity,...
                                source_inner_boundary_index, target_inner_boundary_index,...
                            source_outer_boundary_index,...
                            source_corner_index)

%output Mg and Sg defined one the grid [g3,g4]
%ready to be used for itensity matching

% the rectangular grid (g3, g4) should contain V_t
sx = max(target_vertex(:,1)) - min(target_vertex(:,1));
sy = max(target_vertex(:,2)) - min(target_vertex(:,2));
delta = 0.4;
sx = (1 + delta)*sx;
sy = (1 + delta)*sy;
r = sx/sy; N = round(sqrt(size(target_vertex,1)))*(1 + delta);
[g4,g3] = ndgrid(1*N:-0.4:0, 0:0.4:1*round(r*N)); %natural ordering

%scale and translate the grid
g4 = g4/N*sy; g3 = g3/r/N*sx;
g4 = g4 + min(target_vertex(:,2)) - delta/2*sy; 
g3 = g3 + min(target_vertex(:,1)) - delta/2*sx;

scale = abs((g4(end) - g4(1)) / size(g4,2));
% Ftmp = delaunay(g3(:),g4(:));
% L = cotmatrix([g3(:),g4(:)],Ftmp);

% replicate source and target vertices
nv_source = size(source_vertex,1);
nv_target = size(target_vertex,1);
source_vertex_ext = zeros(nv_source*9,2);
target_vertex_ext = zeros(nv_target*9,2);
source_intensity_ext = repmat(source_intensity, 9, 1);
target_intensity_ext = repmat(target_intensity, 9, 1);
%
[source_outer_boundary_index, source_corner_offset] = genus_one_outer_boundary_segment(source_outer_boundary_index, source_corner_index);
% [target_outer_boundary_index, target_corner_offset] = genus_one_outer_boundary_segment(target_outer_boundary_index, target_corner_index);
source_corner_id = source_outer_boundary_index(source_corner_offset); % ordered
% target_corner_id = target_outer_boundary_index(target_corner_offset); % ordered
basis1 = source_vertex(source_corner_id(3),:) - source_vertex(source_corner_id(2),:);
basis2 = source_vertex(source_corner_id(1),:) - source_vertex(source_corner_id(2),:);
% figure; hold on;
for j = 0:2
    for i = 0:2
        source_vertex_ext((i + 3*j)*nv_source+(1:nv_source),:) = source_vertex + (-1 + j)*repmat(basis1, nv_source, 1) + (-1 + i)*repmat(basis2, nv_source, 1);
        target_vertex_ext((i + 3*j)*nv_target+(1:nv_target),:) = target_vertex + (-1 + j)*repmat(basis1, nv_target, 1) + (-1 + i)*repmat(basis2, nv_target, 1);
%         plot(source_vertex_ext(:,1), source_vertex_ext(:,2), '.'); drawnow; pause;
    end
end

% interpolating to [g3,g4]
intpM = scatteredInterpolant(source_vertex_ext(:,1),source_vertex_ext(:,2), source_intensity_ext);
intpS= scatteredInterpolant(target_vertex_ext(:,1),target_vertex_ext(:,2), target_intensity_ext);
Mg = intpM(g3,g4);
Sg = intpS(g3,g4);
% 
bdy_xym = source_vertex(source_inner_boundary_index,:);
bdy_xyt = target_vertex(target_inner_boundary_index,:);
bin_m = in_or_out(bdy_xym,[g3(:),g4(:)]); bin_m = reshape(bin_m,size(g3,1),size(g3,2));
bin_t = in_or_out(bdy_xyt,[g3(:),g4(:)]); bin_t = reshape(bin_t,size(g3,1),size(g3,2));
mask = isnan(Mg) | isnan(Sg) | bin_m | bin_t;
Mg(mask) = 0; Sg(mask) = 0;
mask = ~mask(:)*1.0;
% for kkk = 1:10
%     mask = mask + 0.1*L*mask;
% end
mask = reshape(mask,size(g3,1),size(g3,2)); 
% figure(10); imshow(mask); title('Overlap domain');
Mg = Mg.*mask; Sg = Sg.*mask;
                            
end