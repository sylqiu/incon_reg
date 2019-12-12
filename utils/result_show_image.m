function result_show_image(source_face, source_vertex, flat_source_vertex, source_vertex_reg,...
                    target_face, target_vertex, flat_target_vertex,...
                    source_intensity, target_intensity,...
                    landmark_source_index, landmark_target_index,...
                    source_vertex_reg_intersect_index, correspondence_mask, ...
                    source_vertex_reg_3D, source_face_reg, displace_2D, dist_2D,...
                    target_intensity_reg, intensity_diff)
% %% find the common domain submesh from inconsistent shape registration
% 
% % assume only one boundary component
% [~, target_boudary_index] = meshboundaries(target_face);
% target_boundary_pos = flat_target_vertex(target_boudary_index,:);
% 
% % find the vertices in the common domain
% source_vertex_reg_intersect_index = in_or_out(target_boundary_pos, source_vertex_reg);   
% source_vertex_reg_intersect_pos = source_vertex_reg(source_vertex_reg_intersect_index,:);
% correspondence_mask = zeros(size(source_vertex_reg,1), 1);
% correspondence_mask(source_vertex_reg_intersect_index, 1) = 1;
% tmp = 1:size(source_vertex, 1);
% source_vertex_reg_intersect_index = tmp(source_vertex_reg_intersect_index);
% % get barycentric coordinate of registration with respect to the target mesh 
% Tri_target = triangulation(target_face, flat_target_vertex);
% tid = pointLocation(Tri_target, source_vertex_reg_intersect_pos);
% source_vertex_reg_bary = cartesianToBarycentric(Tri_target, tid, source_vertex_reg_intersect_pos);
% % use the barycentric coordinate to get 3D registration
% Tri_target_3D = triangulation(target_face, target_vertex);
% source_vertex_reg_3D_sub = barycentricToCartesian(Tri_target_3D, tid, source_vertex_reg_bary);
% source_vertex_reg_3D = arap(source_vertex, source_face, source_vertex_reg_intersect_index', source_vertex_reg_3D_sub);
% % find the faces in the common domain
% f_list = ismember(source_face(:,1),source_vertex_reg_intersect_index') | ismember(source_face(:,2),source_vertex_reg_intersect_index') |...
%     ismember(source_face(:,3), source_vertex_reg_intersect_index');
% source_face_reg = source_face(f_list,:);
% % interpolate target intensity to registered
% target_intensity_reg = griddata(flat_target_vertex(:,1), flat_target_vertex(:,2), target_intensity, source_vertex_reg(:,1), source_vertex_reg(:,2));
% target_intensity_reg(isnan(target_intensity_reg)) = 0;
% % distance from initial alignment
% displace_2D = source_vertex_reg_3D - source_vertex;
% % displace_3D = displace_3D - mean(displace_3D, 2);
% dist_2D = sqrt(sum((displace_2D).^2, 2));
%% plotting figures
figure(1);
subplot(1,2,1); gpp_plot_mesh(source_face, source_vertex, source_intensity); colormap parula; title('Source');
hold on; plot(source_vertex(landmark_source_index,1), source_vertex(landmark_source_index,2), 'ro'); hold off;
subplot(1,2,2); gpp_plot_mesh(target_face, target_vertex, target_intensity); colormap parula; title('Target');
hold on; plot(target_vertex(landmark_target_index,1), target_vertex(landmark_target_index,2), 'ro'); hold off;


figure(2);
gpp_plot_mesh(source_face, source_vertex_reg, 'FaceColor', [36 169 225]/255, 'EdgeColor', [90 169 225]/255, 'FaceAlpha', 0.2); hold on;
gpp_plot_mesh(target_face, flat_target_vertex, 'FaceColor', [216,108,112]/255, 'EdgeColor', [216,108,112]/255, 'FaceAlpha', 0.2);
plot(source_vertex_reg(landmark_source_index,1), source_vertex_reg(landmark_source_index,2), 'bo');
plot(flat_target_vertex(landmark_target_index,1), flat_target_vertex(landmark_target_index,2), 'ro');
title('Overlaid view');
hold off;

figure(3);
% diff_ = abs(source_intensity - target_intensity_reg);
diff_ = dist_2D;
subplot(2,2,1); gpp_plot_mesh(source_face_reg, source_vertex, diff_); colormap parula; title('Registered w/ displacement mag');
subplot(2,2,4); gpp_plot_mesh(target_face, target_vertex, target_intensity); colormap parula; title('Target');
subplot(2,2,3); gpp_plot_mesh(source_face_reg, source_vertex_reg, source_intensity); colormap parula; title('Source deformed');
subplot(2,2,2); gpp_plot_mesh(source_face, source_vertex, correspondence_mask); colormap parula; title('Correspondence region');

target_vertex_3D = [target_vertex, -1*ones(size(target_vertex,1))];
source_vertex_3D = [source_vertex, zeros(size(source_vertex,1))];
displace = [displace_2D, -1*ones(size(displace_2D, 1))];
figure(4);
% subplot(1,2,2); 
hold on;
gpp_plot_mesh(target_face, target_vertex_3D, target_intensity); colormap parula; 
gpp_plot_mesh(source_face_reg, source_vertex_3D, source_intensity, 'EdgeAlpha', 0.2); colormap parula;
sub_ind = source_vertex_reg_intersect_index(1:37:end);
quiver3(source_vertex(sub_ind,1), source_vertex(sub_ind,2), source_vertex_3D(sub_ind,3),...
    displace_2D(sub_ind,1), displace_2D(sub_ind,2), displace(sub_ind, 3), 0, 'Color',...
    [216,108,112]/255, 'LineWidth', 0.2,...
    'MaxHeadSize', 0.05); colorbar;
hold off; title('Displacement field visualization');
figure(5);
% subplot(1,2,1); 
gpp_plot_mesh(source_face_reg, source_vertex, intensity_diff); colormap parula; title('Registered w/ curvature diff');

                
                
end