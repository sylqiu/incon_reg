function gernus_one_result_show_3D(source_face, source_vertex, flat_source_vertex, source_vertex_reg,...
                    target_face, target_vertex, flat_target_vertex,...
                    source_intensity, target_intensity,...
                    landmark_source_index, landmark_target_index,...
                    source_vertex_reg_intersect_index, correspondence_mask, target_correspondence_mask, ...
                    source_vertex_reg_3D, source_face_reg, displace_3D, dist_3D,...
                    target_intensity_reg, intensity_diff,...
                    landmark_err, intensity_err,...
                source_inner_boundary_index, target_inner_boundary_index)
markersize = 15;
%% plotting figures
figure(1);
subplot(2,2,1); gpp_plot_mesh(source_face, source_vertex, source_intensity); colormap parula; title('Source');
hold on; 
% plot3(source_vertex(landmark_source_index,1), source_vertex(landmark_source_index,2), source_vertex(landmark_source_index,3), 'yo');
plot3(source_vertex(source_inner_boundary_index,1), source_vertex(source_inner_boundary_index,2), source_vertex(source_inner_boundary_index,3), 'r');
hold off;
subplot(2,2,2); gpp_plot_mesh(target_face, target_vertex, target_intensity); colormap parula; title('Target');
hold on;
% plot3(target_vertex(landmark_target_index,1), target_vertex(landmark_target_index,2), target_vertex(landmark_target_index,3), 'ro'); 
plot3(target_vertex(target_inner_boundary_index,1), target_vertex(target_inner_boundary_index,2), target_vertex(target_inner_boundary_index,3), 'r');
hold off;
subplot(2,2,3); gpp_plot_mesh(source_face, flat_source_vertex); colormap parula; title('Flattened source');
hold on; 
plot(flat_source_vertex(landmark_source_index,1), flat_source_vertex(landmark_source_index,2), 'b.', 'MarkerSize', markersize); hold off;
subplot(2,2,4); gpp_plot_mesh(target_face, flat_target_vertex); colormap parula; title('Flattened target');
hold on; 
plot(flat_target_vertex(landmark_target_index,1), flat_target_vertex(landmark_target_index,2), 'r.', 'MarkerSize', markersize); hold off;

figure(2);
gpp_plot_mesh(source_face, source_vertex_reg, 'FaceColor', [36 169 225]/255, 'EdgeColor', [90 169 225]/255, 'FaceAlpha', 0.2); hold on;
gpp_plot_mesh(target_face, flat_target_vertex, 'FaceColor', [216,108,112]/255, 'EdgeColor', [216,108,112]/255, 'FaceAlpha', 0.2, 'EdgeAlpha', 0.5);
plot(source_vertex_reg(landmark_source_index,1), source_vertex_reg(landmark_source_index,2), 'b.','MarkerSize', markersize);
plot(flat_target_vertex(landmark_target_index,1), flat_target_vertex(landmark_target_index,2), 'r.', 'MarkerSize', markersize);
title('Overlaid view');
hold off;

figure(3);
% diff_ = abs(source_intensity - target_intensity_reg);
% diff_ = dist_3D;
% subplot(1,2,1); gpp_plot_mesh(source_face_reg, source_vertex, diff_); colormap parula; title('Registered w/ displacement magnitude');
% subplot(2,2,4); gpp_plot_mesh(target_face, target_vertex, target_intensity); colormap parula; title('Target');
% subplot(2,2,3); gpp_plot_mesh(source_face, source_vertex, source_intensity); colormap parula; title('Source');
subplot(1,2,1); gpp_plot_mesh(target_face, target_vertex, target_correspondence_mask*0.4+0.1);colormap parula; caxis([0,1]); title('Target correspondence region');
hold on;
plot3(target_vertex(target_inner_boundary_index,1), target_vertex(target_inner_boundary_index,2), target_vertex(target_inner_boundary_index,3), 'r');
hold off;
subplot(1,2,2); gpp_plot_mesh(source_face, source_vertex, correspondence_mask*0.4+0.1); colormap parula; caxis([0,1]); title('Correspondence region');
hold on;
plot3(source_vertex(source_inner_boundary_index,1), source_vertex(source_inner_boundary_index,2), source_vertex(source_inner_boundary_index,3), 'r');
hold off;
figure(4);
% subplot(1,2,2); 
hold on;
gpp_plot_mesh(target_face, target_vertex, target_intensity, 'EdgeAlpha', 0.2); colormap parula; 
gpp_plot_mesh(source_face_reg, source_vertex, source_intensity); colormap parula;
sub_ind = source_vertex_reg_intersect_index(1:18:end);
quiver3(source_vertex(sub_ind,1), source_vertex(sub_ind,2), source_vertex(sub_ind,3),...
    displace_3D(sub_ind,1), displace_3D(sub_ind,2), displace_3D(sub_ind,3),0, 'Color',...
    [216,108,112]/255, 'LineWidth', 0.2,...
    'MaxHeadSize', 0.05); colorbar;
hold off; title('Displacement field visualization');
figure(5);
% subplot(1,2,1); 
gpp_plot_mesh(source_face_reg, source_vertex, intensity_diff); colormap parula; caxis([0,1]); %title('Registered w/ intensity diff'); colorbar;

figure(6);
subplot(1,2,1); plot(1:length(landmark_err), landmark_err); %title('Landmark Error');
xlabel('Iter') 
ylabel('Landmark Error') 
subplot(1,2,2); plot(1:length(intensity_err), intensity_err); %title('Energy');
xlabel('Iter') 
ylabel('Energy') 
figure(7);
subplot(1,2,1);  hold on;
gpp_plot_mesh(source_face, source_vertex_reg, source_intensity);
plot(source_vertex_reg(landmark_source_index,1), source_vertex_reg(landmark_source_index,2), 'b.','MarkerSize',markersize);
hold off;
subplot(1,2,2); hold on; 
gpp_plot_mesh(target_face, flat_target_vertex, target_intensity);
plot(flat_target_vertex(landmark_target_index,1), flat_target_vertex(landmark_target_index,2), 'r.', 'MarkerSize',markersize);
hold off;


end