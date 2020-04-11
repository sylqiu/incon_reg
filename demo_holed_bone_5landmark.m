clc; close all; clear;
prefix = './data/';
source_name = 'holebone_source';
target_name = 'holebone_target';
save_flag = 0;
[source_face, source_vertex, flat_source_vertex, source_intensity,...
    target_face, target_vertex, flat_target_vertex, target_intensity,...
    landmark_source_index, landmark_target_index, landmark_target_pos,...
    L, V2Fm, F2Vm,...
    source_outer_boundary_index, target_outer_boundary_index, source_inner_boundary_index,...
    target_inner_boundary_index,...
    source_corner_index, target_corner_index] = hole_bone_mesh_preprocess(prefix, source_name, target_name, save_flag);

%%
[~,~,source_intensity_grid, target_intensity_grid, ~,~] = genus_one_combine_to_same_grid(...
                                        flat_source_vertex, source_intensity,...
                                        flat_target_vertex, target_intensity,...
                                        source_inner_boundary_index, target_inner_boundary_index,...
                                        source_outer_boundary_index,...
                                        source_corner_index);
initial_intensity_error = mean(abs(source_intensity_grid(:) - target_intensity_grid(:)).^2);
initial_landmark_error = compute_landamrk_err(flat_source_vertex, landmark_source_index,...
    landmark_target_pos);
%%
clc; warning('off','all')
param.UpperBound = 2;
param.LowerBound = 0.5;
param.alpha = 0.01;
param.beta = 0.1;
param.smooth_iter = 3;
param.intensity_iter = 1;
param.demons_iter = 1; % turn higher if you want stronger intensity matching
param.demons_stepsize = 2;% turn higher if you want stronger intensity matching
param.landmark_iter = 1;
param.overall_iter = 80;
% algo begins
source_vertex_reg_pre = flat_source_vertex;
% landmark_err = [initial_landmark_error];
% intensity_err = [initial_intensity_error];
landmark_err = [];
intensity_err = [];
fprintf('L1 norm intensity difference % f ', initial_intensity_error);
fprintf('Eulidean landmark difference % f \n', initial_landmark_error);
for iter = 1:param.overall_iter
    fprintf('Iter %d \n', iter);
    for l_iter = 1:param.landmark_iter
        tic;
        [source_vertex_reg_pre, ne] = genus_one_reg_landmark(source_face, flat_source_vertex, source_vertex_reg_pre,...
                                    landmark_source_index, landmark_target_pos,...
                                    source_outer_boundary_index, source_inner_boundary_index, source_corner_index, F2Vm, V2Fm, L, param);
%         figure(21); gpp_plot_mesh(source_face, source_vertex_reg_pre); title('Landmark matching step');
        toc;
        drawnow;
    end
    tic;
    [source_vertex_reg, ie] = genus_one_reg_intensity(source_face, flat_source_vertex, source_vertex_reg_pre,...
                                source_intensity, flat_target_vertex, target_intensity,...
                                source_inner_boundary_index, target_inner_boundary_index,...
                                source_outer_boundary_index,...
                                source_corner_index, param);
    toc;
%     figure(31); gpp_plot_mesh(source_face, source_vertex_reg); title('Intensity matching step');
    drawnow;
    % calculate the landmark and intensity difference
    le = compute_landamrk_err(source_vertex_reg, ...
                        landmark_source_index,...
                        landmark_target_pos);
    fprintf('L1 norm intensity difference % f ', ie);
    fprintf('Eulidean landmark difference % f \n', le);
    landmark_err = [landmark_err, le];
    intensity_err = [intensity_err, ie + ne];
    %
    source_vertex_reg_pre = source_vertex_reg;
    
end
% algo ends
%%
[source_vertex_reg_intersect_index, correspondence_mask, target_correspondence_mask,...
    source_vertex_reg_3D, source_face_reg, displace, dist,...
    target_intensity_reg, intensity_diff] = ...
    genus_one_prepare_result(source_face, source_vertex, source_vertex_reg,...
                    target_face, target_vertex, flat_target_vertex,...
                    target_outer_boundary_index, target_inner_boundary_index, target_corner_index,...
                    source_intensity, target_intensity, source_inner_boundary_index);
if save_flag == 1
    file_name = sprintf('%s_%s_%s_results_3landmark.mat', datestr(now,'mm-dd-yyyy HH-MM'), source_name, target_name);
    save(file_name);
    fprintf('Results saved !\n')
end
%% show results 
load('03-30-2020 16-55_holebone_source_holebone_target_results_3landmark.mat')
%%
clc; close all;
% use this one obtained previously if you do not want to run previous part 
%
genus_one_result_show_3D(source_face, source_vertex, flat_source_vertex, source_vertex_reg,...
                    target_face, target_vertex, flat_target_vertex,...
                    source_intensity, target_intensity,...
                    landmark_source_index, landmark_target_index,...
                    source_vertex_reg_intersect_index, correspondence_mask, target_correspondence_mask,...
                    source_vertex_reg_3D, source_face_reg, displace, dist,...
                    target_intensity_reg, intensity_diff,...
                    landmark_err, intensity_err,...
                    source_inner_boundary_index, target_inner_boundary_index);
%%
% param.position_cluster_strength = 0.5;
% param.intensity_cluster_strength = 10;
% param.num_cluster = 10;
% cluster_data = [displace, param.position_cluster_strength * flat_source_vertex,...
%     param.intensity_cluster_strength * source_intensity];
% rng default;
% [cluster_label, clusters] = kmeans(cluster_data, param.num_cluster);
% figure(7); hold on;
% gpp_plot_mesh(source_face, source_vertex, cluster_label); colorbar; caxis([1,param.num_cluster]);
% cluster_direction = clusters(:, 1:3);
% cluster_direction = cluster_direction ./ sqrt(sum(cluster_direction.^2, 2));
% % color_string_map = {'blue', 'shallow blue', 'green', 'yellow'};
% for cluster_id = 1:param.num_cluster
%     tmp_id = cluster_label == cluster_id;
%    
%     mean_dot_dist = mean(sum(cluster_data(tmp_id, 1:3).*cluster_direction(cluster_id,:),2));
%     fprintf('Cluster %d has mean dot-product distance %f \n', cluster_id, mean_dot_dist);
%     tmp_xyz = source_vertex(tmp_id, :);
%     tmp_xyz = tmp_xyz(1:57:end,:);
%     tmp_dir = mean_dot_dist * cluster_direction(cluster_id,:);
%     tmp_dir = repmat(tmp_dir, [size(tmp_xyz,1), 1]);
%     quiver3(tmp_xyz(:,1), tmp_xyz(:,2), tmp_xyz(:,3), tmp_dir(:,1),  tmp_dir(:,2), tmp_dir(:,3));
% end
% hold off


