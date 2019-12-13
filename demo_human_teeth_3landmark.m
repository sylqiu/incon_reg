clc; close all;
% prefix = './data/';
% source_name = 'test_right_0.obj';
% target_name = 'test_left_0.obj';
% 
% source_obj_path = [prefix source_name];
% target_obj_path = [prefix target_name];
% save_flag = 1;
% 
% [source_face, source_vertex, flat_source_vertex, source_intensity,...
%     target_face, target_vertex, flat_target_vertex, target_intensity,...
%     landmark_source_index, landmark_target_index, landmark_target_pos,...
%     L, V2Fm, F2Vm,...
%     source_boundary_index, target_boundary_index] =...
%     data_mesh_preprocess(source_obj_path, target_obj_path,...
%                             source_name, target_name, save_flag);
load('12-09-2019 22-46_test_right_0.obj_test_left_0.obj_workspace_5landmark.mat');
sub_landmark = [1,3,5];
landmark_source_index = landmark_source_index(sub_landmark);
landmark_target_index = landmark_target_index(sub_landmark);
landmark_target_pos = landmark_target_pos(sub_landmark, :);
initial_landmark_error = compute_landamrk_err(flat_source_vertex, landmark_source_index,...
    landmark_target_pos);
[~,~,source_intensity_grid, target_intensity_grid, ~,~] = combine_to_same_grid(...
                                        flat_source_vertex, source_intensity,...
                                        flat_target_vertex, target_intensity,...
                                        source_boundary_index, target_boundary_index);
initial_intensity_error = sum(abs(source_intensity_grid(:) - target_intensity_grid(:)));
%
clc;
param.UpperBound = 1.5;
param.LowerBound = 0.8;
param.alpha = 0.01;
param.beta = 0.1;
param.smooth_iter = 3;
param.intensity_iter = 1;
param.demons_iter = 1; % turn higher if you want stronger intensity matching
param.demons_stepsize = 5;% turn higher if you want stronger intensity matching
param.landmark_iter = 1;
param.overall_iter = 20;
% algo begins
source_vertex_reg_pre = flat_source_vertex;
landmark_err = [initial_landmark_error];
intensity_err = [initial_intensity_error];
fprintf('L1 norm intensity difference % f ', initial_intensity_error);
fprintf('Eulidean landmark difference % f \n', initial_landmark_error);
for iter = 1:param.overall_iter
    fprintf('Iter %d \n', iter);
    for l_iter = 1:param.landmark_iter
        source_vertex_reg_pre = reg_landmark(source_face, flat_source_vertex, source_vertex_reg_pre,...
                                    landmark_source_index, landmark_target_pos,...
                                    source_boundary_index, F2Vm, V2Fm, L, param);
        figure(21); gpp_plot_mesh(source_face, source_vertex_reg_pre); title('Landmark matching step');
        drawnow;
    end
    [source_vertex_reg, ie] = reg_intensity(source_face, flat_source_vertex, source_vertex_reg_pre,...
                                source_intensity, flat_target_vertex, target_intensity,...
                                source_boundary_index, target_boundary_index,...
                                F2Vm, V2Fm, L, param);
    figure(31); gpp_plot_mesh(source_face, source_vertex_reg); title('Intensity matching step');
    drawnow;
    % calculate the landmark and intensity difference
    le = compute_landamrk_err(source_vertex_reg, ...
                        landmark_source_index,...
                        landmark_target_pos);
    fprintf('L1 norm intensity difference % f ', ie);
    fprintf('Eulidean landmark difference % f \n', le);
    landmark_err = [landmark_err, le];
    intensity_err = [intensity_err, ie];
    %
    source_vertex_reg_pre = source_vertex_reg;
    
end
% algo ends
%%
[source_vertex_reg_intersect_index, correspondence_mask, ...
    source_vertex_reg_3D, source_face_reg, displace, dist,...
    target_intensity_reg, intensity_diff] = ...
    prepare_result(source_face, source_vertex, source_vertex_reg,...
                    target_face, target_vertex, flat_target_vertex,...
                    source_intensity, target_intensity);
if save_flag == 1
    file_name = sprintf('%s_%s_%s_results_3landmark.mat', datestr(now,'mm-dd-yyyy HH-MM'), source_name, target_name);
    save(file_name);
    fprintf('Results saved !\n')
end
%% show results 
clc; close all;
% use this one obtained previously if you do not want to run previous part 
clear; 
load('12-13-2019 12-25_test_right_0.obj_test_left_0.obj_results_3landmark.mat');
%
result_show_3D(source_face, source_vertex, flat_source_vertex, source_vertex_reg,...
                    target_face, target_vertex, flat_target_vertex,...
                    source_intensity, target_intensity,...
                    landmark_source_index, landmark_target_index,...
                    source_vertex_reg_intersect_index, correspondence_mask, ...
                    source_vertex_reg_3D, source_face_reg, displace, dist,...
                    target_intensity_reg, intensity_diff,...
                    landmark_err, intensity_err);

