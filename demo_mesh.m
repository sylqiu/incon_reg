%% read obj files
clear; clc; close all;
%%
% source_obj_path = 'MN-Rio-7724_M1223_aligned.obj';
% target_obj_path = 'AMNH-M-38792_M904_aligned.obj';
% save_flag = 1;
% 
% [source_face, source_vertex, flat_source_vertex, source_intensity,...
%     target_face, target_vertex, flat_target_vertex, target_intensity,...
%     landmark_source_index, landmark_target_index, landmark_target_pos,...
%     L, V2Fm, F2Vm,...
%     source_boundary_index, target_boundary_index] =...
%     data_mesh_preprocess(source_obj_path, target_obj_path, save_flag);
%%
% optionally you could load a previously saved workspace
load('11-28-2019 21-20_MN-Rio-7724_M1223_aligned.obj_AMNH-M-38792_M904_aligned.obj_workspace.mat');
%% perform registration
% set params
clc;
param.UpperBound = 1.5;
param.LowerBound = 0.8;
param.alpha = 0.01;
param.beta = 0.1;
param.smooth_iter = 3;
param.intensity_iter = 1;
param.demons_iter = 5;
param.demons_stepsize = 1;
param.landmark_iter = 1;
param.overall_iter = 30;
% algo begins
source_vertex_reg_pre = flat_source_vertex;
for iter = 1:param.overall_iter
    fprintf('Iter %d \n', iter);
    for l_iter = 1:param.landmark_iter
        source_vertex_reg_pre = reg_landmark(source_face, flat_source_vertex, source_vertex_reg_pre,...
                                    landmark_source_index, landmark_target_pos,...
                                    source_boundary_index, F2Vm, V2Fm, L, param);
        figure(21); gpp_plot_mesh(source_face, source_vertex_reg_pre); title('Landmark matching step');
        drawnow;
    end
    source_vertex_reg = reg_intensity(source_face, flat_source_vertex, source_vertex_reg_pre,...
                                source_intensity, flat_target_vertex, target_intensity,...
                                source_boundary_index, target_boundary_index,...
                                F2Vm, V2Fm, L, param);
    figure(31); gpp_plot_mesh(source_face, source_vertex_reg); title('Intensity matching step');
    drawnow;
    source_vertex_reg_pre = source_vertex_reg;
end
% algo ends
%%
result_show_3D(source_face, source_vertex, flat_source_vertex, source_vertex_reg,...
                    target_face, target_vertex, flat_target_vertex,...
                    source_intensity, target_intensity,...
                    landmark_source_index, landmark_target_index)
