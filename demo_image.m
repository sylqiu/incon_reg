clc; close all;
prefix = './data/';
source_name = 'brain_T1_wave_part.png';
target_name = 'brain_T1_part.png';


source_obj_path = [prefix source_name];
target_obj_path = [prefix target_name];
save_flag = 1;
landmark_source_index = [5556,1113,5946,4242,7293];
landmark_target_index = [10850,5591,10760,7040,10807];
[source_face, source_vertex, flat_source_vertex, source_intensity,...
    target_face, target_vertex, flat_target_vertex, target_intensity,...
    landmark_source_index, landmark_target_index, landmark_target_pos,...
    L, V2Fm, F2Vm,...
    source_boundary_index, target_boundary_index] =...
    data_image_preprocess(source_obj_path, target_obj_path,...
                            source_name, target_name, save_flag,...
                            landmark_source_index, landmark_target_index);
%%
clc;
param.UpperBound = 1.5;
param.LowerBound = 0.8;
param.alpha = 0.01;
param.beta = 0.1;
param.smooth_iter = 3;
param.intensity_iter = 1;
param.demons_iter = 2;
param.demons_stepsize = 5;
param.landmark_iter = 1;
% param.overall_iter = 30;
% algo begins
source_vertex_reg_pre = flat_source_vertex;
% source_vertex_reg_pre = reg_landmark(source_face, flat_source_vertex, source_vertex_reg_pre,...
%                                     landmark_source_index, landmark_target_pos,...
%                                     source_boundary_index, F2Vm, V2Fm, L, param);
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
close all
result_show_image(source_face, source_vertex, flat_source_vertex, source_vertex_reg,...
                    target_face, target_vertex, flat_target_vertex,...
                    source_intensity, target_intensity,...
                    landmark_source_index, landmark_target_index)