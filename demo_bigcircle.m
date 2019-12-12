clc; close all; clear;
load('data/QCIR_Bigcircle_example.mat');
save_flag = 0;
%
[source_face, source_vertex, source_intensity] = get_mesh_from_image(moving, 100);
[target_face, target_vertex, target_intensity] = get_mesh_from_image(static, 100);
flat_source_vertex = source_vertex;
flat_target_vertex = target_vertex;
landmark = imgcoord(moving, landmark) / 200;
target = imgcoord(static, target)/ 200;
landmark_source_index = findclstpt(source_vertex, landmark)';
landmark_target_index = findclstpt(target_vertex, target)';
num_landmark = length(landmark_source_index);
fprintf('Number of landmarks %d', num_landmark);
landmark_target_pos = flat_target_vertex(landmark_target_index,1:2);
% build laplacian of the source domain for coefficient smoothing
L = cotmatrix(source_vertex, source_face);
F2Vm = F2V(source_vertex', source_face');
V2Fm = V2F(source_vertex', source_face');
% assume one boundary component                           
[~, source_boundary_index] = meshboundaries(source_face);
[~, target_boundary_index] = meshboundaries(target_face);
% save workspace for later use
figure(1);
subplot(1,2,1); gpp_plot_mesh(source_face, source_vertex, source_intensity); colormap parula; title('Source');
hold on; plot(source_vertex(landmark_source_index,1), source_vertex(landmark_source_index,2), 'ro'); hold off;
subplot(1,2,2); gpp_plot_mesh(target_face, target_vertex, target_intensity); colormap parula; title('Target');
hold on; plot(target_vertex(landmark_target_index,1), target_vertex(landmark_target_index,2), 'ro'); hold off;
if save_flag == 1
    file_name = sprintf('%s_%s_%s_workspace.mat', datestr(now,'mm-dd-yyyy HH-MM'), source_name, target_name);
    save(file_name);
    fprintf('Workspaced saved as %s \n', file_name);
end
%%
clc;
param.UpperBound = 1.5;
param.LowerBound = 0.8;
param.alpha = 0.01;
param.beta = 0.1;
param.smooth_iter = 5;
param.intensity_iter = 1;
param.demons_iter = 2;
param.demons_stepsize = 1;
param.landmark_iter = 1;
param.overall_iter = 10;
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
[source_vertex_reg_intersect_index, correspondence_mask, ...
    source_vertex_reg_3D, source_face_reg, displace, dist,...
    target_intensity_reg, intensity_diff] = ...
    prepare_result(source_face, source_vertex, source_vertex_reg,...
                    target_face, target_vertex, flat_target_vertex,...
                    source_intensity, target_intensity);
if save_flag == 1
    file_name = sprintf('%s_%s_%s_results.mat', datestr(now,'mm-dd-yyyy HH-MM'), source_name, target_name);
    save(file_name);
    fprintf('Results saved !\n')
end
%%
close all; clc;
load('12-12-2019 12-05_brain_T1_wave_part.png_brain_T1_part.png_results.mat');
result_show_image(source_face, source_vertex, flat_source_vertex, source_vertex_reg,...
                    target_face, target_vertex, flat_target_vertex,...
                    source_intensity, target_intensity,...
                    landmark_source_index, landmark_target_index,...
                    source_vertex_reg_intersect_index, correspondence_mask, ...
                    source_vertex_reg_3D, source_face_reg, displace, dist,...
                    target_intensity_reg, intensity_diff);