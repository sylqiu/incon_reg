clear;close all; clc;
load('01-14-2020 11-15_MN-Rio-7724_M1223_aligned.mat_AMNH-M-38792_M904_aligned.mat_results_3landmark.mat');
%%
[error,Reallignedsource,transform]=rigidICP(target_vertex,source_vertex,0,[],[]);
%%
figure; hold on; gpp_plot_mesh(source_face, Reallignedsource, 'EdgeColor', [36 169 225]/255);
gpp_plot_mesh(target_face, target_vertex, 'EdgeColor', [216,108,112]/255); title('Rigid ICP');
%%
[registered,~,~]=nonrigidICPv1(target_vertex,source_vertex,target_face,source_face,10,1);
%%
figure; hold on; gpp_plot_mesh(source_face, registered, 'EdgeColor', [36 169 225]/255);
gpp_plot_mesh(target_face, target_vertex, 'EdgeColor', [216,108,112]/255); title('Non-rigid ICP');
%%
close all; clc;
result_show_3D(source_face, source_vertex, flat_source_vertex, source_vertex_reg,...
                    target_face, target_vertex, flat_target_vertex,...
                    source_intensity, target_intensity,...
                    landmark_source_index, landmark_target_index,...
                    source_vertex_reg_intersect_index, correspondence_mask, ...
                    source_vertex_reg_3D, source_face_reg, displace, dist,...
                    target_intensity_reg, intensity_diff,...
                    landmark_err, intensity_err);
%%
figure; hold on; gpp_plot_mesh(source_face, source_vertex_reg_3D, 'EdgeColor', [36 169 225]/255);
gpp_plot_mesh(target_face, target_vertex, 'EdgeColor', [216,108,112]/255); title('Ours');