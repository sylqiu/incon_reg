function [source_face, source_vertex, flat_source_vertex, source_intensity,...
    target_face, target_vertex, flat_target_vertex, target_intensity,...
    landmark_source_index, landmark_target_index, landmark_target_pos,...
    L, V2Fm, F2Vm,...
    source_outer_boundary_index, target_outer_boundary_index, source_inner_boundary_index,...
    target_inner_boundary_index,...
    source_corner_id, target_corner_id] = hole_bone_mesh_preprocess(path_base, source_name, target_name, save_flag, varargin)


source_mat = load([path_base source_name '.mat']);
target_mat = load([path_base target_name '.mat']);

source_face = source_mat.holed_face;
source_vertex = source_mat.cut_vertex;
% [source_face, source_vertex, source_find] = clean_mesh(source_face, source_vertex);
% [~, ~, ~, ~, source_gauss_curvature, source_mean_curvature,~] = compute_curvature(source_vertex, source_face);
% source_gauss_curvature = source_mat.gauss_curvature(source_find);
% source_mean_curvature = source_mat.mean_curvature(source_find);
source_gauss_curvature = source_mat.gauss_curvature;
source_mean_curvature = source_mat.mean_curvature;


target_face = target_mat.holed_face;
target_vertex = target_mat.cut_vertex;
% [target_face, target_vertex, target_find] = clean_mesh(target_face, target_vertex);
% [~, ~, ~, ~, target_gauss_curvature, target_mean_curvature,~] = compute_curvature(target_vertex, target_face);landmark_target_index = [];
% target_gauss_curvature = target_mat.gauss_curvature(target_find);
% target_mean_curvature = target_mat.mean_curvature(target_find);
target_gauss_curvature = target_mat.gauss_curvature;
target_mean_curvature = target_mat.mean_curvature;

% landmark_source_index = [4362;  6;   2036; 2508; 3718];
% landmark_target_index = [11776; 936; 8886; 9378; 2135];

[target_intensity, min_, max_] = intensity_normalization(target_gauss_curvature);
source_intensity = intensity_normalization(source_gauss_curvature, min_, max_);


% flat_source_vertex = source_mat.flat_vertex(source_find,:);
% flat_target_vertex = target_mat.flat_vertex2(target_find,:);
flat_source_vertex = source_mat.flat_vertex;
flat_target_vertex = target_mat.flat_vertex;

%
if size(varargin) == 0
    disp('Choose landmark by hand \n');
    [landmark_source_index, landmark_target_index] = select_landmark(source_face, source_vertex,...
                                                    target_face, target_vertex, ...
                                                     source_intensity, target_intensity);
else 
    assert(length(varargin) == 2, ...
        'You should give source landmark indices and target landmark vertex indices as 1D lists');
    landmark_source_index = varargin{1};
    landmark_target_index = varargin{2};
end

assert(length(landmark_source_index) == length(landmark_target_index), ...
    'number of landmarks must equal');
%

L = cotmatrix(flat_source_vertex, source_face);
F2Vm = F2V(source_vertex', source_face');
V2Fm = V2F(source_vertex', source_face');
landmark_target_pos = flat_target_vertex(landmark_target_index,1:2);

source_outer_boundary_index = source_mat.total_boundary{1};
source_inner_boundary_index = source_mat.total_boundary{2};
% source_outer_boundary_index = find(ismember(source_find, source_outer_boundary_index));
% source_inner_boundary_index = find(ismember(source_find, source_inner_boundary_index));
% source_outer_boundary_index = orient_boundary(source_face, source_outer_boundary_index);
% source_inner_boundary_index = orient_boundary(source_face, source_inner_boundary_index);

target_outer_boundary_index = target_mat.total_boundary{1};
target_inner_boundary_index = target_mat.total_boundary{2};
% target_outer_boundary_index = find(ismember(target_find, target_outer_boundary_index));
% target_inner_boundary_index = find(ismember(target_find, target_inner_boundary_index));
% target_outer_boundary_index = orient_boundary(target_face, target_outer_boundary_index);
% target_inner_boundary_index = orient_boundary(target_face, target_inner_boundary_index);

source_corner_id = source_mat.new_corner_id;
% source_corner_id = find(ismember(source_find, source_corner_id));
target_corner_id = target_mat.new_corner_id;
% target_corner_id = find(ismember(target_find, target_corner_id));

figure(1);
subplot(2,2,1); gpp_plot_mesh(source_face, source_vertex, source_intensity); colormap parula; title('Source');
hold on; plot3(source_vertex(landmark_source_index,1), source_vertex(landmark_source_index,2), source_vertex(landmark_source_index,3), 'ro'); hold off;
subplot(2,2,2); gpp_plot_mesh(target_face, target_vertex, target_intensity); colormap parula; title('Target');
hold on; plot3(target_vertex(landmark_target_index,1), target_vertex(landmark_target_index,2), target_vertex(landmark_target_index,3), 'ro'); hold off;
subplot(2,2,3); gpp_plot_mesh(source_face, flat_source_vertex, source_intensity); colormap parula; title('Flattened source');
hold on; plot(flat_source_vertex(landmark_source_index,1), flat_source_vertex(landmark_source_index,2), 'ro'); hold off;
subplot(2,2,4); gpp_plot_mesh(target_face, flat_target_vertex, target_intensity); colormap parula; title('Flattened target');
hold on; plot(flat_target_vertex(landmark_target_index,1), flat_target_vertex(landmark_target_index,2), 'ro'); hold off;

if save_flag == 1
    file_name = sprintf('%s_%s_%s_workspace.mat', datestr(now,'mm-dd-yyyy HH-MM'), source_name, target_name);
    save(file_name);
    fprintf('Workspaced saved as %s \n', file_name);
end
end
