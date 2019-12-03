function [source_face, source_vertex, flat_source_vertex, source_intensity,...
    target_face, target_vertex, flat_target_vertex, target_intensity,...
    landmark_source_index, landmark_target_index, landmark_target_pos,...
    L, V2Fm, F2Vm,...
    source_boundary_index, target_boundary_index] = data_mesh_preprocess(source_obj, target_obj,...
    source_name, target_name,save_flag, varargin)

disp('loading obj files ...\n')
[source_face, source_vertex, ~, source_mean_curvature] = load_shape_obj(source_obj);
[target_face, target_vertex, ~, target_mean_curvature] = load_shape_obj(target_obj);
disp('loading finished, here using mean curvature as matching intensity, and are normalized to be in range [0,1] \n')
source_intensity = intensity_normalization(source_mean_curvature);
target_intensity = intensity_normalization(target_mean_curvature);
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

num_landmark = length(landmark_source_index);

disp('furtherst two points in the target will be used for conformal flattening \n');
distance_table = repmat(reshape(target_vertex(landmark_source_index, :), ...
    [1, num_landmark, 3]), [num_landmark, 1, 1]) - ...
    repmat(reshape(target_vertex(landmark_source_index, :), ...
    [num_landmark, 1, 3]), [1, num_landmark, 1]);
distance_table = sum(distance_table.^2, 3);
[~, sorted_] = sort(distance_table(:));
lscm_target_ind1 = mod(sorted_(end), num_landmark);
lscm_target_ind2 = round(sorted_(end)/num_landmark + 0.5);
lscm_target_ind = [landmark_target_index(lscm_target_ind1), landmark_target_index(lscm_target_ind2)];
lscm_source_ind = [landmark_source_index(lscm_target_ind1), landmark_source_index(lscm_target_ind2)];

% conformal flattening
flat_source_vertex = lscm(source_vertex, source_face, lscm_source_ind, [0, 0; 1, 0]);
flat_target_vertex = lscm(target_vertex, target_face, lscm_target_ind, [0, 0; 1, 0]);
disp('Conformal flattening done\n');
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