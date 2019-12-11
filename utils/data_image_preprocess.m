function [source_face, source_vertex, flat_source_vertex, source_intensity,...
    target_face, target_vertex, flat_target_vertex, target_intensity,...
    landmark_source_index, landmark_target_index, landmark_target_pos,...
    L, V2Fm, F2Vm,...
    source_boundary_index, target_boundary_index] = ...
    data_image_preprocess(source_image_path, target_image_path,...
                            source_name, target_name, save_flag, varargin)

% create meshes for the images
disp('loading images ...');

[source_face, source_vertex, source_intensity] = get_mesh_from_image(source_image_path, 100);
[target_face, target_vertex, target_intensity] = get_mesh_from_image(target_image_path, 100);
disp('loading finished');
flat_source_vertex = source_vertex;
flat_target_vertex = target_vertex;

if size(varargin) == 0
    disp('Choose landmark by hand');
    [landmark_source_index, landmark_target_index] = select_landmark(source_face, source_vertex,...
                                                    target_face, target_vertex, ...
                                                    source_intensity, target_intensity);
else 
    assert(length(varargin) == 2, 'You should give source landmarks and target landmarks');
    landmark_source_index = varargin{1};
    landmark_target_index = varargin{2};
end
assert(length(landmark_source_index) == length(landmark_target_index), ...
    'number of landmarks must equal');

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

end