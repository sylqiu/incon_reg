function data_image_preprocess(source_image, target_image, save_flag, varargin)

% create meshes for the images
disp('loading images ...');

[source_face, source_vertex, source_intensity] = get_mesh_from_image(source_image);
[target_face, target_vertex, target_intensity] = get_mesh_from_image(target_image);
disp('loading finished');

if size(varargin) == 0
    disp('Choose landmark by hand');
    [landmark_source_ind, landmark_target_ind] = select_landmark(source_face, source_vertex,...
                                                    target_face, target_vertex);
else 
    assert(size(varagin) == 2, 'You should give source landmarks and target landmarks as pixel positions');
    landmark_source_ind = varargin{1};
    landmark_target_ind = varargin{2};
end

% build laplacian of the source domain for coefficient smoothing
L = cotmatrix(source_vertex, source_face);
F2Vm = F2V(source_vertex', source_face');
V2Fm = V2F(source_vertex', source_face');
% assume one boundary component                           
[~, source_boundary_index] = meshboundaries(source_face);
[~, target_boundary_index] = meshboundaries(target_face);

% save workspace for later use


end