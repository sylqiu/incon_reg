function [face, vertex, gauss_curvature, mean_curvature] = load_shape_obj(obj_path)
disp('loading obj with format from Gabriel Peyr')
[vertex, face] = readOBJ(obj_path);
[face, vertex] = clean_mesh(face, vertex);
[~, ~, ~, ~, gauss_curvature, mean_curvature,~] = compute_curvature(vertex, face);
end