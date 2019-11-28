function [face, vertex, gauss_curvature, mean_curvature] = load_shape_obj(obj_path)
disp('loading obj with format from Gabriel Peyr')
[vertex, face] = readOBJ(obj_path);
[face, vertex] = clean_mesh(face, vertex);
gauss_curvature = discrete_gaussian_curvature(vertex, face);
L = cotmatrix(vertex,face);
mean_curvature = L*vertex;
mean_curvature = real(sqrt(sum(mean_curvature.^2,2)));
end