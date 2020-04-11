function map = genus_one_beltrami_solver(face, vertex, mu, corner_id, boundary_id, target_id, target_xy)
target_id = reshape(target_id, [], 1);
corner_id = reshape(corner_id, [], 1);
% build interior operator
[L, b] = genus_one_beltrami_operator(face, vertex, mu, boundary_id, corner_id);
% set constraint
target_id = [target_id; corner_id];
target_xy = [target_xy; vertex(corner_id,:)];
b = b - L(:, target_id)*target_xy;
b(target_id,:) = target_xy;
%
L(target_id, :) = 0;
L(:, target_id) = 0;
L = L + sparse(target_id, target_id, ones(size(target_id,1),1), size(L,1), size(L,1));
map = L\b;
% map = reshape(map, [], 2);
% figure; hold on;
% gpp_plot_mesh(face, map);
% plot(map(corner_id,1), map(corner_id,2), 'ro');
end