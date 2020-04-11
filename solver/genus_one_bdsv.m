function [map, se] = genus_one_bdsv(Vpre,V,F,b,bc,outer_boundary_id, corner_id, K1,K2)
% QC map from bounds on singular values 0<K2 <= s1,s2 <= K1
% K2 = K1 = 1 equivalent to a version of ARAP
% Initiallize with (V_pre,F), gradient operator of (V,F)
% V_pre, V: #V x 3
% F: #F x 3
% b: #b x1 constraint vertex indices
% bc: #b x 3 target positions

%  Mar 21. 2020 @sylvesterqiu@gamil.com


corner_offset = ismember(outer_boundary_id,corner_id);
tmp = 1:size(outer_boundary_id,1);
corner_offset = tmp(corner_offset);
boundary_id = circshift(outer_boundary_id, -corner_offset(1)+1, 1);
corner_offset = corner_offset - corner_offset(1) + 1;
% get redundant vertices and its corresponding vertices id and shift
% pattern is like following, suppose goes counterclockwise, denote 0 as redundant
% 0 0 0 0 0
% 1       0
% 1       0
% 1       0
% 1 1 1 1 0
redundant_boundary_id = [];
corresponding_boundary_id = [];
redundant_boundary_shift = [];
for i = 1:2
    
    next_i = mod(i, 4)+1;
    start_id = mod(corner_offset(i)-1, size(boundary_id,1))+1;
    end_id = mod(corner_offset(next_i)-1, size(boundary_id,1))+1;
    
    if end_id < start_id
        end_id = end_id + size(boundary_id,1);
        tmp = start_id : end_id;
        tmp(tmp>size(boundary_id,1)) = tmp(tmp>size(boundary_id,1)) - size(boundary_id,1);
        boundary_segment =  boundary_id(tmp);
    else
        boundary_segment =  boundary_id(start_id : end_id);
    end
    

    
    co_i = mod(i+2, 4)+1;
    co_next_i = mod(i+1, 4)+1;
    co_start_id = mod(corner_offset(co_i)-1, size(boundary_id,1))+1;
    co_end_id = mod(corner_offset(co_next_i)-1, size(boundary_id,1))+1;
    if co_start_id < co_end_id
        tmp = co_start_id : -1 : co_start_id - length(boundary_segment)+1;
        tmp(tmp<=0) = tmp(tmp<=0) + size(boundary_id,1);
        co_boundary_segment = boundary_id(tmp);
    else
        co_boundary_segment = boundary_id(co_start_id :-1: co_end_id);
    end
    
    
    % gather corner redundant firtst
    if i == 1
       redundant_boundary_id = [redundant_boundary_id; boundary_segment(1); co_boundary_segment(1); co_boundary_segment(end)]; 
       corresponding_boundary_id = [corresponding_boundary_id; boundary_segment(end);  boundary_segment(end); boundary_segment(end)]; 
       redundant_boundary_shift = [redundant_boundary_shift;... 
           V(boundary_segment(1), :) - V(boundary_segment(end), :); ...
           V(co_boundary_segment(1), :) - V(boundary_segment(end), :); ...
           V(co_boundary_segment(end), :) - V(boundary_segment(end), :)];
    end
    
    redundant_boundary_id = [redundant_boundary_id; co_boundary_segment(2:end-1)];
    corresponding_boundary_id = [corresponding_boundary_id; boundary_segment(2:end-1)];
    redundant_boundary_shift = [redundant_boundary_shift; ...
        repmat(V(co_boundary_segment(end), :) - V(boundary_segment(end), :), size(boundary_segment,1)-2,1)];
    

end
solve_vertex_find = setdiff((1:size(V,1))', redundant_boundary_id); % new vertices indexing trace back
Vpre_solve = Vpre(solve_vertex_find,:);
% reset constraint index
% fix one corner point
b = reshape(b, [], 1);
b = [b; boundary_id(corner_offset(2))];
bc = [bc; V(boundary_id(corner_offset(2)), :)];

% sanity check no other boundary point prescription
[b, tmp] = setdiff(b, setdiff(boundary_id, boundary_id(corner_offset(2))));
bc_solve = bc(tmp,:);

b_solve = ismember(solve_vertex_find, b);
tmp = 1:size(solve_vertex_find,1);
b_solve = tmp(b_solve);

% solve_vertex_id = (1:size(solve_vertex_find,1))';

% build gradient operator
% grad(Xijk) = Xj * (Vi - Vk)^R90 / 2A + Xk * (Vj - Vi)^R90 / 2A + 
%             -Xi * (Vi - Vk)^R90 / 2A - Xi * (Vj - Vi)^R90 / 2A
% where Xi is the scalar value at vertex i, Vi is the 3D position of vertex
% i, and A is the area of triangle (i,j,k). ^R90 represent a rotation of 
% 90 degrees
% renaming indices of vertices of triangles for convenience
i1 = F(:,1); i2 = F(:,2); i3 = F(:,3); 
V = [V, zeros(size(V,1),1)];
% #F x 3 matrices of triangle edge vectors, named after opposite vertices
v32 = V(i3,:) - V(i2,:);  v13 = V(i1,:) - V(i3,:); v21 = V(i2,:) - V(i1,:);
% area of parallelogram is twice area of triangle
% area of parallelogram is || v1 x v2 || 
n  = cross(v32,v13,2); 
% This does correct l2 norm of rows, so that it contains #F list of twice
% triangle areas
dblA = normrow(n);
% now normalize normals to get unit normals
u = normalizerow(n);
eperp21 = bsxfun(@times,cross(u,v21),1./dblA);
eperp13 = bsxfun(@times,cross(u,v13),1./dblA);

G_pre = zeros(2*size(F,1), 3);
II = zeros(2*size(F,1), 3);
JJ = zeros(2*size(F,1), 3);
b_pre = zeros(2*size(F,1), 2);
num_F = size(F,1);

% set boundary condition for grad operator
for i = 1:size(F,1)
    G_pre(i,:) = [-eperp13(i, 1)-eperp21(i,1), eperp13(i, 1), eperp21(i,1)];
    G_pre(i+num_F, :) = [-eperp13(i, 2)-eperp21(i,2), eperp13(i, 2), eperp21(i,2)];
    II(i,:) = [i, i, i];
    II(i+num_F,:) = [i+num_F,i+num_F,i+num_F];
    
    for j = 1:3        
        redundancy_check = ismember(redundant_boundary_id, F(i, j));
        if sum(redundancy_check(:)) > 0
            nFi_j = find(solve_vertex_find == corresponding_boundary_id(redundancy_check));
            b_pre(i,:) = b_pre(i,:) + G_pre(i, j) * redundant_boundary_shift(redundancy_check,:);
            b_pre(i+num_F,:) = b_pre(i+num_F,:) + G_pre(i+num_F, j) * redundant_boundary_shift(redundancy_check,:);
        else
            nFi_j = find(solve_vertex_find == F(i,j));            
        end
         JJ(i, j) = nFi_j;
         JJ(i+num_F,j) = nFi_j;
    end
    
end
G = sparse(II(:), JJ(:), G_pre(:));

%
Gu = G*Vpre_solve(:,1) + b_pre(:,1); 
Gv = G*Vpre_solve(:,2) + b_pre(:,2); 
%
Gv = reshape(Gv, num_F, 2);
Gu = reshape(Gu, num_F, 2);
S = [Gu,Gv]'; S = reshape(S,2,[]);
R = zeros(2,2*num_F);
se = 0;
for i = 1:num_F
    [uu,ss,vv] = svd(S(:,2*i-1 : 2*i));
    sv1 = min(ss(1,1),K1);
    sv2 = max(ss(2,2),K2);
    se = se + ((ss(1,1) - sv1)^2 + (ss(2,2) - sv2)^2) * face_area(F(i,:), V);
    R(:,2*i-1:2*i) = uu*diag([sv1,sv2])*vv';
end
R = reshape(R,4,[]); R = R';
R = [reshape(R(:,1:2),[],1),reshape(R(:,3:4),[],1)];
A = sparse(1:2*num_F, 1:2*num_F, repmat(dblA/2,2,1));
rhs = G'*A*(R - b_pre);
L = G'*A*G;
%

% 
if ~isempty(b_solve)
    rhs = -L(:,b_solve)*bc_solve + rhs;
    rhs(b_solve,:) = bc_solve;
    L(b_solve,:) = 0;
    L(:,b_solve) = 0;
    for ind = 1:length(b_solve)
        L(b_solve(ind),b_solve(ind)) = 1;
    end
end

map_solve = L\rhs;

map = zeros(size(V,1),2);
map(solve_vertex_find, :) = map_solve;
map(redundant_boundary_id,:) = map(corresponding_boundary_id, :) + redundant_boundary_shift;


end