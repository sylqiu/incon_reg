function [A, b] = genus_one_beltrami_operator(face, vertex, mu, boundary_id, corner_id)
boundary_id = reshape(boundary_id, [], 1);
corner_id = reshape(corner_id, [], 1);
%% parameterization portion
aaf = -(1-2*real(mu)+abs(mu).^2)./(1.0-abs(mu).^2);
bbf = 2*imag(mu)./(1.0-abs(mu).^2);
ggf = -(1+2*real(mu)+abs(mu).^2)./(1.0-abs(mu).^2);

aaf(abs(mu)>1) = -aaf(abs(mu)>1);
bbf(abs(mu)>1) = -bbf(abs(mu)>1);
ggf(abs(mu)>1) = -ggf(abs(mu)>1);

af = aaf;
bf = bbf;
gf = ggf;
f0 = face(:,1);
f1 = face(:,2);
f2 = face(:,3);

uxv0 = vertex(f1,2) - vertex(f2,2);
uyv0 = vertex(f2,1) - vertex(f1,1);
uxv1 = vertex(f2,2) - vertex(f0,2);
uyv1 = vertex(f0,1) - vertex(f2,1); 
uxv2 = vertex(f0,2) - vertex(f1,2);
uyv2 = vertex(f1,1) - vertex(f0,1);
%lengths
l = [sqrt(sum(uxv0.^2 + uyv0.^2,2)) ...
    sqrt(sum(uxv1.^2 + uyv1.^2,2)) ...
    sqrt(sum(uxv2.^2 + uyv2.^2,2))];
%half the perimeter of each triangle, denoted s
s = sum(l,2)*0.5;
%area of triangle with sides (a,b,c) = sqrt(s(s-a)(s-b)(s-c))
area = sqrt( s.*(s-l(:,1)).*(s-l(:,2)).*(s-l(:,3))).*2;
%
v00 = (af.*uxv0.*uxv0 + 2*bf.*uxv0.*uyv0 + gf.*uyv0.*uyv0)./area;
v11 = (af.*uxv1.*uxv1 + 2*bf.*uxv1.*uyv1 + gf.*uyv1.*uyv1)./area;
v22 = (af.*uxv2.*uxv2 + 2*bf.*uxv2.*uyv2 + gf.*uyv2.*uyv2)./area;

v01 = (af.*uxv1.*uxv0 + bf.*uxv1.*uyv0 + bf.*uxv0.*uyv1 + gf.*uyv1.*uyv0)./area;
v12 = (af.*uxv2.*uxv1 + bf.*uxv2.*uyv1 + bf.*uxv1.*uyv2 + gf.*uyv2.*uyv1)./area;
v20 = (af.*uxv0.*uxv2 + bf.*uxv0.*uyv2 + bf.*uxv2.*uyv0 + gf.*uyv0.*uyv2)./area;

I = [f0;f1;f2;f0;f1;f1;f2;f2;f0];
J = [f0;f1;f2;f1;f0;f2;f1;f0;f2];
V = [v00;v11;v22;v01;v01;v12;v12;v20;v20];
A = sparse(I,J,V);
A0 = A;
%% build rhs
b = zeros(size(vertex,1), 2);
%% boundary condition
M_adj = adjacency_matrix(face);
%
corner_offset = ismember(boundary_id,corner_id);
tmp = 1:size(boundary_id,1);
corner_offset = tmp(corner_offset);
boundary_id = circshift(boundary_id, -corner_offset(1)+1, 1);
corner_offset = corner_offset - corner_offset(1) + 1;
for i = 1:4
    
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
        tmp = co_start_id : -1 : co_start_id - length(boundary_segment) + 1;
        tmp(tmp<=0) = tmp(tmp<=0) + size(boundary_id,1);
        co_boundary_segment = boundary_id(tmp);
    else
        co_boundary_segment = boundary_id(co_start_id :-1: co_end_id);
    end
    
    basis = vertex(boundary_id(start_id),:) - vertex(boundary_id(co_start_id),:);
    
%     figure(1); 
%     xlim([-1 1])
%     ylim([-1 1]); hold on;
%     % ind = boundary_id;
%     ind = [boundary_segment; co_boundary_segment];
%     for ii = 1:size(ind,1)
%        plot(vertex(ind(ii),1), vertex(ind(ii),2), 'o');
%        drawnow;
%     end
    
%     figure(2); hold on; gpp_plot_mesh(face,vertex);
    for j = 2:size(boundary_segment,1)-1
%         neighbors = find(M_adj(:, boundary_segment(j)) == 1);
%         plot(vertex(boundary_segment(j), 1), vertex(boundary_segment(j),2),'ro');
        co_neighbors = find(M_adj(:, co_boundary_segment(j)) == 1);
        for k = 1:size(co_neighbors)
            
            if co_neighbors(k) ~= co_boundary_segment(j-1) && co_neighbors(k) ~= co_boundary_segment(j+1) 
                
%                 plot(vertex(co_neighbors(k), 1), vertex(co_neighbors(k),2),'go');
                
                w = A0(co_boundary_segment(j), co_neighbors(k));
                A(boundary_segment(j), co_neighbors(k)) = w;
%                 A(co_neighbors(k), boundary_segment(j)) = w;
                b(boundary_segment(j),:) = b(boundary_segment(j),:) - w*basis;
           
            end
        end
        A(boundary_segment(j), boundary_segment(j)) = A0(boundary_segment(j), boundary_segment(j))+A0(co_boundary_segment(j), co_boundary_segment(j));
        A(boundary_segment(j), boundary_segment(j-1)) = A0(boundary_segment(j), boundary_segment(j-1))+A0(co_boundary_segment(j), co_boundary_segment(j-1));
        A(boundary_segment(j), boundary_segment(j+1)) = A0(boundary_segment(j), boundary_segment(j+1))+A0(co_boundary_segment(j), co_boundary_segment(j+1));
        
    end
    
    
end


end