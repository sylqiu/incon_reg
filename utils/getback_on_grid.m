function [Mg,mask] = getback_on_grid(vertex1,V_t,Mv,g3,g4,correct_bdym,correct_bdyt)
% interpolating to [g3,g4]
Ftmp = delaunay(g3(:),g4(:));
L = cotmatrix([g3(:),g4(:)],Ftmp);
intpM = scatteredInterpolant(vertex1(:,1),vertex1(:,2),Mv);
Mg = intpM(g3,g4);
bdy_xym = vertex1(correct_bdym,:);
bdy_xyt = V_t(correct_bdyt,:);
bout_m = ~in_or_out(bdy_xym,[g3(:),g4(:)]); bout_m = reshape(bout_m,size(g3,1),size(g3,2));
bout_t = ~in_or_out(bdy_xyt,[g3(:),g4(:)]); bout_t = reshape(bout_t,size(g3,1),size(g3,2));
mask = isnan(Mg) | bout_m | bout_t;
Mg(mask) = 0; 
mask = ~mask(:)*1.0;
% for kkk = 1:100
%     mask = mask + 0.2*L*mask;
% end
mask = reshape(mask,size(g3,1),size(g3,2)); 
% Mg = Mg.*mask;
end