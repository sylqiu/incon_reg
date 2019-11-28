function [Tx,Ty] =  intensity_fitting(M,S,Niter)
% assume S is already interpolated to the correct grid
alpha=3;
Tx=zeros(size(S)); Ty=zeros(size(S));
[Sx,Sy] = gradient(S);
Hsmooth=fspecial('gaussian',[60 60],10);
initial_M = M;
% imshow(M);
for itt=1:Niter
        % Difference image between moving and static image
        Idiff = (M-S);
        [Mx,My] = gradient(M);
        % Default demon force, (Thirion 1998)
        %Ux = -(Idiff.*Sx)./((Sx.^2+Sy.^2)+Idiff.^2);
        %Uy = -(Idiff.*Sy)./((Sx.^2+Sy.^2)+Idiff.^2);

        % Extended demon force. With forces from the gradients from both
        % moving as static image. (Cachier 1999, He Wang 2005)
%       figure(1);quiver(vertex2(:,1),vertex2(:,2),Mx,My);drawnow;
%       igure(2);quiver(vertex2(:,1),vertex2(:,2),Sx,Sy);drawnow;
%       figure(3);quiver(vertex2(:,1),vertex2(:,2),Idiff,Idiff);drawnow;
        Ux = -Idiff.*  ((Sx./((Sx.^2+Sy.^2)+alpha^2*Idiff.^2))+(Mx./((Mx.^2+My.^2)+alpha^2*Idiff.^2)));
        Uy = -Idiff.*  ((Sy./((Sx.^2+Sy.^2)+alpha^2*Idiff.^2))+(My./((Mx.^2+My.^2)+alpha^2*Idiff.^2)));
    
        % When divided by zero
        Ux(isnan(Ux))=0; Uy(isnan(Uy))=0;

        % Smooth the transformation field
        Uxs=imfilter(Ux,Hsmooth);
        Uys=imfilter(Uy,Hsmooth);

%         figure(15);quiver(x(1:16:end),y(1:16:end),Uxs(1:16:end),Uys(1:16:end));drawnow;
        % Add the new transformation field to the total transformation field.
        Tx=Tx+Uxs;
        Ty=Ty+Uys;
%         M=movepixels(M1,Ty,Tx);
        D(:,:,1) = Tx;
        D(:,:,2) = Ty;
        M = imwarp(initial_M,D);

        figure(20);
%             plot(V_m(:,1)+Tx,V_m(:,2)+Ty,'.');
        subplot(1,2,1);imshow(M,'InitialMagnification', 800); title('Source intensity')
        subplot(1,2,2);imshow(S,'InitialMagnification', 800); title('Target intensity')
        drawnow;        
end

end