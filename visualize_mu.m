f = [1,2, 3];
v = [-1, 0;
    1, 0;
    0, sqrt(3)];
[x, y] = meshgrid(linspace(-5,5,400), linspace(10, 0, 400));
I = zeros(400, 400);
Ia = zeros(400, 400);
for i = 1:400
    for j = 1:400
        vtmp = [v(1:2,:); x(i,j), y(i,j)];
        mu = compute_bc(f, v, vtmp, 2);
        I(i,j) = abs(mu);
        Ia(i,j) = angle(mu);
    end
end
% M = max(I(:)); m = 0;
% I2 = round(I / M * 64) + 1;
%%
figure(1);
imshow(I); colormap parula; colorbar;
%%
figure(2);
imshow(Ia); colormap hsv; colorbar;