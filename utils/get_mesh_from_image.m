function [face, vertex, intensity] = get_mesh_from_image(image_path, num_p_x)
if isstring(image_path)
    I = im2double(imread(image_path));
else
    I = image_path;
if length(size(I)) == 3
   I = mean(I, 3); 
end
[H, W] = size(I, 1, 2);
num_p_y = round(num_p_x * (H-1) / (W-1))+1;
[y_ind, x_ind] = ndgrid(linspace(1, H, num_p_y),...
                linspace(1, W, num_p_x));
y = (H-1) / (W-1) - (y_ind-1) / (W-1);
x = (x_ind-1) / (W-1);
vertex = [x(:), y(:)];
face = delaunay(vertex);
% bilinear sample the image to the vertex of the mesh
[y_im, x_im] = ndgrid(1:H,1:W);
intensity = zeros(num_p_y * num_p_x, size(I,3));
for ch = 1:size(I,3)
    I_ch = I(:,:,ch);
    intensity(:,ch) = griddata(x_im(:), y_im(:), I_ch(:), x_ind(:), y_ind(:));

end
end