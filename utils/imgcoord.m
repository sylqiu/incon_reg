function q = imgcoord(I,p)
% convert from pixel index to coodinate, assume left down corner is orgin
q = p;
q(:,2) = size(I,1)+1 - q(:,2);
end