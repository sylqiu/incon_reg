function index = get_index_by_location(point, vertex)
d = dist(vertex, point);
[~,index] = min(d);

end

function d = dist(P, q)
if size(P, 2) == 3
    Pq = [P(:,1)-q(1),P(:,2)-q(2),P(:,3)-q(3)];
elseif size(P, 2) == 2
    Pq = [P(:,1)-q(1),P(:,2)-q(2)];
end
d = dot(Pq,Pq,2);
end