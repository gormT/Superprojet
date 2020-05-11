function [coor_coll_points,coor_coll_points_pol] = coord_coll(Nbtri,Numtri,Coorneu)

coor_coll_points = zeros(Nbtri,2);
coor_coll_points_pol = zeros(Nbtri,2);

for t = 1:Nbtri
    points_number = Numtri(t,:);
    point_sommet1_coor = Coorneu(points_number(1),:);
    point_sommet2_coor = Coorneu(points_number(2),:);
    point_sommet3_coor = Coorneu(points_number(3),:);
    coor_coll_points(t,:) = point_sommet1_coor/3 + point_sommet2_coor/3 + point_sommet3_coor/3;
end

% Transform the cartesian coordinates of the collocation points in the mesh triangles to the polar coordinates
for t =1:size(coor_coll_points,1)
    [coor_coll_points_pol(t,1),coor_coll_points_pol(t,2)] = cart2pol(coor_coll_points(t,1),coor_coll_points(t,2));
end
end