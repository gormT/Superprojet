function [baseline_Gauss_points, coor_Gauss_points, Nbtri, areas] = Gauss_points_coordinates(Nbtri,Coorneu,Numtri,nbg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              Project                              %
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function computes the coordinates of the Gauss points in a mesh, as
% well as the number of triangles in said mesh and their areas. The
% function also returns the coordinates of the Gauss points in the baseline
% triangle and their weight coefficients used in the discreet sum.

%% Input arguments:
%   MeshFileName: Name of the mesh file to use.
%   nbg: The number of Gauss points to use. Amongst {3, 6, 7, 16, 19}

%% Output arguments:
%   coor_Gauss_points: Coordinates of the Gauss points in the mesh.
%   baseline_Gauss_points: Coordinates and corresponding weights of the
% Gauss points in the baseline triangle ( (0,0) , (0,1) , (1,0) )
%   Nbtri: Number of triangles in the mesh.
%   areas: Areas of the triangles in the mesh.


Coorneu(:,3) = 0;                                           % Working in 3D, 3rd coordinate is zero

[baseline_Gauss_points] = coord_Gauss(nbg);           % Coordinates of the Gauss points in the baseline triangle
coor_Gauss_points = zeros(Nbtri*nbg,3);       % Initialization of the cartesian coordinates of the Gauss points in each one triangle
areas = zeros(Nbtri,1);                           % Initialization of the areas of each mesh triangle

% Computes the Cartesian Coordinates of the Gauss points in the mesh triangles
for t = 1:Nbtri
    points_number = Numtri(t,:);
    point_sommet1_coor = Coorneu(points_number(1),:);
    point_sommet2_coor = Coorneu(points_number(2),:);
    point_sommet3_coor = Coorneu(points_number(3),:);
    areas(t,1) = area(point_sommet1_coor,point_sommet2_coor,point_sommet3_coor);
    for q = 1:nbg
        coor_Gauss_points(t*nbg+q-nbg,:) = baseline_Gauss_points(q,1)*point_sommet1_coor...
            + baseline_Gauss_points(q,2)*point_sommet2_coor...
            + (1-baseline_Gauss_points(q,1)-baseline_Gauss_points(q,2))*point_sommet3_coor;
    end
end

end