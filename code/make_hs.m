function coord_h = make_hs(coor_triangles,coor_orthocenters)
% Donne les longueurs h dans les triangles pour la méthode de Matsumoto
% hs = [h1 h2 h3]

% pente des côtés
mA = (coor_triangles(:,[6 2 4])-coor_triangles(:,[4 6 2]))./...
     (coor_triangles(:,[5 1 3])-coor_triangles(:,[3 5 1]));

% ordonnée à l'origine des côtés
bA = coor_triangles(:,[6 2 4])-mA.*coor_triangles(:,[5 1 3]);

% pente des hauteurs
mB = -1./mA;

% ordonnée à l'origine des hauteurs
bB = coor_orthocenters(:,2)-mB.*coor_orthocenters(:,1);

% abscisses des h
xh = (bB-bA)./(mA-mB);

% ordonnées des h
yh = mB.*xh+bB;

% coordonnées des points h
coord_h = [xh(:,1) yh(:,1) xh(:,2) yh(:,2) xh(:,3) yh(:,3)];

end