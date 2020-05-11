function C = orthocenter(coor_triangles)
% Donne les coordonnées des orthocentres C pour une liste de triangles 
% [x1 y1 x2 y2 x3 y3]
% J'utilise les hauteurs issues des sommets 1 et 3, ça marche mais il
% faudrait que je remplace le sommet 3 par le sommet 2, pour améliorer la
% lisibilité de l'ensemble.
    
    % p =[p1 p3] les pentes des hauteurs issues des sommets 1 et 3
    p = [(coor_triangles(:,3)-(coor_triangles(:,5)))./(coor_triangles(:,6)-(coor_triangles(:,4)))...
         (coor_triangles(:,1)-(coor_triangles(:,3)))./(coor_triangles(:,4)-(coor_triangles(:,2)))];
    
    % q = [q1 q3] les ordonnées à l'origine des hauteurs issues de 1 et 3
    q = [coor_triangles(:,2)-p(:,1).*coor_triangles(:,1)...
         coor_triangles(:,6)-p(:,2).*coor_triangles(:,5)];
    
    % abscisses des orthocentres
    C = (q(:,1)-q(:,2))./(p(:,2)-p(:,1));
    
    % ordonnées des orthocentres
    C = [C p(:,1).*C+q(:,1)];

end