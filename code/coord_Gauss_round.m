function [points_Gauss_coor_norm] = coord_Gauss_round(nb_points)

% This function can output the Gauss points in the normal
% triangle(x>0, y>0, 1-x-y>0) and corresponding weights.
% points_Gauss_coor_norm(:,1) is coordinate x,
% points_Gauss_coor_norm(:,2) is coordinate y,
% points_Gauss_coor_norm(:,3) is the corresponding weight.
switch nb_points
    case 3
        points_Gauss_coor_norm = [ 0.17 0.17 0.1667;
                                   0.17 0.66 0.1667;
                                   0.66 0.17 0.1667];
    case 6
        points_Gauss_coor_norm = [ 0.45 0.45 0.1117;
                                   0.45 0.10 0.1117;
                                   0.10 0.45 0.1117
                                   0.09 0.09 0.0550;
                                   0.82 0.09 0.0550;
                                   0.09 0.82 0.0550];
    case 7
        points_Gauss_coor_norm = [0.33 0.33 0.1125;
                                  0.47 0.47 0.0662;
                                  0.47 0.06 0.0662;
                                  0.06 0.47 0.0662;
                                  0.10 0.10 0.0630;
                                  0.10 0.80 0.0630;
                                  0.80 0.10 0.0630];
    case 16
        points_Gauss_coor_norm = [0.33 0.33 0.0722;
                                  0.08 0.46 0.0475;
                                  0.46 0.08 0.0475;
                                  0.46 0.08 0.0475;
                                  0.90 0.05 0.0162;
                                  0.05 0.05 0.0162;
                                  0.05 0.90 0.0162;
                                  0.66 0.17 0.0516;
                                  0.17 0.66 0.0516;
                                  0.17 0.17 0.0516;
                                  0.008 0.728 0.0136;
                                  0.728 0.008 0.0136;
                                  0.264 0.008 0.0136;
                                  0.008 0.264 0.0136;
                                  0.728 0.264 0.0136;
                                  0.264 0.728 0.0136];
end