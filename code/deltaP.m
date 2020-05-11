%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                               DELTAP                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Setup
% This part of the code is here for clearing variables, selecting the
% desired case (embedded or general), specifying the output folder name,
% importing the mesh file, adding the necessary files and folders to the
% path, specifying the names of the parameter files to use and changing
% some parameters if needed.

clear variables;
clc;

condition = 'GENERAL';

MshFileName = 'cercle0057.msh';          % Name of the msh file
[~,FileName] = fileparts(MshFileName);

Project_path = 'C:\Users\gormt\Superprojet\Project'; % Path to the project folder
addpath(genpath(Project_path));         % Adds the project folder and its subfolders to the path

VKG_path = 'C:\Users\gormt\Superprojet\VK-Gong';     % Path to the VK-Gong folder
addpath(genpath(VKG_path));             % Adds the VK-Gong folder and its subfolders to the path


PlateCharacteristicsFileName = 'PlateCharacteristics.mat';  % Physical characteristics of the plate: Dimensions, imperfection profile, material and boundary conditions.  
SimulationParametersFileName = 'SimulationParameters.mat';  % Parameters related to the simulation: Time length, scheme, number of modes, output points, accuracy.
GammaFileName = 'GammaCircular-Nphi_NPHICircular.mat';      % Name of the file containing the Gamma Tensor.  
ScoreFileName = 'ScoreParameters.mat';                      % Characteristics of the excitation.
PressureParametersFileName = 'PressureParameters.mat';      % The position of the output point(out_points) and the size of the meshes (lc) in the mesh file and the number of points 
% Gauss in a maillage(nbg) and the time of the simulation(Tsd) and the density of the fluide(rhof).


fsd = 0;                    % sampling frequency, is updated later on
save(SimulationParametersFileName,'fsd','-append');
Tsd = 5;                    % time interval simulated
save(SimulationParametersFileName,'Tsd','-append');
Nphi = 100;                 % Nphi
save(SimulationParametersFileName,'Nphi','-append');
Npsi = 100;                 % Npsi
save(SimulationParametersFileName,'Npsi','-append');
nbg = 6;                    % number of Gauss points: choose amongst {3, 6, 7, 16, 19}
save(PressureParametersFileName,'nbg','-append');
sing = 0;                  % Singluarity : 2 for nothing, 0 for zero on diag, 1 for eq areas method, 3 for keeping only the singularity, 4 for eq areas method but only the singularity
save(SimulationParametersFileName,'sing','-append');
%out_points = [zeros(21,2) [-5:0.5:5]'];      % output points
out_points = [0 0 2];
save(PressureParametersFileName,'nbg','out_points','-append');

% Structure results folder name:
% Results_<case>_<output points format>_<msh file element size>_<number of Gauss points>_<singularity>
ResultsFolderName = sprintf('test_aire_0057_%d_%d',nbg,sing);      % Name of the results folder

out_nd_filename = strcat(FileName,sprintf('_Tsd=%d_out_nd.mat',Tsd));



load(PressureParametersFileName);
load(PlateCharacteristicsFileName, 'Rd');

%% Mean element size and collocation points
% Uses another function to determine the mean element size in the mesh. It
% also allows the user to plot the whole mesh, if needed.
% This function also computes the collocation points (centroids) of the
% triangles. This part of the calculation might change at some point to
% somwhere else.

% [lc,Nb_segments,segments,X1,X2,Y1,Y2,len_seg,coor_coll_points,coor_coll_points_pol,~] = avg_msh_size(MshFileName);
% save(PressureParametersFileName,'lc','-append');    % Redefining the element size int the mesh. This creates causality in the results, which is more accurate
% Normalization of the Ploar Coordinates for VK-gong
[~,Nbtri,Coorneu,~,Numtri] = lecture_msh(MshFileName);
[coor_coll_points,coor_coll_points_pol] = coord_coll(Nbtri,Numtri,Coorneu);
coor_coll_points_pol_norm(:,2) = coor_coll_points_pol(:,2)/Rd;
%% Gauss points
% This part of the code gives us the coordinates of the Gauss points and
% other useful data. See Gauss_points_coordinates for more info.

[baseline_Gauss_points, coor_Gauss_points, Nbtri, areas] = Gauss_points_coordinates(Nbtri,Coorneu,Numtri,nbg);
nbg_tot = Nbtri*nbg;    % Total number of Gauss points

coor_Gauss_points_pol = zeros(Nbtri*nbg,3);     % Initialization of the Polar Coordinates of the Gauss points in each mesh triangle

% Transform the Cartesian Coordinates of the Gauss points in the mesh triangles to the Polar Coordinates
for t =1:size(coor_Gauss_points,1)
    [coor_Gauss_points_pol(t,1),coor_Gauss_points_pol(t,2)] = cart2pol(coor_Gauss_points(t,1),coor_Gauss_points(t,2));
end

% Normalization of the polar coordinates for VK-gong
coor_Gauss_points_pol_norm = coor_Gauss_points_pol.*[1 1/Rd 1];     % Adimensions the radiuses for VK-gong

%% Building the ouput points on the plate

% THE CASE IS GENERAL
op = coor_coll_points_pol_norm(:,1:2);
save(SimulationParametersFileName,'op','-append');

%% Simulation setup
% This part sets up the necessary data for the simulation. It uses VK-Gong
% functions.

cd(VKG_path);

[Rd, hd, E, BC, e, Nphi, Npsi, scheme, H0, H1, H2, filename, Ai, C, C1, C2, k_t, c_t, xkn, JJ, II, Kkn, rp, tad, fs, Tsd] = plate_def_circ(PlateCharacteristicsFileName, SimulationParametersFileName, ResultsFolderName, GammaFileName);

[ f_time, Tn ] = score_circ(ScoreFileName, Rd, hd, E, BC, e, Nphi, scheme, C, k_t, c_t, xkn, JJ, II, Kkn, tad, fs, Tsd);

%% Computation of the movement of each Gauss point using VK-Gong

% switch scheme
%     case 'ECS'
%         [out_nd] = ftime_imperfect_ECS( Nphi, Npsi, Ai, H0, H1, H2, C, C1,C2, Tn, e, f_time, rp);
%         
%     case 'verlet'
%         [out_nd] = ftime_imperfect_verlet( Nphi, Npsi, Ai, H1, C, C1, C2, Tn, e, f_time, rp);
%     
%     otherwise
%         disp('Unknown scheme');
% end
cd(Project_path);


load('out_nd_folder/cercle0057_Tsd=5_out_nd.mat');

fsd = round(fs/tad);            % New sampling frequency
fprintf('The new sampling frequency is %d\n', fsd);
out = out_nd*hd;                % Dimensioned output displacement at the Gauss points.
out_vel = diff(out,1,1)*fsd;    % Dimensioned output velocity at the Gauss points.

cd(Project_path);
cd('DeltaP');

% carac_out_nd = sprintf('_Tsd=%d_nbg=%d_',Tsd,nbg);
% out_nd_filename = strcat(FileName,carac_out_nd,'out_nd.mat');
% % saving out_nd for later use
% fprintf('Saving out_nd...   ');
% save(out_nd_filename,'out_nd','-v7.3');
% fprintf('saved !\n');
%% Computation of the pressure

c = 340; % Speed of Sound


nb_out_points = size(out_points,1);   % nombre de points de sortie
eq_radius = sqrt(areas/pi);         % Equivalent radius of each mesh

%% Calculate the pression
nt = size(out_vel,1);

if mod(nt,2) == 1           % on s'assure que out_vel a une taille paire
    out_vel(nt+1,:) = 0;
end

nt = size(out_vel,1);

out_vel_fre = fft(out_vel);
% clear('out_vel');
% Calculate the distances of each Gauss point and each output point

DeltaP = zeros(nt/2+1,Nbtri);           % Il y a bien autant de colonnes que de points de collocation, OK !
dw = 2*pi/(nt/fsd);                     % Pas de fréquence
Dw = [0:dw:(nt/2-1)*dw,nt/2*dw];            % J'ai corrigé le signe du dernier terme
% W = repmat(areas',Nbtri,1);           % Cette matrice est a priori inutile du coup
R_coll_Gauss = new_distances(coor_coll_points,coor_Gauss_points);
% R_coll_Gauss = new_distances([5*ones(Nbtri,1) zeros(Nbtri,1)],coor_Gauss_points);   % A COMMENTER

Vn = out_vel_fre(1:nt/2+1,:);

% % vectorize the areas of meshes
% An = repmat(baseline_Gauss_points(:,3),Nbtri,1);
% Areas = repmat(areas,1,nbg);
% AN = reshape(Areas',Nbtri*nbg,1);
% AN = 2*AN.*An;





% Fréquences d'étude de DeltaP
f1 = 140;    % Pour f = 140Hz
f2 = 1400;   % Pour f = 1400Hz
% Pulsations
w1 = 2*pi*f1;
w2 = 2*pi*f2;
% Indices
N1 = round(w1/dw)+1;
N2 = round(w2/dw)+1;

SolFileName1 = strcat('DeltaP_',FileName,sprintf('_nbg=%d_%dHz_sing=%d.sol',nbg,f1,sing));
SolFileName2 = strcat('DeltaP_',FileName,sprintf('_nbg=%d_%dHz_sing=%d.sol',nbg,f2,sing));


SolFileID1 = fopen(SolFileName1,'w');
SolFileID2 = fopen(SolFileName2,'w');

fprintf(SolFileID1,'MeshVersionFormatted 3\n\nDimension 3\n\nSolAtTriangles\n%d\n1 1\n',Nbtri);
fprintf(SolFileID2,'MeshVersionFormatted 3\n\nDimension 3\n\nSolAtTriangles\n%d\n1 1\n',Nbtri);


fprintf('Calculating and writing in files the DeltaP...   ');
k1 = Dw(N1)/c;
[M1] = MM(Nbtri, nbg, k1, eq_radius, R_coll_Gauss, baseline_Gauss_points,areas,sing);
% DP1 = -1i*rhof*Dw(i1)*(M1\(Vn(i1,:).'));
% DP1 = M1*ones(Nbtri,1);
% dp1 = abs(DP1);
dp1 = M1*ones(Nbtri,1);
% DeltaP1(i1,:) = DP.';
for i=1:size(dp1)
    fprintf(SolFileID1,'%g\n',dp1(i));
end

fprintf('first one done...   ');

k2 = Dw(N2)/c;
[M2] = MM(Nbtri, nbg, k2, eq_radius, R_coll_Gauss, baseline_Gauss_points,areas,sing);
% DP2 = -1i*rhof*Dw(i2)*(M2\(Vn(i2,:).'));
% DP2 = M2*ones(Nbtri,1);
% dp2 = abs(DP2);
dp2 = M2*ones(Nbtri,1);
% DeltaP2(i2,:) = DP.';
for i=1:size(dp2)
    fprintf(SolFileID2,'%g\n',dp2(i));
end
fprintf('second one done !\n');


fprintf(SolFileID1,'\nEnd\n');
fprintf(SolFileID2,'\nEnd\n');

fclose('all');

M_out_filename = strcat('M_',FileName,sprintf('_f1=%d_f2=%d_nbg=%d_sing=%d.mat',f1,f2,nbg,sing));
% saving M1 and M2 for later use
fprintf('Saving M1 and M2...   ');
save(M_out_filename,'M1','M2','-v7.3');
fprintf('saved !\n');

fprintf('All done!\n');

