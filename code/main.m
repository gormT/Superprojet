%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                               MAIN                                % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is the main file used to compute the acoustic pressure in both the
% embedded and general cases.
% Simulating the vibration of a circular plate

%  !! Run this file while being in the 'Project' folder !!



%% Setup
% This part of the code is here for clearing variables, selecting the
% desired case (embedded or general), specifying the output folder name,
% importing the mesh file, adding the necessary files and folders to the
% path, specifying the names of the parameter files to use and changing
% some parameters if needed.

clear variables;
clc;

% condition = 'EMBEDDED';            % Select the case : EMBEDDED or GENERAL
condition = 'GENERAL';

MshFileName = 'cercle0037.msh';          % Name of the msh file
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
nbg = 7;                    % number of Gauss points: choose amongst {3, 6, 7, 16, 19}
save(PressureParametersFileName,'nbg','-append');
sing = 1;                  % Singluarity : 2 for nothing, 0 for zero on diag, 1 for eq areas method, 3 for keeping only the singularity
save(SimulationParametersFileName,'sing','-append');
%out_points = [zeros(21,2) [-5:0.5:5]'];      % output points
out_points = [0 0 2];
save(PressureParametersFileName,'nbg','out_points','-append');

% Structure results folder name:
% Results_<case>_<output points format>_<msh file element size>_<number of Gauss points>_<singularity>
ResultsFolderName = sprintf('Matsumoto_0037_%d_%d',nbg,sing);      % Name of the results folder

out_nd_filename = strcat(FileName,sprintf('_Tsd=%d_out_nd.mat',Tsd));



load(PressureParametersFileName);
load(PlateCharacteristicsFileName, 'Rd');

%% Collocation points

[~,Nbtri,Coorneu,~,Numtri] = lecture_msh(MshFileName);      % Reading the msh file
[coor_coll_points,coor_coll_points_pol] = coord_coll(Nbtri,Numtri,Coorneu);
coor_coll_points_pol_norm(:,2) = coor_coll_points_pol(:,2)/Rd;          % Normalization of the Ploar Coordinates for VK-gong
%% Gauss points

[baseline_Gauss_points, coor_Gauss_points, Nbtri, areas] = Gauss_points_coordinates(Nbtri,Coorneu,Numtri,nbg);
nbg_tot = Nbtri*nbg;    % Total number of Gauss points

coor_Gauss_points_pol = zeros(Nbtri*nbg,3);     % Initialization of the Polar Coordinates of the Gauss points in each mesh triangle

% Transform the Cartesian Coordinates of the Gauss points in the mesh triangles to the Polar Coordinates
for t =1:size(coor_Gauss_points,1)
    [coor_Gauss_points_pol(t,1),coor_Gauss_points_pol(t,2)] = cart2pol(coor_Gauss_points(t,1),coor_Gauss_points(t,2));
end

% Normalization of the polar coordinates for VK-gong
coor_Gauss_points_pol_norm = coor_Gauss_points_pol.*[1 1/Rd 1];     % Adimensions the radiuses for VK-gong

%% Orthocenters, their angles and the whole shebang

% coordinates of the triangles [x1 y1 x2 y2 x3 y3]
coor_triangles = [Coorneu(Numtri(:,1),1:2) Coorneu(Numtri(:,2),1:2) Coorneu(Numtri(:,3),1:2)];

% coordinates of the orthocenters in each triangles
% coor_orthocenters = orthocenter(coor_triangles);

% R_theta_ref in the triangles for the Matsumoto method [R1 R2 R3]
R_theta_ref = [sqrt((coor_triangles(:,1)-coor_coll_points(:,1)).^2+(coor_triangles(:,2)-coor_coll_points(:,2)).^2)...
               sqrt((coor_triangles(:,3)-coor_coll_points(:,1)).^2+(coor_triangles(:,4)-coor_coll_points(:,2)).^2)...
               sqrt((coor_triangles(:,5)-coor_coll_points(:,1)).^2+(coor_triangles(:,6)-coor_coll_points(:,2)).^2)];

% coordinates of the h points in the triangles for the Matsumoto method [xh1 yh1 xh2 yh2 xh3 yh3]
coord_h = make_hs(coor_triangles,coor_coll_points);

% lengths of the h points
lengths_h = [sqrt((coord_h(:,1)-coor_coll_points(:,1)).^2+(coord_h(:,2)-coor_coll_points(:,2)).^2)...
               sqrt((coord_h(:,3)-coor_coll_points(:,1)).^2+(coord_h(:,4)-coor_coll_points(:,2)).^2)...
               sqrt((coord_h(:,5)-coor_coll_points(:,1)).^2+(coord_h(:,6)-coor_coll_points(:,2)).^2)];

% angles in the triangles for the Matsumoto method [1_1 1_2 2_1 2_2 3_1 3_2]
angles_triangles = [-acos(lengths_h(:,1)./R_theta_ref(:,2))...
                      acos(lengths_h(:,1)./R_theta_ref(:,3))...
                     -acos(lengths_h(:,2)./R_theta_ref(:,3))...
                      acos(lengths_h(:,2)./R_theta_ref(:,1))...
                     -acos(lengths_h(:,3)./R_theta_ref(:,1))...
                      acos(lengths_h(:,3)./R_theta_ref(:,2))];

% number of segments per sector for the numerical integration
N_seg_int = 10;

% angle steps for the numerical integration
angle_steps = (angles_triangles(:,[2 4 6])- angles_triangles(:,[1 3 5]))/N_seg_int;

% the angles for the numerical integration
angles = make_angles(angles_triangles,angle_steps,Nbtri,N_seg_int);

% the R_theta for the numerical integration
R_theta = [repmat(lengths_h(:,1),1,N_seg_int+1)...
            repmat(lengths_h(:,2),1,N_seg_int+1)...
            repmat(lengths_h(:,2),1,N_seg_int+1)];
R_theta = R_theta./angles;



%% Building the ouput points on the plate

% A switch depending on the selected condition for op, the output points on
% the plate change in each condition.
switch condition
    case 'EMBEDDED'     % Rayleigh, embedded case
        op = coor_Gauss_points_pol_norm(:,1:2);
        save(SimulationParametersFileName,'op','-append');
    case 'GENERAL'      % General case
        op = coor_coll_points_pol_norm(:,1:2);
        save(SimulationParametersFileName,'op','-append');
end

%% Simulation setup
% This part sets up the necessary data for the simulation. It uses VK-Gong
% functions.

cd(VKG_path)

[Rd, hd, E, BC, e, Nphi, Npsi, scheme, H0, H1, H2, filename, Ai, C, C1, C2, k_t, c_t, xkn, JJ, II, Kkn, rp, tad, fs, Tsd] = plate_def_circ(PlateCharacteristicsFileName, SimulationParametersFileName, ResultsFolderName, GammaFileName);

[ f_time, Tn ] = score_circ(ScoreFileName, Rd, hd, E, BC, e, Nphi, scheme, C, k_t, c_t, xkn, JJ, II, Kkn, tad, fs, Tsd);

%% Computation of the movement of each Gauss point using VK-Gong

cd(Project_path)
if isfile(strcat('./out_nd_folder/',out_nd_filename))
    fprintf('out_nd file found, loading...     ');
    load(strcat('./out_nd_folder/',out_nd_filename));
    fprintf('done.\n');
else
    cd(VKG_path)
    fprintf('No out_nd file found, calculating...     ');
    switch scheme
        case 'ECS'
            [out_nd] = ftime_imperfect_ECS( Nphi, Npsi, Ai, H0, H1, H2, C, C1,C2, Tn, e, f_time, rp);
            
        case 'verlet'
            [out_nd] = ftime_imperfect_verlet( Nphi, Npsi, Ai, H1, C, C1, C2, Tn, e, f_time, rp);
            
        otherwise
            disp('Unknown scheme');
    end
    fprintf('done.\n');
    cd(Project_path)
    cd('out_nd_folder')
    carac_out_nd = sprintf('_Tsd=%d_nbg=%d_',Tsd,nbg);
    out_nd_filename = strcat(FileName,carac_out_nd,'out_nd.mat');
    % saving out_nd for later use
    fprintf('Saving out_nd...   ');
    save(out_nd_filename,'out_nd','-v7.3');
    fprintf('saved !\n');
end





fsd = round(fs/tad);            % New sampling frequency
fprintf('The new sampling frequency is %d\n', fsd);
out = out_nd*hd;                % Dimensioned output displacement at the Gauss points.
out_vel = diff(out,1,1)*fsd;    % Dimensioned output velocity at the Gauss points.

% Clearing some possibly large arrays
clear out_nd out

%% Computation of the pressure

switch condition
    
    case 'EMBEDDED'      % Rayleigh, embedded case
        
        disp('Case chosen : Rayleigh');
        P = pressure_Rayleigh(baseline_Gauss_points, Nbtri, areas, nbg, out_points, out_vel, coor_Gauss_points, fsd, rhof);

        
    case 'GENERAL'      % General case
        
        disp('Case chosen : general case (3D)');
        
        % We redefine the output points in the 3D case
%         x = -0.3:0.01:0.3;
%         z = -0.2:0.01:0.2;
%         [X,Z] = meshgrid(x,z);
%         out_points_x = reshape(X,[],1);
%         out_points_z = reshape(Z,[],1);
%         out_points = zeros(size(out_points_x,1),3);
%         out_points(:,1) = out_points_x;
%         out_points(:,3) = out_points_z;
%         out_points(:,2) = 0;
        
%         save('PressureParameters-CaseP3.mat','out_points','-append');
        
        
        [P,dw,Dw] = pressure_general(baseline_Gauss_points, Nbtri, areas,...
                    nbg, out_points, out_vel, coor_Gauss_points, fsd, rhof,...
                    coor_coll_points, FileName, sing,R_theta,angle_steps,N_seg_int);
        
end

% Clearing another possibly large array
clear out_vel

%% Computation of the IFFT of P

P_t = ifft(P,'symmetric');


%% Saving and displaying the results
% The user can save as many of the results as needed.

cd(Project_path);
cd('results/');

% The following lines of code create a folder with an added suffix in the
% case of the folder already existing. The user is then advised to change
% the results folder name.

Nb_folder = 0;        % The copy number for the results folder
NewFolderName = ResultsFolderName;
while isfolder(NewFolderName)     % Checking whether or not the folder for the results exists
    Nb_folder = Nb_folder+1;            % We keep increasing the Nb_folder number until there is no folder conflict
    NewFolderName = strcat(ResultsFolderName,'_number_',num2str(Nb_folder));
end
real_Nb_folder = Nb_folder-1;
if Nb_folder == 1
    fprintf(2,'The folder ''%s'' already exists in ''results/''.\nPlease consider changing folder name.\n', ResultsFolderName);
elseif Nb_folder > 1
    fprintf(2,'The folder ''%s'' already exists in ''results/'', together with %d other folders with the same beginning.\nPlease consider changing folder name.\n', ResultsFolderName, real_Nb_folder);
end
mkdir(NewFolderName);
cd(NewFolderName);



% Saving the results, the user can add data to save
save(sprintf('Results-Tsd_%d-nbg_%d.mat',Tsd,nbg),'P','P_t','fsd','nbg_tot','-v7.3');


