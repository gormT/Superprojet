function [P,dw,Dw] = pressure_general(baseline_Gauss_points, Nbtri, areas, nbg, out_points, out_vel, coor_Gauss_points, fsd, rhof, coor_coll_points, FileName, sing,R_theta,angle_steps,N_seg_int)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              Project                              %
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function computes the FFT of the pressure in the general case.

%% Input arguments:
%   baseline_Gauss_points: Coordinates and corresponding weights of the
% Gauss points in the baseline triangle ( (0,0) , (0,1) , (1,0) )
%   Nbtri: Number of triangles in the mesh.
%   areas: Areas of the triangles in the mesh.
%   nbg: The number of Gauss used. Amongst {3, 6, 7, 16, 19}
%   out_points: The coordinates of the output points.
%   out_vel: Dimensioned output velocity at the Gauss points.
%   coor_Gauss_points: Coordinates of the Gauss points in the mesh.
%   fsd: Sampling frequency used.
%   rhof: Density of the fluid.

%% Output arguments:
%   P: FFT of the pressure at the output points.

%% Computation

c = 340; % Speed of Sound



% cd(path);
% [F, out_vel, fsd, areas, points_Gauss_coor_local, points_Gauss_coor_norm, Nbtri] = Found_Displace_circ(OutputFileName,lc,nbg,Tsd);


nb_out_points = size(out_points,1);   % nombre de points de sortie
eq_radius = sqrt(areas/pi);         % Equivalent radius of each mesh

%% Calculate the pression
nt = size(out_vel,1);

if mod(nt,2) == 1           % on s'assure que out_vel a une taille paire
    out_vel(nt+1,:) = 0;
end

nt = size(out_vel,1);

out_vel_fre = fft(out_vel);
clear('out_vel');
% Calculate the distances of each Gauss point and each output point

DeltaP = zeros(nt/2+1,Nbtri);           % Il y a bien autant de colonnes que de points de collocation, OK !
dw = 2*pi/(nt/fsd);                     % Pas de fréquence
Dw = [0:dw:(nt/2-1)*dw,nt/2*dw];        % J'ai corrigé le signe du dernier terme
% W = repmat(areas',Nbtri,1);           % Cette matrice est a priori inutile du coup
R_coll_Gauss = new_distances(coor_coll_points,coor_Gauss_points);

Vn = out_vel_fre(1:nt/2+1,:);

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

cd('DeltaP')

SolFileID1 = fopen(SolFileName1,'w');
SolFileID2 = fopen(SolFileName2,'w');

fprintf(SolFileID1,'MeshVersionFormatted 3\n\nDimension 3\n\nSolAtTriangles\n%d\n1 1\n',Nbtri);
fprintf(SolFileID2,'MeshVersionFormatted 3\n\nDimension 3\n\nSolAtTriangles\n%d\n1 1\n',Nbtri);


fprintf('\n--------------Begin to calculate the DeltaP-----------------\n');

% WaitMessage = parfor_wait(nt/2+1, 'Waitbar', true);

tLoop = tic;

parfor i = 1:nt/2+1
    k = Dw(i)/c;
    [M] = make_M(Nbtri, nbg, k, R_coll_Gauss, baseline_Gauss_points, areas, sing, R_theta, angle_steps, N_seg_int);
    DP = -1i*rhof*Dw(i)*(M\(Vn(i,:).'));
    DeltaP(i,:) = DP.';
    
    if (i == round(nt/40))
        Time5 = toc(tLoop);
        EstimatedTime = Time5*20;
        hours = floor(EstimatedTime / 3600);
        EstimatedTime = EstimatedTime - hours * 3600;
        mins = floor(EstimatedTime / 60);
        secs = EstimatedTime - mins * 60;
        disp(['Estimated time for the time integration: ', num2str(hours), ' hours ', num2str(mins), ' minutes ', num2str(secs), ' seconds ']);
    end
    
    if mod(i,round(nt/40)) == 0
        disp(['Progress: ', num2str(round(200*i/nt)), '%']);
    end
%     WaitMessage.Send;
%     pause(0.002);
    
end

% WaitMessage.Destroy

dp1 = abs(DeltaP(N1,:)');
dp2 = abs(DeltaP(N2,:)');

fprintf('Writing first DeltaP sol file (f = %dHz)...     ',f1);
for j=1:size(dp1)
    fprintf(SolFileID1,'%g\n',dp1(j));
end
fprintf('done.\n');

fprintf('Writing second DeltaP sol file (f = %dHz)...     ',f1);
for j=1:size(dp1)
    fprintf(SolFileID2,'%g\n',dp2(j));
end
fprintf('done.\n');

fprintf(SolFileID1,'\nEnd\n');
fprintf(SolFileID2,'\nEnd\n');

fclose('all');

cd ..

% mkdir('DeltaP');
% cd('DeltaP');
% save(sprintf('DeltaP-Tsd_%f-lc_%f.mat',Tsd,lc), 'DeltaP', 'fsd', 'points_Gauss_coor_local','areas','Nbtri','-v7.3');
% cd ..;

% Clearing another possibly large array
clear out_vel_fre

fprintf('\n------------The DeltaP calculation is done----------------\n');


P = zeros(nt,nb_out_points); % Initialization of the pression

% Gauss Quadrature
fprintf('\n--------------Beginning to calculate the pressure-----------------\n');

R_out_coll = new_distances(out_points,coor_coll_points);
n = out_points(:,3);
n = repmat(n,1,Nbtri);
% W = repmat(areas',nb_out_points,1);       % Idem, osef

for l = 1:nt/2+1
    k = Dw(l)/c; % wave number
    L = -exp(-1i*k.*R_out_coll)./(4*pi*R_out_coll).*(1./R_out_coll + 1i*k).*(n./R_out_coll);
    PP = L*DeltaP(l,:).';
    P(l,:) = PP.';
    if mod(l,round(nt/40)) == 0
        disp(['Progress: ', num2str(round(200*l/nt)), '%']);
    end
end

% Because P is Conjugate symmetry, we can get the rest half part of P simply
P(nt/2+2:end,:) = conj(flipud(P(2:nt/2,:))); 

fprintf('\n------------The pressure calculation is done----------------\n');


end

