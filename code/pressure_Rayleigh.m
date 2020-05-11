function [P] = pressure_Rayleigh(baseline_Gauss_points, Nbtri, areas, nbg, out_points, out_vel, coor_Gauss_points, fsd, rhof)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              Project                              %
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function computes the FFT of the pressure in the embedded case,
% using Rayleigh's integral.

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


% vectorize the areas of meshes
An = repmat(baseline_Gauss_points(:,3),Nbtri,1);
Areas = repmat(areas,1,nbg);
AN = reshape(Areas',Nbtri*nbg,1);
AN = 2*AN.*An;

N = Nbtri*nbg; % the total number of point Gauss
rop = size(out_points,1); % the number of the output points

nt = size(out_vel,1); % the total number of the sampling of the velocity

% If the number of sampling is not even, we add a 0 at the end to make it
% even. This step is just for simplifing the code.
if mod(nt,2) == 1
    out_vel(nt+1,:) = 0;
    nt = nt + 1;
end

%% Take fft to the velocity
out_vel_fre = fft(out_vel);

% Clear some large data not needed any more
clear('out_vel','C','C1','C2','H0','H1','H2');

% Calculate the distances of each Gauss point and each output point
R = distances(out_points,coor_Gauss_points);

dw = 2*pi/(nt/fsd); % Frequency interval of the velocity in friquency domain.
c = 340; % Speed of Sound
P = zeros(nt,rop); % Initialization of the pression

fprintf('\n--------------Beginning to calculate the pressure-----------------\n');

Dw = [0:dw:(nt/2-1)*dw,-nt/2*dw];
DW = repmat(Dw',1,rop); % Vectorize the frequency according to each element of the velocity in friquency domain

tLoop = tic;

for t = 1:N
    
    RR = repmat(R(t,:),nt/2+1,1); % Vectorize the distances
    Vn = repmat(out_vel_fre(1:(nt/2+1),t),1,rop);
    
%% Pressure in frequency domain
    P(1:(nt/2+1),:) = 1i*rhof*DW/(2*pi).*Vn.*AN(t)./RR.*exp(-1i*DW./c.*RR) + P(1:(nt/2+1),:);
    
    if (t == round(N/20))
        Time5 = toc(tLoop);
        EstimatedTime = Time5*20;
        hours = floor(EstimatedTime / 3600);
        EstimatedTime = EstimatedTime - hours * 3600;
        mins = floor(EstimatedTime / 60);
        secs = EstimatedTime - mins * 60;
        disp(['Estimated time for the time integration: ', num2str(hours), ' hours ', num2str(mins), ' minutes ', num2str(secs), ' seconds ']);
    end
    
    if mod(t,round(N/20)) == 0
        disp(['Progress: ', num2str(round(100*t/(N))), '%']);
    end
end

P(nt/2+2:end,:) = conj(flipud(P(2:nt/2,:)));

fprintf('\n------------Finishing the pressure calculation----------------\n');

end