function [M] = MM(Nbtri, nbg, k, eq_radius, R, baseline_Gauss_points, areas, sing)

% This fonction is used to calculate the matrix M.

% save('args.mat', 'N', 'k', 'B', 'R');

% M = W.*exp(-1i*k.*R)./(4*pi.*R.^2).*(1./R + 1i*k);

M = zeros(Nbtri);

coef_Gauss = baseline_Gauss_points(:,3);

% for t=1:Nbtri
%     for q=1:nbg
%         aux_M = aux_M + coef_Gauss(q)*exp(-1i*k.*R(:,t*nbg+q-nbg))./(4*pi.*R(:,t*nbg+q-nbg).^2).*(1./R(:,t*nbg+q-nbg) + 1i*k);
% %         aux_M = aux_M + coef_Gauss(q);
%     end
%     M(:,t) = aux_M*areas(t)*2;
%     aux_M = zeros(Nbtri,1);
% end

% Traitement de la singularité
switch sing
    case 0          % 0 on the diagonal
        aux_M = zeros(Nbtri,1);
        for t=1:Nbtri
            for q=1:nbg
%                 aux_M = aux_M + coef_Gauss(q)*exp(-1i*k.*R(:,t*nbg+q-nbg))./(4*pi.*R(:,t*nbg+q-nbg).^2).*(1./R(:,t*nbg+q-nbg) + 1i*k);
                aux_M = aux_M + coef_Gauss(q)*exp(1i*k.*R(:,t*nbg+q-nbg))./(4*pi.*R(:,t*nbg+q-nbg).^3).*(1 - 1i*k*R(:,t*nbg+q-nbg));
                %         aux_M = aux_M + coef_Gauss(q);
                if t==5
                    coef_Gauss(q)
                    R(1,t*nbg+q-nbg)
                end
            end
            M(:,t) = aux_M*areas(t)*2;
            if t==5
                M(1,5)
            end
            aux_M = zeros(Nbtri,1);
        end
%         M(logical(eye(Nbtri))) = 0; 
    case 1          % circular equivalent patches
        aux_M = zeros(Nbtri,1);
        for t=1:Nbtri
            for q=1:nbg
%                 aux_M = aux_M + coef_Gauss(q)*exp(-1i*k.*R(:,t*nbg+q-nbg))./(4*pi.*R(:,t*nbg+q-nbg).^2).*(1./R(:,t*nbg+q-nbg) + 1i*k);
                aux_M = aux_M + coef_Gauss(q)*exp(1i*k.*R(:,t*nbg+q-nbg))./(4*pi.*R(:,t*nbg+q-nbg).^3).*(1 - 1i*k*R(:,t*nbg+q-nbg));
                %         aux_M = aux_M + coef_Gauss(q);
            end
            M(:,t) = aux_M*areas(t)*2;
            aux_M = zeros(Nbtri,1);
        end
%         M(logical(eye(Nbtri))) = -exp(-1i*k*eq_radius)./(2*eq_radius) - 1i*k/2; % Assign values to diagonal elements
        M(logical(eye(Nbtri))) = -exp(-1i*k*eq_radius)./(2*eq_radius) + 1i*k/2; % Assign values to diagonal elements
    case 2          % Nothing : keeping it as is
        aux_M = zeros(Nbtri,1);
        for t=1:Nbtri
            for q=1:nbg
                aux_M = aux_M + coef_Gauss(q)*exp(-1i*k.*R(:,t*nbg+q-nbg))./(4*pi.*R(:,t*nbg+q-nbg).^2).*(1./R(:,t*nbg+q-nbg) + 1i*k);
                %         aux_M = aux_M + coef_Gauss(q);
            end
            M(:,t) = aux_M*areas(t)*2;
            aux_M = zeros(Nbtri,1);
        end
    case 3          % Keeping only the singularity (diagonal)
        aux = 0;
        for t=1:Nbtri
            for q=1:nbg
                aux = aux + coef_Gauss(q)*exp(-1i*k.*R(t,t*nbg+q-nbg))./(4*pi.*R(t,t*nbg+q-nbg).^2).*(1./R(t,t*nbg+q-nbg) + 1i*k);
            end
            M(t,t) = aux*areas(t)*2;
            aux = 0;
        end
    case 4          % Eq areas method but keeping only the singularity
        M(logical(eye(Nbtri))) = -exp(-1i*k*eq_radius)./(2*eq_radius) - 1i*k/2; % Assign values to diagonal elements
    case 5
        aux_M = zeros(Nbtri,1);
        for t=1:Nbtri
            for q=1:nbg
%                 aux_M = aux_M + coef_Gauss(q)*exp(-1i*k.*R(:,t*nbg+q-nbg))./(4*pi.*R(:,t*nbg+q-nbg).^2).*(1./R(:,t*nbg+q-nbg) + 1i*k);
                aux_M = aux_M + coef_Gauss(q);
            end
            M(:,t) = aux_M*areas(t)*2;
            aux_M = zeros(Nbtri,1);
        end
    case 6
        aux_M = zeros(Nbtri,1);
        for t=1:Nbtri
            for q=1:nbg
%                 aux_M = aux_M + coef_Gauss(q)*exp(-1i*k.*R(:,t*nbg+q-nbg))./(4*pi.*R(:,t*nbg+q-nbg).^2).*(1./R(:,t*nbg+q-nbg) + 1i*k);
                aux_M = aux_M + coef_Gauss(q);
            end
            M(:,t) = aux_M*areas(t)*2;
            aux_M = zeros(Nbtri,1);
        end
end

end

