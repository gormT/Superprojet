function [L] = LL(Nbtri, nbg, k, eq_radius, R, baseline_Gauss_points)

% This fonction is used to calculate the matrix M.

% save('args.mat', 'N', 'k', 'B', 'R');

% M = W.*exp(-1i*k.*R)./(4*pi.*R.^2).*(1./R + 1i*k);

L = zeros(Nbtri);
aux_M = zeros(Nbtri,1);
coef_Gauss = baseline_Gauss_points(:,3);

for t=1:Nbtri
    for q=1:nbg
        aux_M = aux_M + coef_Gauss(q)*exp(-1i*k.*R(:,t*nbg+q-nbg))./(4*pi.*R(:,t*nbg+q-nbg).^2).*(1./R(:,t*nbg+q-nbg) + 1i*k);
    end
    L(:,t) = aux_M;
    aux_M = zeros(Nbtri,1);
end