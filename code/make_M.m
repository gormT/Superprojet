function M = make_M(Nbtri, nbg, k, R, baseline_Gauss_points, areas, sing, R_theta, angle_steps, N_seg_int)

M = zeros(Nbtri);
coef_Gauss = baseline_Gauss_points(:,3);


switch sing
    case -1         % case in which the function used it the matrix M is one, for test purposes
        aux_M = zeros(Nbtri,1);
        for t=1:Nbtri
            for q=1:nbg
                aux_M = aux_M + coef_Gauss(q);
            end
            M(:,t) = aux_M*areas(t)*2;
            aux_M = zeros(Nbtri,1);
        end
    case 0          % default case : don't change anything
        aux_M = zeros(Nbtri,1);
        for t=1:Nbtri
            for q=1:nbg
                aux_M = aux_M + coef_Gauss(q)*exp(1i*k.*R(:,t*nbg+q-nbg))./(4*pi.*R(:,t*nbg+q-nbg).^3).*(1 - 1i*k*R(:,t*nbg+q-nbg));
            end
            M(:,t) = aux_M*areas(t)*2;
            aux_M = zeros(Nbtri,1);
        end
    case 1          % case in which the singularity is taken care of : Matsumoto method
        aux_M = zeros(Nbtri,1);
        for t=1:Nbtri
            for q=1:nbg
%                 current_R = [R([1:t-1],(t-1)*nbg+q); 5 ; R([t+1:Nbtri],(t-1)*nbg+q)];
%                 aux_M = aux_M + coef_Gauss(q)*exp(1i*k.*current_R)./(4*pi.*current_R.^3).*(1 - 1i*k*current_R);
                  aux_M = aux_M + coef_Gauss(q)*exp(1i*k.*R(:,t*nbg+q-nbg))./(4*pi.*R(:,t*nbg+q-nbg).^3).*(1 - 1i*k*R(:,t*nbg+q-nbg));
            end
            M(:,t) = aux_M*areas(t)*2;
            aux_M = zeros(Nbtri,1);
        end
        M(logical(eye(Nbtri))) = make_diag_M(R_theta,angle_steps,k,N_seg_int,Nbtri);
end







end