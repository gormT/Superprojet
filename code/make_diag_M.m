function diag_M = make_diag_M(R_theta,angle_steps,k,N,Nbtri)
% Donne la diagonale de M en suivant la méthode de Matsumoto pour le
% traitement de la singularité.

%% Initialisation de la diagonale et du vecteur auxiliaire pour les calculs
diag_M = zeros(Nbtri,1);
aux = zeros(Nbtri,1);

%% Intégration 1er secteur : 2 à N
for p=1:N-1
    aux = aux+exp(1i*k*R_theta(:,p+1))./(4*pi*R_theta(:,p+1));
end

% Intégration 1er secteur : 1 et N+1
aux = aux+0.5*(exp(1i*k*R_theta(:,1))./(4*pi*R_theta(:,1))+...
                    exp(1i*k*R_theta(:,N+1))./(4*pi*R_theta(:,N+1)));

% Multiplication par le pas angulaire d'intégration 1er secteur
aux = aux.*angle_steps(:,1);

% Ajout sur la diagonale de l'intégration du 1er secteur
diag_M = diag_M + aux;

% Remise à zéro du vecteur auxiliaire
aux = zeros(Nbtri,1);




%% Intégration 2eme secteur : N+3 à 2N+1
for p=1:N-1
    aux = aux+exp(1i*k*R_theta(:,p+N+2))./(4*pi*R_theta(:,p+N+2));
end

% Intégration 2eme secteur : N+2 et 2N+2
aux = aux+0.5*(exp(1i*k*R_theta(:,N+2))./(4*pi*R_theta(:,N+2))+...
                    exp(1i*k*R_theta(:,2*N+2))./(4*pi*R_theta(:,2*N+2)));

% Multiplication par le pas angulaire d'intégration 2eme secteur
aux = aux.*angle_steps(:,2);

% Ajout sur la diagonale de l'intégration du 2eme secteur
diag_M = diag_M + aux;

% Remise à zéro du vecteur auxiliaire
aux = zeros(Nbtri,1);


%% Intégration 3eme secteur : 2N+4 à 3N+2
for p=1:N-1
    aux = aux+exp(1i*k*R_theta(:,p+2*N+3))./(4*pi*R_theta(:,p+2*N+3));
end

% Intégration 3eme secteur : 2N+3 et 3N+3
aux = aux+0.5*(exp(1i*k*R_theta(:,2*N+3))./(4*pi*R_theta(:,2*N+3))+...
                    exp(1i*k*R_theta(:,3*N+3))./(4*pi*R_theta(:,3*N+3)));

% Multiplication par le pas angulaire d'intégration 3eme secteur
aux = aux.*angle_steps(:,3);

% Ajout sur la diagonale de l'intégration du 2eme secteur
diag_M = diag_M + aux;


% Structure de l'intégration : e pour extrémités et i pour intermédiaire
% [1 2 ... N N+1 N+2 N+3 ... 2N+1 2N+2 2N+3 2N+4 ... 3N+2 3N+3]
%  e iiiiiii  e   e   iiiiiiiiii    e    e    iiiiiiiiii    e




end