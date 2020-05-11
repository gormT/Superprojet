function R = distances(A,B)

% Input: A,B are matrix of n*3
% na,nb are the numbres of the points in A and B
% To calculate the distances
% Output: distances is a Matrix of nb*na. distances_ij is the distance
% between A_j and B_i.

% R(A,B) = R(A,B).'

na = size(A,1);
nb = size(B,1);
R = zeros(nb,na);

for i = 1:nb
    R(i,:) = sqrt((B(i,1)-A(:,1)).^2 + (B(i,2)-A(:,2)).^2 + (B(i,3)-A(:,3)).^2);
end