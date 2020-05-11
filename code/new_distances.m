function R = new_distances(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% new_distances.m:
% Computes the euclidean distances between points in two arrays. The rows
% are the coordinates of one point. The number of rows determine the number
% of points, and the number of columns determine the dimension of the
% points. If A and B contain points of different dimensions, the function
% fills the smallest dimension matrix with zeros to equal the dimension of
% the biggest dimension matrix.
% If given only one matrix, computes the distances between the 
% points of said matrix, thus making R a symmetric matrix, with a diagonal
% full of zeros.
%
%
% INPUT  - A,B : arrays containing the cartesian coordinates of the points
%
% OUTPUT - R : array containing the euclidean distances beteween the points
%        in each arrays (A and B). R(i,j) = distance between point i of A
%        and point j of B. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 2
        [mA, nA] = size(A);
        [mB, nB] = size(B);
        
        size_diff = nA-nB;
        
        if size_diff<0                      % Equalizing the dimensions
            A = [A zeros(mA,-size_diff)];
        else
            B = [B zeros(mB,size_diff)];
        end
        
        R = zeros(mA,mB);
        for t=1:mA
            R(t,:) = sqrt(sum((A(t,:)-B).^2,2))';   % Euclidean distance
        end
    case 1
        mA = size(A,1);
        R = zeros(mA);
        for t=1:mA
            R(t,:) = sqrt(sum((A(t,:)-A).^2,2))';   % Euclidean distance
        end
end