function [D] = Norm(A, B)
%NORM Calculates the pairwise Euclidean distance matrix.
%   A: N x d matrix of coordinates (e.g., N x 2 for 2D)
%   B: M x d matrix of coordinates
%   Outputs an N x M matrix of distances.

arguments (Input)
    A (:,:) double        % N x d matrix of coordinates (d = dimension)
    B (:,:) double        % M x d matrix of coordinates 
end
arguments (Output)
    D (:,:) double        % N x M 
end

%% FUNCTION BODY
A_sq = sum(A.^2, 2);                     % Find A^2 and B^2
B_sq = sum(B.^2, 2)'; 
AB = A * B';                             % Find AB
D = sqrt(max(A_sq + B_sq - 2 * AB, 0));  % Expansion of Norm Formula

end