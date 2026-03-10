function [Nx_idx_out] = IP_Distance(Nx_idx_in, Ax, K)
%IP_DISTANCE Select K initial points based on Euclidean distance.
%   Inputs: Nx_idx_in, Ax, K
%   Uses Farthest Point Sampling (FPS) to append K new indices 
%   to the existing Nx_idx_in array.

arguments (Input)
    Nx_idx_in            % Current List of Initial Point Indices (e.g., LE and TE)
    Ax                   % Complete List of Aerofoil Points
    K                    % Number of Initial Points to Select
end
arguments (Output)
    Nx_idx_out           % New List of Initial Point Indices
end

%% FUNCTION BODY
% Preallocate output matrix
N = length(Nx_idx_in);
Nx_idx_out = zeros(N + K, 1);       
Nx_idx_out(1:N) = Nx_idx_in;
% Calc distance from surf mesh points to intial control points in list
Ax_MinDist = min(Norm(Ax, Ax(Nx_idx_in, :)), [], 2);
% Iterative selection loop
for i = 1:K
    [~, next_idx] = max(Ax_MinDist);                           % Find idx that maximizes the min dist
    Nx_idx_out(N + i) = next_idx;                              % Add to Nx list
    Ax_MinDist = min(Ax_MinDist, Norm(Ax, Ax(next_idx, :)));   % update running min distances
end

end