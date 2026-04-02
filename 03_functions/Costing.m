function [cost_fps, cost_greedy_old, cost_greedy_new, total_old, total_new] = Costing(M, N_IP, N_G)
%COSTING Computes theoretical operations for RBF point selection.
%
% Inputs:
%   M    - Number of candidate nodes in the surface mesh
%   N_IP - Number of Initial Points selected via FPS
%   N_G  - Number of points selected via the Greedy Algorithm
%
% Outputs:
%   cost_fps        - Operations for Farthest Point Sampling
%   cost_greedy_old - Operations for the standard O(N^4) greedy solver
%   cost_greedy_new - Operations for the O(N^3) Cholesky greedy solver
%   total_old       - FPS + Old Greedy
%   total_new       - FPS + New Greedy

    % 1. FPS Cost: O(M * N_IP)
    cost_fps = M * N_IP;
    
    N = N_IP + N_G;

    % 2. Greedy Costs
    if N_G > 0
        k = (N_IP + 1):N;
        
        % Original Method: k^3 (Matrix Solve) + M*k (Full Field Eval)
        cost_greedy_old = sum(k.^3 + M .* k);
        
        % New Method: Setup chol() + k^2 (Cholesky Update) + (M-k)*k (Subset Field Eval)
        cost_greedy_new = (N_IP^3) + sum(k.^2 + (M - k) .* k);
    else
        % If no greedy iterations run, original cost is zero.
        cost_greedy_old = 0;
        
        % New method still incurs the one-off cost to build the final matrix.
        cost_greedy_new = (N_IP^3); 
    end
    
    % 3. Total Costs
    total_old = cost_fps + cost_greedy_old;
    total_new = cost_fps + cost_greedy_new;
end