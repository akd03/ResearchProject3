function [total_std, total_inc] = CostTot(N_s, N_ip, N_g)
% COSTTOT Calculates the total cumulative FLOPs for the selection process.
% Inputs:
%   N_s  - Total number of candidate surface mesh points
%   N_ip - Number of initial points selected via FPS
%   N_g  - Number of greedy points selected
% Outputs:
%   total_std - Total FLOPs using the Standard Greedy method
%   total_inc - Total FLOPs using the Incremental Greedy method

    % Get the per-step costs
    [cost_fps, cost_std, cost_inc] = CostStep(N_s, N_ip, N_g);

    % Total Cost = Sum(FPS) + Init_Chol + Sum(Greedy)
    total_std = sum(cost_fps) + sum(cost_std);
    total_inc = sum(cost_fps) + sum(cost_inc);
end