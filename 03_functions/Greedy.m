function [xN_idx_out, max_err_history] = Greedy(xN_idx, DAx_act, Norm_Ax, Phi, SF_R, options)
% GREEDY_SELECTION Adds control points based on maximum deformation error.

arguments (Input)
    xN_idx               % Initial Point Indices
    DAx_act              % True Deformation at all nodes
    Norm_Ax              % Precomputed Distance Matrix
    Phi                  % Basis Function handle
    SF_R                 % Support Radius
    
    % Optional Name-Value pairs with default fail-safes
    options.err_tol (1,1) double = 0     % Defaults to 0 so it never triggers if uncalled
    options.K       (1,1) double = Inf   % Defaults to Inf (max iterations) if uncalled
end
arguments (Output)
    xN_idx_out            % Final list of Control Point Indices
    max_err_history      % Array tracking the maximum error at each step
end