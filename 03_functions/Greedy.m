function [Nx_idx_out, max_err_history] = Greedy(Nx_idx, DAx, Phi_Ax, options)
%GREEDY_SELECTION Adds control points based on maximum deformation error.
%   Inputs: Nx_idx, DAx, Phi_Ax, options (err_tol, K)
arguments (Input)
    Nx_idx               % Initial Point Indices
    DAx                  % True Deformation at all nodes
    Phi_Ax               % Precomputed Influence Matrix
    options.err_tol (1,1) double = 0     % Defaults to 0 so it never triggers if uncalled
    options.K       (1,1) double = Inf   % Defaults to Inf (max iterations) if uncalled
end
arguments (Output)
    Nx_idx_out            % Final list of Control Point Indices
    max_err_history      % Array tracking the maximum error at each step
end

%% FUNCTION BODY
% Initialise Counter/Output
Nx_idx_out = Nx_idx;
max_err_history = [];
iterations = 0;
% Greedy Loop
while true
    % Def, influence from Norm Matrix and True Def
    DxN = DAx(Nx_idx_out, :);              % Nc x 2                         
    CEx = Phi_Ax(:, Nx_idx_out);           % Ne x Nc
    CNx = CEx(Nx_idx_out, :);              % Nc x Nc
    % Weights of Known xN
    Gamma_x = CNx \ DxN(:, 1);             % Nc x 1
    Gamma_y = CNx \ DxN(:, 2);             % Nc x 1 
    % Interpolated Def
    Sx_x = CEx * Gamma_x;                  % Ne x 1
    Sx_y = CEx * Gamma_y;                  
    % Positional Error
    Err_x = DAx(:, 1) - Sx_x;              % Ne x 1
    Err_y = DAx(:, 2) - Sx_y;              
    Err_mag = sqrt(Err_x.^2 + Err_y.^2);
    % Find Max Error and Location
    [max_err, next_idx] = max(Err_mag);
    max_err_history(end+1, 1) = max_err;
    % End Loop At Stopping Criteria
    if max_err <= options.err_tol || iterations >= options.K
        break;
    end
    % Add Point w/ Max Error To Nx
    Nx_idx_out(end+1, 1) = next_idx;
    iterations = iterations + 1;
end
end