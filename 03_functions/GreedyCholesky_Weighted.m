function [Nx_idx_out, max_err_history, rmse_history] = GreedyCholesky_Weighted(Nx_idx_in, Ax, RAx, SF_R, kappa_norm, omega, options)
%GREEDYCHOLESKY_WEIGHTED Adds control points using a mixed Error-Curvature objective.
%
% This solver differs from the standard GreedyCholesky by utilizing a weighted 
% objective function to select the next point, balancing the maximum physical 
% displacement error against the normalized geometric curvature.
%
% The err_tol stopping criteria is strictly evaluated against the true physical
% displacement error, independent of the curvature weighting.

arguments (Input)
    Nx_idx_in            % Initial Point Indices
    Ax                   % Complete List of Aerofoil/Wing Points
    RAx                  % Applied Deformation (Displacement vectors)
    SF_R                 % Support Radius
    kappa_norm           % Pre-calculated, normalized curvature array (M x 1)
    omega                % Weighting parameter (0 = Pure Error, 1 = Pure Curvature)
    
    options.err_tol (1,1) double = 0     
    options.K       (1,1) double = Inf   
end
arguments (Output)
    Nx_idx_out           % Final list of Control Point Indices
    max_err_history      % Array tracking the maximum physical error at each step
    rmse_history         % Array tracking the physical RMSE at each step
end

%% FUNCTION BODY
Nx_idx_out = Nx_idx_in;
max_err_history = [];
rmse_history = [];
iterations = 0;
M = size(Ax, 1); 

% Initialize the evaluation matrix for the starting points
Dist_initial = Norm(Ax, Ax(Nx_idx_in, :));
CEx = Phi_WC2(Dist_initial / SF_R);
current_max_err = Inf; 

% --- CHOLESKY INITIALIZATION ---
CNx = CEx(Nx_idx_out, :);
L = chol(CNx, 'lower');
% -------------------------------

% Greedy Loop
while (current_max_err > options.err_tol) && (iterations < options.K)
    
    % Extract the known displacement values
    RNx = RAx(Nx_idx_out, :);                                       
    
    % --- VECTORIZED FORWARD & BACKWARD SUBSTITUTION ---
    Y_mat = L \ RNx;
    Gamma = L' \ Y_mat;
    % ---------------------------------------------
    
    % --- LOGICAL INDEXING OPTIMIZATION ---
    is_eval = true(M, 1);
    is_eval(Nx_idx_out) = false; % Strictly exclude active control points
    
    eval_idx = find(is_eval); % Get global indices of the evaluation points
    
    % Vectorized interpolation
    CEx_eval = CEx(is_eval, :);
    Sx = CEx_eval * Gamma;                  
    
    % Calculate the true physical error vector
    Err = RAx(is_eval, :) - Sx;    
    Err_mag = sqrt(sum(Err.^2, 2));
    
    % Pure displacement metrics (Used strictly for tracking and stopping)
    current_max_err = max(Err_mag);
    rmse_val = sqrt(sum(Err.^2, 'all') / M);
    % -------------------------------------
    
    max_err_history(end+1, 1) = current_max_err;
    rmse_history(end+1, 1) = rmse_val;
    
    % Check stopping condition before adding a new point
    if (current_max_err > options.err_tol) && (iterations < options.K)
        
        % --- THE MIXED OBJECTIVE FUNCTION ---
        % 1. Normalize the current iteration's error field
        % (eps prevents division by zero if error perfectly hits 0)
        E_norm = Err_mag / max(current_max_err, eps);
        
        % 2. Extract the curvature specifically for the unselected points
        K_eval = kappa_norm(is_eval);
        
        % 3. Calculate Objective: F = (1-w)*E + w*k
        Objective = (1 - omega) * E_norm + omega * K_eval;
        
        % 4. Select the next point based on the highest objective score
        [~, local_idx] = max(Objective);
        next_idx = eval_idx(local_idx);
        % ------------------------------------
        
        % --- INCREMENTAL CHOLESKY UPDATE ---
        new_dist = Norm(Ax, Ax(next_idx, :));
        new_CEx_col = Phi_WC2(new_dist / SF_R);
        
        c = new_CEx_col(Nx_idx_out); 
        v = 1.0; 
        
        l_vec = L \ c;
        gamma = sqrt(v - (l_vec' * l_vec));
        
        k_current = size(L, 1);
        L = [L, zeros(k_current, 1); l_vec', gamma];
        % -----------------------------------
        
        Nx_idx_out(end+1, 1) = next_idx;
        CEx = [CEx, new_CEx_col];
        
        iterations = iterations + 1;
    end
end
end