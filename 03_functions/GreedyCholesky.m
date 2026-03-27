function [Nx_idx_out, max_err_history, rmse_history] = GreedyCholesky(Nx_idx_in, Ax, RAx, SF_R, options)
%GREEDYCHOLESKY Adds control points using an incremental Cholesky factorization.
arguments (Input)
    Nx_idx_in            % Initial Point Indices
    Ax                   % Complete List of Aerofoil Points
    RAx                  % Applied Deformation (Displacement vectors) at all nodes
    SF_R                 % Support Radius
    
    options.err_tol (1,1) double = 0     
    options.K       (1,1) double = Inf   
end
arguments (Output)
    Nx_idx_out           % Final list of Control Point Indices
    max_err_history      % Array tracking the maximum error at each step
    rmse_history         % Array tracking the root mean square error at each step
end

%% FUNCTION BODY
Nx_idx_out = Nx_idx_in;
max_err_history = [];
rmse_history = [];
iterations = 0;
M = size(Ax, 1); % Total number of nodes

% Initialize the evaluation matrix for the starting points
Dist_initial = Norm(Ax, Ax(Nx_idx_in, :));
CEx = Phi_WC2(Dist_initial / SF_R);
max_err = Inf; 

% --- CHOLESKY INITIALIZATION ---
% Extract the initial k x k interpolation matrix
CNx = CEx(Nx_idx_out, :);
% Perform a standard Cholesky decomposition for the initial points (L * L^T = CNx)
L = chol(CNx, 'lower');
% -------------------------------

% Greedy Loop
while (max_err > options.err_tol) && (iterations < options.K)
    
    % Extract the known displacement values for the control points
    RNx = RAx(Nx_idx_out, :);                                       
    
    % --- FORWARD & BACKWARD SUBSTITUTION SOLVE ---
    % MATLAB's '\' operator automatically detects lower/upper triangular matrices 
    % and applies optimal forward/backward substitution.
    
    % 1. Forward Substitution: Solve L * Y = R
    Y_x = L \ RNx(:, 1);
    Y_y = L \ RNx(:, 2);
    
    % 2. Backward Substitution: Solve L^T * Gamma = Y
    Gamma_x = L' \ Y_x;
    Gamma_y = L' \ Y_y;
    % ---------------------------------------------
    
    % --- LOGICAL INDEXING OPTIMIZATION ---
    is_eval = true(M, 1);
    is_eval(Nx_idx_out) = false; % Exclude active control points
    
    eval_idx = find(is_eval); % Get global indices of the evaluation points
    
    % Interpolate the displacement field ONLY at the unselected indices
    CEx_eval = CEx(is_eval, :);
    Sx_x = CEx_eval * Gamma_x;                  
    Sx_y = CEx_eval * Gamma_y;                  
    
    % Calculate the error ONLY at the evaluation indices
    Err_x = RAx(is_eval, 1) - Sx_x;              
    Err_y = RAx(is_eval, 2) - Sx_y;              
    Err_mag = sqrt(Err_x.^2 + Err_y.^2);
    
    % Find maximum error among the unselected points
    [max_err, local_idx] = max(Err_mag);
    
    % Map the local index from Err_mag back to the global mesh array index
    next_idx = eval_idx(local_idx);
    
    % Calculate global RMSE
    rmse_val = sqrt(sum(Err_x.^2 + Err_y.^2) / M);
    % -------------------------------------
    
    max_err_history(end+1, 1) = max_err;
    rmse_history(end+1, 1) = rmse_val;
    
    if (max_err > options.err_tol) && (iterations < options.K)
        
        % Incrementally calculate the single new influence column
        new_dist = Norm(Ax, Ax(next_idx, :));
        new_CEx_col = Phi_WC2(new_dist / SF_R);
        
        % --- INCREMENTAL CHOLESKY UPDATE ---
        % Extract the influence vector 'c' mapping the new point to existing control points
        c = new_CEx_col(Nx_idx_out); 
        
        % v is the self-influence scalar. Phi_WC2(0) = 1.
        v = 1.0; 
        
        % Solve L * l = c using forward substitution
        l_vec = L \ c;
        
        % Calculate gamma
        gamma = sqrt(v - (l_vec' * l_vec));
        
        % Expand the lower triangular matrix L
        k_current = size(L, 1);
        L = [L, zeros(k_current, 1); l_vec', gamma];
        % -----------------------------------
        
        % Append the new point and column to the main arrays
        Nx_idx_out(end+1, 1) = next_idx;
        CEx = [CEx, new_CEx_col];
        
        iterations = iterations + 1;
    end
end
end