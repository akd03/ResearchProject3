function [Nx_idx_out, max_err_history, rmse_history] = GreedyCholesky(Nx_idx_in, Ax, RAx, SF_R, options)
%GREEDYCHOLESKY Adds control points using an incremental Cholesky factorization.
arguments (Input)
    Nx_idx_in            % Initial Point Indices
    Ax                   % Complete List of Aerofoil/Wing Points
    RAx                  % Applied Deformation (Displacement vectors)
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
CNx = CEx(Nx_idx_out, :);
L = chol(CNx, 'lower');
% -------------------------------

% Greedy Loop
while (max_err > options.err_tol) && (iterations < options.K)
    
    % Extract the known displacement values (works for any number of columns)
    RNx = RAx(Nx_idx_out, :);                                       
    
    % --- VECTORIZED FORWARD & BACKWARD SUBSTITUTION ---
    % Solves for all spatial dimensions (X, Y, Z) simultaneously
    Y_mat = L \ RNx;
    Gamma = L' \ Y_mat;
    % ---------------------------------------------
    
    % --- LOGICAL INDEXING OPTIMIZATION ---
    is_eval = true(M, 1);
    is_eval(Nx_idx_out) = false; % Exclude active control points
    
    eval_idx = find(is_eval); % Get global indices of the evaluation points
    
    % Vectorized interpolation for all dimensions
    CEx_eval = CEx(is_eval, :);
    Sx = CEx_eval * Gamma;                  
    
    % Calculate the multidimensional error vector
    Err = RAx(is_eval, :) - Sx;    
    
    % Vectorized magnitude calculation (sqrt(x^2 + y^2 + z^2))
    Err_mag = sqrt(sum(Err.^2, 2));
    
    % Find maximum error among the unselected points
    [max_err, local_idx] = max(Err_mag);
    next_idx = eval_idx(local_idx);
    
    % Calculate global RMSE across all points and all dimensions
    rmse_val = sqrt(sum(Err.^2, 'all') / M);
    % -------------------------------------
    
    max_err_history(end+1, 1) = max_err;
    rmse_history(end+1, 1) = rmse_val;
    
    if (max_err > options.err_tol) && (iterations < options.K)
        
        new_dist = Norm(Ax, Ax(next_idx, :));
        new_CEx_col = Phi_WC2(new_dist / SF_R);
        
        % --- INCREMENTAL CHOLESKY UPDATE ---
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