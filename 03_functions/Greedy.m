function [Nx_idx_out, max_err_history, rmse_history] = Greedy(Nx_idx_in, Ax, RAx, SF_R, options)
%GREEDY_SELECTION Adds control points based on maximum deformation error.
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

% Greedy Loop
while (max_err > options.err_tol) && (iterations < options.K)
    
    % Extract the known displacement values for the control points
    RNx = RAx(Nx_idx_out, :);                                       
    CNx = CEx(Nx_idx_out, :);              
    
    % Solve for the weights (k x k system)
    Gamma_x = CNx \ RNx(:, 1);             
    Gamma_y = CNx \ RNx(:, 2);             
    
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
        Nx_idx_out(end+1, 1) = next_idx;
        
        % Incrementally calculate the single new influence column
        new_dist = Norm(Ax, Ax(next_idx, :));
        new_CEx_col = Phi_WC2(new_dist / SF_R);
        
        % Append the new column to the existing evaluation matrix
        CEx = [CEx, new_CEx_col];
        
        iterations = iterations + 1;
    end
end
end