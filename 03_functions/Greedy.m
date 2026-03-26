function [Nx_idx_out, max_err_history, rmse_history] = Greedy(Nx_idx_in, Ax, RAx, SF_R, options)
%GREEDY_SELECTION Adds control points based on maximum deformation error.
arguments (Input)
    Nx_idx_in            % Initial Point Indices
    Ax                   % Complete List of Aerofoil Points
    RAx                  % Applied Deformation (Displacement vectors) at all nodes
    SF_R                 % Support Radius
    
    options.err_tol (1,1) double = 0       % Error Tolerance Cutoff
    options.K       (1,1) double = Inf     % Number of Iterations Cutoff
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

% Total number of candidate nodes
M = size(Ax, 1);
all_idx = (1:M)';

% Initialize the evaluation matrix for the starting points
Dist_initial = Norm(Ax, Ax(Nx_idx_in, :));
CEx = Phi_WC2(Dist_initial / SF_R);
max_err = Inf; 

% Greedy Loop
while (max_err > options.err_tol) && (iterations < options.K)
    
    % Extract the known displacement values for the control points
    RNx = RAx(Nx_idx_out, :);                                       
    CNx = CEx(Nx_idx_out, :);              
    
    % Solve for the weights
    Gamma_x = CNx \ RNx(:, 1);             
    Gamma_y = CNx \ RNx(:, 2);             
    
    % --- MODIFIED SECTION START ---
    % Identify indices that are NOT currently control points
    eval_idx = setdiff(all_idx, Nx_idx_out);
    
    % Interpolate the displacement field ONLY at the unselected evaluation indices
    CEx_eval = CEx(eval_idx, :);
    Sx_x = CEx_eval * Gamma_x;                  
    Sx_y = CEx_eval * Gamma_y;                  
    
    % Calculate the error ONLY at the evaluation indices
    Err_x = RAx(eval_idx, 1) - Sx_x;              
    Err_y = RAx(eval_idx, 2) - Sx_y;              
    Err_mag = sqrt(Err_x.^2 + Err_y.^2);
    
    % Find maximum error among the unselected points
    [max_err, local_idx] = max(Err_mag);
    
    % Map the local index from Err_mag back to the global mesh array index
    next_idx = eval_idx(local_idx);
    
    % Calculate RMSE for the whole mesh.
    % Since error at control points is mathematically 0, the sum of squared 
    % errors for the evaluation points represents the total sum of squares.
    rmse_val = sqrt(sum(Err_x.^2 + Err_y.^2) / M);
    % --- MODIFIED SECTION END ---
    
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