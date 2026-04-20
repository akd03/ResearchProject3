function [Nx_idx_out, max_err_history, rmse_history] = Greedy_Standard(Nx_idx_in, Ax, RAx, SF_R, options)
%GREEDY_STANDARD Adds control points using full matrix recalculation (Baseline).
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
max_err = Inf; 

% Greedy Loop
while (max_err > options.err_tol) && (iterations < options.K)
    
    % 1. Extract coordinates and displacements for current control points
    Ax_c = Ax(Nx_idx_out, :);
    RNx = RAx(Nx_idx_out, :);                                       
    
    % 2. BUILD INTERPOLATION MATRIX FROM SCRATCH (k x k)
    Dist_c = Norm(Ax_c, Ax_c);
    CNx = Phi_WC2(Dist_c / SF_R);
    
    % 3. SOLVE SYSTEM FROM SCRATCH (Vectorized for all dimensions)
    Gamma = CNx \ RNx;             
    
    % 4. BUILD EVALUATION MATRIX FROM SCRATCH (M x k)
    Dist_eval = Norm(Ax, Ax_c);
    CEx = Phi_WC2(Dist_eval / SF_R);
    
    % 5. EVALUATE ENTIRE MESH AND CALCULATE ERROR
    Sx = CEx * Gamma;                  
    Err = RAx - Sx;              
    
    % Calculate magnitude of the error vector (works for 2D or 3D)
    Err_mag = sqrt(sum(Err.^2, 2));
    
    % Zero out the error at existing control points to prevent re-selection
    Err_mag(Nx_idx_out) = 0;
    
    % Find maximum error among the available points
    [max_err, next_idx] = max(Err_mag);
    
    % Calculate global RMSE across all points and all dimensions
    rmse_val = sqrt(sum(Err.^2, 'all') / M);
    
    max_err_history(end+1, 1) = max_err;
    rmse_history(end+1, 1) = rmse_val;
    
    if (max_err > options.err_tol) && (iterations < options.K)
        Nx_idx_out(end+1, 1) = next_idx;
        iterations = iterations + 1;
    end
end
end