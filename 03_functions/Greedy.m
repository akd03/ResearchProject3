function [Nx_idx_out, max_err_history] = Greedy(Nx_idx_in, Ax, RAx, SF_R, options)
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
end

%% FUNCTION BODY
Nx_idx_out = Nx_idx_in;
max_err_history = [];
iterations = 0;

% Initialize the evaluation matrix for the starting points
% Directly pass the 2D coordinate arrays
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
    
    % Interpolate the displacement field across the aerofoil
    Sx_x = CEx * Gamma_x;                  
    Sx_y = CEx * Gamma_y;                  
    
    % Calculate the error between the true displacement and the interpolated displacement
    Err_x = RAx(:, 1) - Sx_x;              
    Err_y = RAx(:, 2) - Sx_y;              
    Err_mag = sqrt(Err_x.^2 + Err_y.^2);
    
    % Find maximum error and its index
    [max_err, next_idx] = max(Err_mag);
    max_err_history(end+1, 1) = max_err;
    
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