function [S_x] = Def_RBF_Interpolator(x_N, Gamma, E_N, R, Phi, Norm)
%RFB_Interpolator_NoPoly  - Interpolates Between Points With An RBF
%   Inputs: x_N, F_x_N, N, C_x, Gamma, E_N, R, Phi, Norm
%   Must First Run Script To Load phi and norm Functions
%   This function interpolates between known points 
arguments (Input)
    x_N          % Position of Known Points
    Gamma        % Influence Coefficients
    E_N          % Evaluation Points
    R            % Support Function Radius
    Phi          % Basis Function
    Norm         % Euclidean Norm
end
arguments (Output)
    S_x          % Value of Evaluation Points
end

%% Function Body
% Finding Dependence Matrix A
[x_N_N,x_N_D] = size(x_N);
[E_N_N, E_N_D] = size(E_N);                    % E_N_D should equal x_N_D
A_rows = reshape(x_N, [x_N_N, 1, x_N_D]);
A_cols = reshape(E_N, [1, E_N_N, E_N_D]);
norm = Norm(A_rows, A_cols);
A = Phi(norm/R);
% Solving for Interpolated Function S(x)
S_x = A' * Gamma;
end