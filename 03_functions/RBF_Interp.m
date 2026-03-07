function [Sx] = RBF_Interp(xN, Gamma, EN, R, Phi, Norm)
%RFB_Interp - Interpolates Between Points With An RBF
%   Inputs: x_N, F_x_N, N, C_x, Gamma, E_N, R, Phi, Norm
%   Must First Run Script To Load phi and norm Functions
%   This function interpolates between known points 
arguments (Input)
    xN           % Position of Known Points
    Gamma        % Influence Coefficients
    EN           % Evaluation Points
    R            % Support Function Radius
    Phi          % Basis Function
    Norm         % Euclidean Norm
end
arguments (Output)
    Sx          % Value of Evaluation Points
end

%% Function Body
% Finding Dependence Matrix A
[xN_N,xN_D] = size(xN);
[EN_N, EN_D] = size(EN);                    % E_N_D should equal x_N_D
A_rows = reshape(xN, [xN_N, 1, xN_D]);
A_cols = reshape(EN, [1, EN_N, EN_D]);
norm = Norm(A_rows, A_cols);
A = Phi(norm/R);
% Solving for Interpolated Function S(x)
Sx = A' * Gamma;
end