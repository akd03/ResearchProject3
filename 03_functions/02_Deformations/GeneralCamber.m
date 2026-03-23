function DAx = GeneralCamber(Ax, yc_func, dyc_func)
% GENERALCAMBER Deforms a symmetric aerofoil mesh along an arbitrary camber line
% while strictly preserving orthogonal thickness.
%
% Inputs:
%   Ax: N x 2 matrix of original symmetric aerofoil coordinates [x, y]
%   yc_func: Function handle for the camber line y_c(x)
%   dyc_func: Function handle for the camber line derivative dy_c/dx
%
% Output:
%   DAx: N x 2 matrix of deformed coordinates [\hat{x}, \hat{y}]

    % Extract original coordinates
    x_orig = Ax(:, 1);
    y_orig = Ax(:, 2); % This inherently contains the +/- thickness

    % 1. Evaluate camber line and gradient at the specific x coordinates
    y_c = yc_func(x_orig);
    dy_c = dyc_func(x_orig);

    % 2. Calculate local tangent angle
    theta = atan(dy_c);

    % 3. Apply orthogonal thickness projection
    x_def = x_orig - y_orig .* sin(theta);
    y_def = y_c + y_orig .* cos(theta);

    % Return deformed mesh
    DAx = [x_def, y_def];
end