function DAx = NACACamber(Ax, m, p)
% NACACAMBER Deforms a symmetric aerofoil mesh using the exact NACA 4-digit
% camber line equations and orthogonal thickness projection.
%
% Inputs:
%   Ax: N x 2 matrix of symmetric aerofoil coordinates [x, y]
%   m:  Maximum camber (e.g., 0.05 for NACA 5xxx, 0.02 for NACA 2xxx)
%   p:  Position of maximum camber (e.g., 0.4 for NACA x4xx)
%
% Output:
%   DAx: N x 2 matrix of deformed coordinates [x_def, y_def]

arguments
    Ax (:,2) double
    m (1,1) double
    p (1,1) double
end

%% SYMMETRIC BYPASS
% Prevent division by zero if a symmetric aerofoil is requested
if m == 0 || p == 0
    DAx = Ax;
    return;
end

%% PIECEWISE CAMBER FUNCTIONS
% Camber line equation
yc_func = @(x) (x < p) .* (m/p^2 .* (2*p.*x - x.^2)) + ...
               (x >= p) .* (m/(1-p)^2 .* ((1-2*p) + 2*p.*x - x.^2));

% Gradient equation
dyc_func = @(x) (x < p) .* (2*m/p^2 .* (p - x)) + ...
                (x >= p) .* (2*m/(1-p)^2 .* (p - x));

%% ORTHOGONAL PROJECTION
% Pass the functions to the general camber mapping routine
DAx = GeneralCamber(Ax, yc_func, dyc_func);

end