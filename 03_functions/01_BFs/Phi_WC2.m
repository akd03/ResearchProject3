function [Phi] = Phi_WC2(r)
%PHI_WC2 Evaluates the Wendland C2 radial basis function on Matrix
arguments (Input)
    r double % Normalized distance array (scalar, vector, or matrix)
end
arguments (Output)
    Phi      % Evaluated basis function values of the same size
end

%% FUNCTION BODY
Phi = (max(0, 1 - abs(r))).^4 .* (4 * abs(r) + 1);
end