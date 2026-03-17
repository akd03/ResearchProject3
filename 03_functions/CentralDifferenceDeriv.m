function [grad1, grad2] = CentralDifferenceDeriv(Ax)
%DERIVATIVESLOOP Calculates the 1st and 2nd spatial derivatives w.r.t arc length
% for a closed 2D loop using a non-uniform central difference scheme.
%
% Input:
%   Ax: N x 2 matrix of coordinates [x, y]
% Outputs:
%   grad1: N x 2 matrix of first derivatives [dx/ds, dy/ds]
%   grad2: N x 2 matrix of second derivatives [d2x/ds2, d2y/ds2]

arguments
    Ax (:,2) double
end

%% FUNCTION BODY
% Forward and Back shifts, open loop
Ax_f = [Ax(2:end, :); Ax(end, :)];     
Ax_b = [Ax(1, :); Ax(1:end-1, :)]; 
% Forward and Backward arc length step sizes
hf = sqrt((Ax_f(:, 1) - Ax(:, 1)).^2 + (Ax_f(:, 2) - Ax(:, 2)).^2);
hb = sqrt((Ax(:, 1) - Ax_b(:, 1)).^2 + (Ax(:, 2) - Ax_b(:, 2)).^2);
% Handle the boundary nodes (Index 1 and End)
den1 = hf .* hb .* (hf + hb) + 1e-12;
den2 = 0.5 .* den1; 

dx = (hb.^2 .* Ax_f(:, 1) - hf.^2 .* Ax_b(:, 1) + (hf.^2 - hb.^2) .* Ax(:, 1)) ./ den1;
dy = (hb.^2 .* Ax_f(:, 2) - hf.^2 .* Ax_b(:, 2) + (hf.^2 - hb.^2) .* Ax(:, 2)) ./ den1;

% 5. Second Derivative (grad2)
d2x = (hb .* Ax_f(:, 1) + hf .* Ax_b(:, 1) - (hf + hb) .* Ax(:, 1)) ./ den2;
d2y = (hb .* Ax_f(:, 2) + hf .* Ax_b(:, 2) - (hf + hb) .* Ax(:, 2)) ./ den2;

% 6. Format the outputs
grad1 = [dx, dy];
grad2 = [d2x, d2y];

% 7. Zero out the boundaries so they are ignored by the sorter
grad1([1, end], :) = 0;
grad2([1, end], :) = 0;