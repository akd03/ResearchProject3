function [grad1, grad2] = CentralDifferenceDeriv(Ax, gamma)
%CENTRALDIFFERENCEDERIV Calculates the 1st and 2nd spatial derivatives w.r.t arc length
% for a closed 2D loop using a non-uniform central difference scheme.
%
% Input:
%   Ax: N x 2 matrix of coordinates [x, y]
%   gamma: (Optional) Relaxation factor for Laplacian smoothing. 
%          If omitted, the smoothing step is skipped.
% Outputs:
%   grad1: N x 2 matrix of first derivatives [dx/ds, dy/ds]
%   grad2: N x 2 matrix of second derivatives [d2x/ds2, d2y/ds2]

arguments
    Ax (:,2) double
    gamma double = [] % Empty default makes it optional and detectable
end

%% FUNCTION BODY
% Shift forward and back, double start/end point
Ax_f = [Ax(2:end, :); Ax(end, :)];     
Ax_b = [Ax(1, :); Ax(1:end-1, :)]; 
% Step size, forward and back
hf = sqrt((Ax_f(:, 1) - Ax(:, 1)).^2 + (Ax_f(:, 2) - Ax(:, 2)).^2);
hb = sqrt((Ax(:, 1) - Ax_b(:, 1)).^2 + (Ax(:, 2) - Ax_b(:, 2)).^2);
% Calc Denominators
den1 = hf .* hb .* (hf + hb) + 1e-12;
den2 = 0.5 .* den1;
% 1st Derivative
dx = (hb.^2 .* Ax_f(:, 1) - hf.^2 .* Ax_b(:, 1) + (hf.^2 - hb.^2) .* Ax(:, 1)) ./ den1;
dy = (hb.^2 .* Ax_f(:, 2) - hf.^2 .* Ax_b(:, 2) + (hf.^2 - hb.^2) .* Ax(:, 2)) ./ den1;
% 2nd Derivative
d2x = (hb .* Ax_f(:, 1) + hf .* Ax_b(:, 1) - (hf + hb) .* Ax(:, 1)) ./ den2;
d2y = (hb .* Ax_f(:, 2) + hf .* Ax_b(:, 2) - (hf + hb) .* Ax(:, 2)) ./ den2;
% Non-smoothed output
grad1 = [dx, dy];
grad2 = [d2x, d2y];


%% LAPLACIAN SMOOTHING
if ~isempty(gamma)
    % shift gradients forward and back
    grad1_f = [grad1(2:end, :); grad1(end, :)];
    grad1_b = [grad1(1, :); grad1(1:end-1, :)];
    grad2_f = [grad2(2:end, :); grad2(end, :)];
    grad2_b = [grad2(1, :); grad2(1:end-1, :)];
    % Apply the smoothing function
    grad1 = grad1 + gamma .* (grad1_f - 2.*grad1 + grad1_b);
    grad2 = grad2 + gamma .* (grad2_f - 2.*grad2 + grad2_b);
end


%% FINAL FORMATTING
grad1([1, end], :) = 0;
grad2([1, end], :) = 0;

end