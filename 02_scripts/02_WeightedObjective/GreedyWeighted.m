%% PARAMETER SWEEP: CURVATURE WEIGHTED GREEDY METHOD (2D)
clear; clc; close all;

%% 1. LOAD 2D MESH
MeshFile = fopen("05_meshes\NACA0012257x129.xyz", "r");                    
Header = fscanf(MeshFile, '%d', 2);                                        
Ni = Header(1);
Nj = Header(2);
MeshData = textscan(MeshFile, '%f %f %f');
fclose(MeshFile);

% Extract just the Aerofoil Surface Points (2D: X and Y)
Mx = [reshape(MeshData{1}, [Ni*Nj, 1]), reshape(MeshData{3}, [Ni*Nj, 1])];
Ax = [Mx(1:Ni, 1), Mx(1:Ni, 2)];
M = size(Ax, 1);

% Base Initial Point - Seed with the Trailing Edge (Index 1)                                       
Nx_idx_base = [1];
 
%% 2. DEFORMATION
% Calls your untouched NACACamber function
DAx = NACACamber(Ax, 0.02, 0.5);
RAx = DAx - Ax;

%% 3. CALCULATE NORMALIZED CURVATURE (CENTRAL DIFFERENCES)
% MATLAB's gradient function computes 2nd-order central differences for interior points
dX = gradient(Ax(:, 1));
dY = gradient(Ax(:, 2));
ddX = gradient(dX);
ddY = gradient(dY);

% Curvature formula: k = |x'y'' - y'x''| / (x'^2 + y'^2)^(3/2)
Numerator = abs(dX .* ddY - dY .* ddX);
Denominator = (dX.^2 + dY.^2).^(3/2);

% Prevent division by zero at perfectly straight segments
Denominator(Denominator == 0) = eps;

kappa = Numerator ./ Denominator;
kappa(isnan(kappa)) = 0; % Strip NaN noise at any stagnant points

% Normalize curvature strictly between 0 and 1 for the objective function
kappa_norm = kappa / max(kappa);

%% 4. EXPERIMENTAL SETUP
SF_R = 3;                  % 2D Support Radius
omega_vals = 0:0.01:1.0;    % Sweep weighting from 0 (Pure Error) to 1 (Pure Curvature)
num_tests = length(omega_vals);

% Test 1 Parameters: Fixed Node Count
fixed_nodes = 50;
final_rmse_vals = zeros(1, num_tests);

% Test 2 Parameters: Fixed Error Threshold
err_threshold = 1e-4;
final_node_counts = zeros(1, num_tests);

fprintf('Starting Omega Sweep (0 to 1)...\n');

%% 5. EXECUTE SWEEP
for i = 1:num_tests
    omega = omega_vals(i);
    
    % --- TEST 1: Fixed Nodes (Track RMSE) ---
    [~, ~, rmse_history] = GreedyCholesky_Weighted(Nx_idx_base, Ax, RAx, SF_R, kappa_norm, omega, "K", fixed_nodes);
    final_rmse_vals(i) = rmse_history(end);
    
    % --- TEST 2: Fixed Error Threshold (Track Node Count) ---
    [Nx_idx_out, ~, ~] = GreedyCholesky_Weighted(Nx_idx_base, Ax, RAx, SF_R, kappa_norm, omega, "err_tol", err_threshold);
    final_node_counts(i) = length(Nx_idx_out);
    
    fprintf('Omega = %.1f completed.\n', omega);
end

%% 6. PLOT RESULTS
figure;
tiledlayout(1, 2, "TileSpacing", "compact");

% Plot 1: Omega vs Final RMSE (Fixed Nodes)
ax1 = nexttile;
hold(ax1, 'on'); grid(ax1, 'on');
plot(ax1, omega_vals, final_rmse_vals, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 2);
xlabel(ax1, 'Curvature Weighting (\omega)');
ylabel(ax1, sprintf('Final RMSE (at N = %d)', fixed_nodes));
title(ax1, 'Effect of \omega on Interpolation Accuracy');
set(ax1, 'YScale', 'log');
hold(ax1, 'off');

% Plot 2: Omega vs Final Node Count (Fixed Error Threshold)
ax2 = nexttile;
hold(ax2, 'on'); grid(ax2, 'on');
plot(ax2, omega_vals, final_node_counts, 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'MarkerSize', 2);
xlabel(ax2, 'Curvature Weighting (\omega)');
ylabel(ax2, sprintf('Total Nodes Required (Error < %g)', err_threshold));
title(ax2, 'Effect of \omega on Data Reduction Efficiency');
hold(ax2, 'off');
figure; hold on; axis equal;
plot(Ax(:,1), Ax(:,2), 'b-');
plot(DAx(:,1), DAx(:,2), 'r-');
