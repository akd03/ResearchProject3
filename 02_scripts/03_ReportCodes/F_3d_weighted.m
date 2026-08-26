%% VARIABLE NAMING CONVENTION                                              
%  Nx       Known (Control) Points                                         
%  Nx_idx   Nx indices in Ax                                               
%  Ax       Aerofoil Mesh Points (3D)                                           
%  RAx      Applied Deformation (3D)                                            
%  DAx      Actual Aerofoil Position After Rx                              
clear; clc; close all;

%% 1. LOAD AND MERGE 3D MESH
filename = '05_meshes\01_3d_meshes\surfacepoints140K.plt'; % Update path if needed
fid = fopen(filename, 'r');
headerLines = 0;
while ~feof(fid)
    line = fgetl(fid);
    if ~isempty(regexp(line, '^\s*[-0-9.eE]+', 'once'))
        break; 
    end
    headerLines = headerLines + 1;
end
fclose(fid);
fprintf('Loading numeric data...\n');
RawData = readmatrix(filename, 'NumHeaderLines', headerLines, 'FileType', 'text');
Points = RawData(:, 1:3); 

% Split and Reshape (Bottom First, then Top)
bot_points = Points(1:2009, :);
top_points = Points(2010:end, :);

X_bot_full = reshape(bot_points(:,1), 41, 49);
Y_bot_full = reshape(bot_points(:,2), 41, 49);
Z_bot_full = reshape(bot_points(:,3), 41, 49);

X_top = reshape(top_points(:,1), 39, 48);
Y_top = reshape(top_points(:,2), 39, 48);
Z_top = reshape(top_points(:,3), 39, 48);

% Merge into a continuous C-grid (TE Bottom -> LE -> TE Top)
X_merged = [flipud(X_bot_full(:, 1:48)); X_top(2:end, :)]; 
Y_merged = [flipud(Y_bot_full(:, 1:48)); Y_top(2:end, :)];
Z_merged = [flipud(Z_bot_full(:, 1:48)); Z_top(2:end, :)];

% Flatten for the Greedy Solver
Ax = [X_merged(:), Y_merged(:), Z_merged(:)];
M = size(Ax, 1);

% Base Initial Point - Seed with a Trailing Edge point
Nx_idx_base = [1];
 
%% 2. DEFORMATION
% Assuming BendTwist3D takes (Points, RootOffset, TipOffset, span_length, chord_length)
DAx = BendTwist3D(Ax, [24.9, 0, -3.3], [47, 35.4, 1.6], 35.4, 47-24.9);
RAx = DAx - Ax;

%% 3. CALCULATE NORMALIZED CURVATURE
fprintf('Calculating 3D Principal Curvature...\n');
[k1_merged, k2_merged] = Curvature3D_Continuous(X_merged, Y_merged, Z_merged);
K_max_merged = max(abs(k1_merged), abs(k2_merged));

% Flatten and Normalize
kappa_flat = K_max_merged(:);
kappa_flat(isnan(kappa_flat)) = 0;
kappa_norm = kappa_flat / max(kappa_flat);

%% 4. PARAMETER SWEEP SETUP
SF_R = 100; % Reduced Support Radius for 3D Mesh Density
omega_vals = 0:0.01:1.0;    
num_omega = length(omega_vals);

% Setup Sweep 1: Fixed Node Count
N_vals = [100, 250, 500];
num_N_tests = length(N_vals);
s1_final_max_err = zeros(num_N_tests, num_omega);
s1_final_rmse = zeros(num_N_tests, num_omega); 

% Setup Sweep 2: Fixed Error Threshold 
err_threshold = 1e-4;
s2_total_nodes = zeros(1, num_omega);

% Setup Sweep 3: Visualizations
vis_N = 100;
vis_omega = [0, 0.5, 0.75];
vis_Nx_idx = cell(3, 1);

%% 5. DATA GENERATION LOOP
fprintf('Generating Data for Curvature Sweeps...\n');
for j = 1:num_omega
    omega = omega_vals(j);
    
    if mod(j, 10) == 0
        fprintf('  Evaluating Omega = %.2f\n', omega);
    end
    
    % --- SWEEP 1: Fixed Nodes ---
    for i = 1:num_N_tests
        N_target = N_vals(i);
        N_G = N_target - length(Nx_idx_base); 
        
        [~, max_err_hist, rmse_hist] = GreedyCholesky_Weighted(Nx_idx_base, Ax, RAx, SF_R, kappa_norm, omega, 'K', N_G);
        s1_final_max_err(i, j) = max_err_hist(end);
        s1_final_rmse(i, j) = rmse_hist(end); 
    end
    
    % --- SWEEP 2: Fixed Error Threshold ---
    [Nx_idx_out, ~, ~] = GreedyCholesky_Weighted(Nx_idx_base, Ax, RAx, SF_R, kappa_norm, omega, 'err_tol', err_threshold, 'K', Inf);
    s2_total_nodes(j) = length(Nx_idx_out);
end

% --- SWEEP 3: Visualization Points ---
fprintf('Generating Visualization Points (N=%d)...\n', vis_N);
N_G_vis = vis_N - length(Nx_idx_base);
for k = 1:3
    [Nx_idx_final, ~, ~] = GreedyCholesky_Weighted(Nx_idx_base, Ax, RAx, SF_R, kappa_norm, vis_omega(k), 'K', N_G_vis);
    vis_Nx_idx{k} = Nx_idx_final;
end
fprintf('Data generation complete. Plotting...\n');

%% 6. PLOT 1: PARAMETRIC SWEEPS (1x3 Layout)
fig_sweeps = figure('Name', '3D Curvature Weighting Sweeps', 'Position', [100, 100, 1400, 450]);
tiledlayout(1, 3, "TileSpacing", "compact");
colors = colororder;

% Find visualization indices
N_idx = find(N_vals == vis_N, 1);
omega_idx = zeros(1, length(vis_omega));
for v = 1:length(vis_omega)
    omega_idx(v) = find(abs(omega_vals - vis_omega(v)) < 1e-5, 1);
end
vis_color_idx = mod(N_idx-1, size(colors,1)) + 1;
vis_color = colors(vis_color_idx, :);

% --- TILE 1: Max Error vs Omega (Fixed N) ---
ax1 = nexttile;
hold(ax1, 'on'); grid(ax1, 'on');
set(ax1, 'YScale', 'log');
xlabel(ax1, 'Curvature Weighting (\omega)');
ylabel(ax1, 'Final Maximum Deformation Error');
title(ax1, 'Effect of \omega on Max Error');
for i = 1:num_N_tests
    c_idx = mod(i-1, size(colors,1)) + 1;
    plot(ax1, omega_vals, s1_final_max_err(i, :), '.-', 'LineWidth', 1.5, ...
        'MarkerSize', 8, 'Color', colors(c_idx, :), 'DisplayName', sprintf('N = %d', N_vals(i)));
end
vis_max_errors = s1_final_max_err(N_idx, omega_idx);
scatter(ax1, vis_omega, vis_max_errors, 60, vis_color, 'filled', 'DisplayName', sprintf('Visualised (N=%d)', vis_N));
legend(ax1, 'Location', 'best');
hold(ax1, 'off');

% --- TILE 2: RMSE vs Omega (Fixed N) ---
ax2 = nexttile;
hold(ax2, 'on'); grid(ax2, 'on');
set(ax2, 'YScale', 'log');
xlabel(ax2, 'Curvature Weighting (\omega)');
ylabel(ax2, 'Final RMSE');
title(ax2, 'Effect of \omega on Global RMSE');
for i = 1:num_N_tests
    c_idx = mod(i-1, size(colors,1)) + 1;
    plot(ax2, omega_vals, s1_final_rmse(i, :), '.-', 'LineWidth', 1.5, ...
        'MarkerSize', 8, 'Color', colors(c_idx, :), 'DisplayName', sprintf('N = %d', N_vals(i)));
end
vis_rmse_vals = s1_final_rmse(N_idx, omega_idx);
scatter(ax2, vis_omega, vis_rmse_vals, 60, vis_color, 'filled', 'DisplayName', sprintf('Vis Points (N=%d)', vis_N));
legend(ax2, 'Location', 'best');
hold(ax2, 'off');

% --- TILE 3: Required Nodes vs Omega (Fixed Tol) ---
ax3 = nexttile;
hold(ax3, 'on'); grid(ax3, 'on');
plot(ax3, omega_vals, s2_total_nodes, 'k.-', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel(ax3, 'Curvature Weighting (\omega)');
ylabel(ax3, 'Total Nodes Required');
title(ax3, sprintf('Nodes Required to Reach E_{tol} = %g', err_threshold));
hold(ax3, 'off');

%% 7. PLOT 2: VISUALIZATION OF SELECTED POINTS (X-Y Projection)
fig_vis = figure('Name', 'Spatial Clustering Visualization (X-Y Plane)', 'Position', [150, 200, 1200, 400]);
tiledlayout(1, 3, "TileSpacing", "compact");

for k = 1:3
    ax_curr = nexttile;
    hold(ax_curr, 'on'); grid(ax_curr, 'on'); axis(ax_curr, 'equal');
    
    % Extract the chosen node coordinates directly from the flattened array
    chosen_idx = vis_Nx_idx{k};
    control_points = Ax(chosen_idx, :); % Using undeformed coordinates for clarity
    
    % Pure scatter plot (top-down X-Y) without surface mesh
    scatter(ax_curr, control_points(:,1), control_points(:,2), 40, 'r', 'filled');
    
    title(ax_curr, sprintf('Control Nodes: N = %d, \\omega = %g', vis_N, vis_omega(k)));
    xlabel(ax_curr, 'Chordwise (X)'); ylabel(ax_curr, 'Spanwise (Y)');
    hold(ax_curr, 'off');
end

%% 8. SAVE PLOTS
if ~exist('06_results', 'dir')
    mkdir('06_results');
end
date_str = datestr(now, 'yyyymmdd');
save_filename_1 = sprintf('Test12_3D_OmegaSweeps_%s', date_str);
save_filename_2 = sprintf('Test13_3D_OmegaVis_%s', date_str);
savefig(fig_sweeps, fullfile('06_results', [save_filename_1, '.fig']));
savefig(fig_vis, fullfile('06_results', [save_filename_2, '.fig']));
fprintf('Plots successfully saved to /06_results/\n');

%% LOCAL FUNCTION: CONTINUOUS CURVATURE
function [k1, k2] = Curvature3D_Continuous(X, Y, Z)
    [Ni, Nj] = size(X);
    k1 = zeros(Ni, Nj);
    k2 = zeros(Ni, Nj);
    
    for i = 2:(Ni-1)
        for j = 2:(Nj-1)
            % Extract the 3x3 local coordinate neighborhood
            P = [reshape(X(i-1:i+1, j-1:j+1), 9, 1), ...
                 reshape(Y(i-1:i+1, j-1:j+1), 9, 1), ...
                 reshape(Z(i-1:i+1, j-1:j+1), 9, 1)];
            
            % Mean-shift
            P_shifted = P - [X(i,j), Y(i,j), Z(i,j)];
            
            % PCA Alignment
            [coeff, ~, ~] = pca(P_shifted);
            if det(coeff) < 0
                coeff(:, 2) = -coeff(:, 2);
            end
            
            % Rotate to local tangent space
            P_local = P_shifted * coeff;
            x_loc = P_local(:, 1);
            y_loc = P_local(:, 2);
            z_loc = P_local(:, 3);
            
            % Fit quadric surface
            A_mat = [0.5 * x_loc.^2, x_loc .* y_loc, 0.5 * y_loc.^2, x_loc, y_loc];
            c = A_mat \ z_loc;
            
            % Principal Curvatures
            A_coeff = c(1); B_coeff = c(2); C_coeff = c(3);
            T1 = (A_coeff + C_coeff) / 2;
            T2 = sqrt(T1^2 - A_coeff*C_coeff + B_coeff^2);
            
            k1(i, j) = T1 + T2;
            k2(i, j) = T1 - T2;
        end
    end
    
    % Extrapolate boundaries
    k1(1, :) = k1(2, :); k1(end, :) = k1(end-1, :);
    k2(1, :) = k2(2, :); k2(end, :) = k2(end-1, :);
    k1(:, 1) = k1(:, 2); k1(:, end) = k1(:, end-1);
    k2(:, 1) = k2(:, 2); k2(:, end) = k2(:, end-1);
    
    % Corners
    k1(1, 1) = k1(2, 2); k1(1, end) = k1(2, end-1);
    k1(end, 1) = k1(end-1, 2); k1(end, end) = k1(end-1, end-1);
    k2(1, 1) = k2(2, 2); k2(1, end) = k2(2, end-1);
    k2(end, 1) = k2(end-1, 2); k2(end, end) = k2(end-1, end-1);
end