%% VARIABLE NAMING CONVENTION                                              
%  Nx       Known (Control) Points                                         
%  Nx_idx   Nx indices in Ax                                               
%  Ax       Aerofoil Mesh Points                                           
%  RAx      Applied Deformation                                            
%  DAx      Actual Aerofoil Position After Rx                              

clear; clc; close all;

%% 1. LOAD 2D MESH
MeshFile = fopen("05_meshes\NACA0012257x129.xyz", "r");                    
Header = fscanf(MeshFile, '%d', 2);                                        
Ni = Header(1);
Nj = Header(2);
MeshData = textscan(MeshFile, '%f %f %f');
fclose(MeshFile);

% Extract Aerofoil Surface Points (2D)
Mx = [reshape(MeshData{1}, [Ni*Nj, 1]), reshape(MeshData{3}, [Ni*Nj, 1])];
Ax = [Mx(1:Ni, 1), Mx(1:Ni, 2)];
M = size(Ax, 1);

% Base Initial Point - Seed with the Trailing Edge (Index 1)                                       
Nx_idx_base = [1];
 
%% 2. DEFORMATION (Rigid Placeholder)
DAx = NACACamber(Ax, 0.15, 0.5);
RAx = DAx - Ax;

%% 3. CALCULATE NORMALIZED CURVATURE (CENTRAL DIFFERENCES)
dX = gradient(Ax(:, 1));
dY = gradient(Ax(:, 2));
ddX = gradient(dX);
ddY = gradient(dY);

Numerator = abs(dX .* ddY - dY .* ddX);
Denominator = (dX.^2 + dY.^2).^(3/2);
Denominator(Denominator == 0) = eps; 

kappa = Numerator ./ Denominator;
kappa(isnan(kappa)) = 0; 
kappa_norm = kappa / max(kappa); 

%% 4. PARAMETER SWEEP SETUP
SF_R = 3;                  
omega_vals = 0:0.01:1.0;    
num_omega = length(omega_vals);

% Setup Sweep 1: Fixed Node Count (Track Max Error and RMSE)
N_vals = [25, 50, 100];
num_N_tests = length(N_vals);
s1_final_max_err = zeros(num_N_tests, num_omega);
s1_final_rmse = zeros(num_N_tests, num_omega); % NEW ARRAY

% Setup Sweep 2: Fixed Error Threshold (Track Node Count)
err_threshold = 1e-6;
s2_total_nodes = zeros(1, num_omega);

% Setup Sweep 3: Visualizations
vis_N = 25;
vis_omega = [0, 0.4, 0.75];
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
        s1_final_rmse(i, j) = rmse_hist(end); % STORE RMSE
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
% Increased figure width to accommodate the third plot cleanly
fig_sweeps = figure('Name', 'Curvature Weighting Sweeps', 'Position', [100, 100, 1400, 450]);
tiledlayout(1, 3, "TileSpacing", "compact");

colors = colororder;

% --- Find the visualization points indices once for both plots ---
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
scatter(ax1, vis_omega, vis_max_errors, 60, vis_color, 'filled', ...
        'LineWidth', 1, 'DisplayName', sprintf('Visualised (N=%d)', vis_N));

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
scatter(ax2, vis_omega, vis_rmse_vals, 60, vis_color, 'filled', ...
    'LineWidth', 1, 'DisplayName', sprintf('Vis Points (N=%d)', vis_N));

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

%% 7. PLOT 2: VISUALIZATION OF SELECTED POINTS (1x4 Layout)
fig_vis = figure('Name', 'Spatial Clustering Visualization', 'Position', [150, 200, 1600, 350]);
tiledlayout(1, 4, "TileSpacing", "compact");

% --- TILE 1: Undeformed Geometry ---
ax4 = nexttile;
hold(ax4, 'on'); grid(ax4, 'on'); axis(ax4, 'equal');
plot(ax4, Ax(:,1), Ax(:,2), 'k-', 'LineWidth', 1.5);
title(ax4, 'Undeformed Aerofoil');
xlabel(ax4, 'X'); ylabel(ax4, 'Y');
hold(ax4, 'off');
xlim([-0.1, 1.1]); ylim([-0.3, 0.3]);

% --- TILES 2, 3, 4: Deformed with Control Points ---
for k = 1:3
    ax_curr = nexttile;
    hold(ax_curr, 'on'); grid(ax_curr, 'on'); axis(ax_curr, 'equal');
    
    plot(ax_curr, DAx(:,1), DAx(:,2), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, 'DisplayName', 'Deformed Profile');
    
    chosen_idx = vis_Nx_idx{k};
    control_points = DAx(chosen_idx, :);
    
    scatter(ax_curr, control_points(:,1), control_points(:,2), 30, 'r', 'filled', 'DisplayName', 'Control Points');
    
    title(ax_curr, sprintf('Deformed: N = %d, \\omega = %g', vis_N, vis_omega(k)));
    xlabel(ax_curr, 'X'); ylabel(ax_curr, 'Y');
    xlim([-0.1, 1.1]); ylim([-0.3, 0.3]);
    if k == 3
        legend(ax_curr, 'Location', 'northeast');
    end
    hold(ax_curr, 'off');
end

%% 8. SAVE PLOTS
if ~exist('06_results', 'dir')
    mkdir('06_results');
end

date_str = datestr(now, 'yyyymmdd');
save_filename_1 = sprintf('Test10_2DcamberOmegaSweeps_%s', date_str);
save_filename_2 = sprintf('Test11_2DcamberOmegaVis_%s', date_str);

savefig(fig_sweeps, fullfile('06_results', [save_filename_1, '.fig']));
savefig(fig_vis, fullfile('06_results', [save_filename_2, '.fig']));

%exportgraphics(fig_sweeps, fullfile('06_results', [save_filename_1, '.pdf']), 'ContentType', 'vector');
%exportgraphics(fig_vis, fullfile('06_results', [save_filename_2, '.pdf']), 'ContentType', 'vector');

fprintf('Plots successfully saved to /06_results/\n');