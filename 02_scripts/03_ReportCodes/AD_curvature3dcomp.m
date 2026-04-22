%% 3D CURVATURE ESTIMATION TEST (PARABOLOID & SADDLE)
clear; clc; close all;

% --- 1. SETUP TEST GRID ---
% Create a 50x50 structured grid
[X, Y] = meshgrid(linspace(-1, 1, 50));

% --- 2. DEFINE ANALYTICAL SURFACES ---
% Surface 1: Elliptic Paraboloid (Bowl) -> z = 0.5*x^2 + 0.5*y^2
Z_bowl = 0.5 * X.^2 + 0.5 * Y.^2;

% Surface 2: Hyperbolic Paraboloid (Saddle) -> z = 0.5*x^2 - 0.5*y^2
Z_sad = 0.5 * X.^2 - 0.5 * Y.^2;

% --- 3. CALCULATE NUMERICAL CURVATURE ---
% Run your new function on both surfaces
[k1_bowl, k2_bowl] = Curvature3DTest(X, Y, Z_bowl);
[k1_sad,  k2_sad]  = Curvature3DTest(X, Y, Z_sad);

% Calculate Numerical Maximum Absolute Principal Curvature
k_max_num_bowl = max(abs(k1_bowl), abs(k2_bowl));
k_max_num_sad  = max(abs(k1_sad), abs(k2_sad));

% --- 4. CALCULATE ANALYTICAL MAXIMUM ABSOLUTE PRINCIPAL CURVATURE ---
% Shared geometric denominator: W^2 = 1 + Zx^2 + Zy^2
W2 = 1 + X.^2 + Y.^2;

% 1. Analytical Mean (H) and Gaussian (K) for Bowl (z = 0.5*x^2 + 0.5*y^2)
H_bowl_exact = ( (1 + X.^2).*(1) + (1 + Y.^2).*(1) ) ./ (2 .* W2.^1.5);
K_bowl_exact = (1) ./ (W2.^2);

% 2. Analytical Mean (H) and Gaussian (K) for Saddle (z = 0.5*x^2 - 0.5*y^2)
H_sad_exact = ( (1 + X.^2).*(-1) + (1 + Y.^2).*(1) ) ./ (2 .* W2.^1.5);
K_sad_exact = (-1) ./ (W2.^2);

% 3. Extract Analytical Principal Curvatures (k = H +/- sqrt(H^2 - K))
k1_ana_bowl = H_bowl_exact + sqrt(H_bowl_exact.^2 - K_bowl_exact);
k2_ana_bowl = H_bowl_exact - sqrt(H_bowl_exact.^2 - K_bowl_exact);

k1_ana_sad = H_sad_exact + sqrt(H_sad_exact.^2 - K_sad_exact);
k2_ana_sad = H_sad_exact - sqrt(H_sad_exact.^2 - K_sad_exact);

% 4. Evaluate Maximum Absolute Principal Curvature
k_max_ana_bowl = max(abs(k1_ana_bowl), abs(k2_ana_bowl));
k_max_ana_sad  = max(abs(k1_ana_sad),  abs(k2_ana_sad));

% --- 5. CALCULATE ERRORS ---
% Absolute spatial error for plotting
Err_bowl = abs(k_max_num_bowl - k_max_ana_bowl);
Err_sad  = abs(k_max_num_sad - k_max_ana_sad);

% Root Mean Square Error (RMSE) for the whole surface
RMSE_bowl = sqrt(mean(Err_bowl(:).^2));
RMSE_sad  = sqrt(mean(Err_sad(:).^2));

% --- 6. PLOT RESULTS ---
fig_3d = figure('Name', '3D Curvature Test', 'Units', 'centimeters', 'Position', [5 5, 40, 10]);
tiledlayout(1, 4, "TileSpacing", "compact");

% Color limits to keep visuals consistent
max_K = max(max(k_max_num_bowl(:)), max(k_max_num_sad(:)));
max_E = max(max(Err_bowl(:)), max(Err_sad(:)));

% Tile 1: Bowl Curvature
ax1 = nexttile; hold on; grid on;
surf(ax1, X, Y, Z_bowl, k_max_num_bowl, 'EdgeColor', 'none');
title(ax1, 'Bowl: Max Abs Curvature |\kappa_{max}|');
xlabel(ax1, 'X'); ylabel(ax1, 'Y'); zlabel(ax1, 'Z');
colormap(ax1, 'jet');
clim(ax1, [0, max_K]);
colorbar(ax1);
view(ax1, 3); axis(ax1, 'equal');

% Tile 2: Saddle Curvature
ax2 = nexttile; hold on; grid on;
surf(ax2, X, Y, Z_sad, k_max_num_sad, 'EdgeColor', 'none');
title(ax2, 'Saddle: Max Abs Curvature |\kappa_{max}|');
xlabel(ax2, 'X'); ylabel(ax2, 'Y'); zlabel(ax2, 'Z');
colormap(ax2, 'jet');
clim(ax2, [0, max_K]);
colorbar(ax2);
view(ax2, 3); axis(ax2, 'equal');

% Tile 3: Bowl Absolute Error
ax3 = nexttile; hold on; grid on;
surf(ax3, X, Y, Z_bowl, Err_bowl, 'EdgeColor', 'none');
title(ax3, sprintf('Bowl Error (RMSE: %.2e)', RMSE_bowl));
xlabel(ax3, 'X'); ylabel(ax3, 'Y'); zlabel(ax3, 'Z');
colormap(ax3, 'hot'); % Switch colormap to distinguish error plots
clim(ax3, [0, max_E]);
cb3 = colorbar(ax3);
cb3.Label.String = 'Absolute Error';
view(ax3, 3); axis(ax3, 'equal');

% Tile 4: Saddle Absolute Error
ax4 = nexttile; hold on; grid on;
surf(ax4, X, Y, Z_sad, Err_sad, 'EdgeColor', 'none');
title(ax4, sprintf('Saddle Error (RMSE: %.2e)', RMSE_sad));
xlabel(ax4, 'X'); ylabel(ax4, 'Y'); zlabel(ax4, 'Z');
colormap(ax4, 'hot'); 
clim(ax4, [0, max_E]);
cb4 = colorbar(ax4);
cb4.Label.String = 'Absolute Error';
view(ax4, 3); axis(ax4, 'equal');

% --- APPLY GLOBAL STYLES ---
set(findall(fig_3d, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(fig_3d, '-property', 'FontWeight'), 'FontWeight', 'bold');
set(findall(fig_3d, '-property', 'FontSize'), 'FontSize', 12);
set(findall(fig_3d, 'Type', 'axes'), 'LineWidth', 1.5);

% --- 7. SAVE FIGURES ---
output_dir = '04_figures';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
date_str = datestr(now, 'yyyymmdd');
save_filename = sprintf('F4_3DCurvatureTest_%s', date_str);
save_filepath = fullfile(output_dir, [save_filename, '.fig']);

savefig(fig_3d, save_filepath);
%exportgraphics(fig_3d, fullfile(output_dir, [save_filename, '.pdf']), 'ContentType', 'vector');

fprintf('Figure successfully saved to: %s\n', save_filepath);