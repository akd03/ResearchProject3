%% NACA 0012 2D CURVATURE ESTIMATION (RAW VS SMOOTHED)
clear; clc; close all;

%% 1. NACA 0012 ANALYTICAL FUNCTIONS
NACA = @(x) 5 * 0.12 * (0.2969.*x.^0.5 - 0.1260.*x - 0.3516.*x.^2 + 0.2843.*x.^3 - 0.1036.*x.^4);
dNACAdx = @(x) 5 * 0.12 * (0.5*0.2969.*x.^-0.5 - 0.1260 - 2*0.3516.*x + 3*0.2843.*x.^2 - 4*0.1036.*x.^3);
d2NACAdx2 = @(x) 5 * 0.12 * (-0.25*0.2969.*x.^-1.5 - 2*0.3516 + 6*0.2843.*x - 12*0.1036.*x.^2);

%% 2. LOAD TRUE NACA MESH
MeshFile = fopen("05_meshes\NACA0012257x129.xyz", "r");
Header = fscanf(MeshFile, '%d', 2);
Ni = Header(1);
Nj = Header(2);
MeshData = textscan(MeshFile, '%f %f %f');
fclose(MeshFile);

% Extract Aerofoil Points
Mx = [MeshData{1}, MeshData{3}];
Ax = Mx(1:Ni, 1:2); 

%% 3. NUMERICAL CURVATURE (CDS)
% Calculate the full loop twice: without smoothing (gamma = 0) and with smoothing (gamma = 0.25)
[grad1_raw, grad2_raw]       = CentralDifferenceDeriv(Ax, 0);
[grad1_smooth, grad2_smooth] = CentralDifferenceDeriv(Ax, 0.25);

% Numerical kappa using the parametric equation
kappa_num_raw    = abs(grad1_raw(:,1) .* grad2_raw(:,2) - grad1_raw(:,2) .* grad2_raw(:,1)) ./ ...
                   ((grad1_raw(:,1).^2 + grad1_raw(:,2).^2).^1.5);
               
kappa_num_smooth = abs(grad1_smooth(:,1) .* grad2_smooth(:,2) - grad1_smooth(:,2) .* grad2_smooth(:,1)) ./ ...
                   ((grad1_smooth(:,1).^2 + grad1_smooth(:,2).^2).^1.5);

%% 4. ANALYTICAL CURVATURE (TOP SURFACE EXACT MATCH)
% Extract upper surface (Leading Edge to Trailing Edge)
idx_top = 129:257;
x_top = Ax(idx_top, 1);

% Evaluate analytical equations at the EXACT x-coordinates of the mesh
y_prime = dNACAdx(x_top);
y_double_prime = d2NACAdx2(x_top);

% Analytical kappa
kappa_ana = abs(y_double_prime) ./ (1 + y_prime.^2).^1.5;

%% 5. RELATIVE ERROR CALCULATION
rel_err_raw    = abs(kappa_num_raw(idx_top) - kappa_ana) ./ kappa_ana;
rel_err_smooth = abs(kappa_num_smooth(idx_top) - kappa_ana) ./ kappa_ana;

%% 6. RADIUS OF CURVATURE CALCULATION
% Numerical Radius
r_num_raw    = 1 ./ kappa_num_raw(idx_top);
r_num_smooth = 1 ./ kappa_num_smooth(idx_top);

% Analytical Radius with LE Singularity Override (x = 0)
r_ana = 1 ./ kappa_ana;
t = 0.12;
r_ana(1) = 1.1019 * t^2; 

%% 7. PLOT COMPARISON
fig_curv = figure('Name', 'Curvature Estimation', 'Units', 'centimeters', 'Position', [5 5, 25, 10]);
tiledlayout(1, 2, "TileSpacing", "compact");

% --- TILE 1: RADIUS OF CURVATURE ---
ax1 = nexttile; hold on; grid on;
plot(x_top, r_ana, "r.-", "LineWidth", 2, "DisplayName", "Analytical");
plot(x_top, r_num_raw, "g.-", "LineWidth", 1.5, "DisplayName", "Numerical (Raw)");
plot(x_top, r_num_smooth, "b.-", "LineWidth", 1.5, "DisplayName", "Numerical (Smoothed)");

xlabel('Chord (x)');
ylabel('Radius of Curvature (r)');
title('Radius of Curvature Comparison');
legend('Location', 'northwest');
ylim([0, 8.5]);
xlim([0, 1]);
set(ax1, "Layer", "Top");

% --- TILE 2: NORMALIZED ERROR ---
ax2 = nexttile; hold on; grid on;
plot(x_top, rel_err_raw*100, "g.-", "LineWidth", 1.5, "DisplayName", "Raw Error");
plot(x_top, rel_err_smooth*100, "b.-", "LineWidth", 1.5, "DisplayName", "Smoothed Error");

xlabel('Chord (x)');
ylabel('Percentage Error (%)');
title('Normalized Curvature Error');
legend('Location', 'best');
xlim([0, 1]);

% Set reasonable y-limits in case the raw error spikes excessively at the LE
max_smooth_err = max(rel_err_smooth*100);
if max_smooth_err < 100
    ylim([0, max_smooth_err * 2]); % Buffer the top of the plot
end
set(ax2, "Layer", "Top");

% --- APPLY GLOBAL STYLES ---
set(findall(fig_curv, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(fig_curv, '-property', 'FontWeight'), 'FontWeight', 'bold');
set(findall(fig_curv, '-property', 'FontSize'), 'FontSize', 12);
set(findall(fig_curv, 'Type', 'axes'), 'LineWidth', 2);

%% 8. SAVE FIGURES
output_dir = '04_figures';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

date_str = datestr(now, 'yyyymmdd');
save_filename = sprintf('F3_2DCurvatureEstimation_%s', date_str);
save_filepath = fullfile(output_dir, [save_filename, '.fig']);

savefig(fig_curv, save_filepath);
%exportgraphics(fig_curv, fullfile(output_dir, [save_filename, '.pdf']), 'ContentType', 'vector');

fprintf('Figure successfully saved to: %s\n', save_filepath);