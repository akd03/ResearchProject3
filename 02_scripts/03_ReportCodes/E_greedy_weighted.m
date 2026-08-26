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

% Extract Full Mesh and Surface Points
Mx = [reshape(MeshData{1}, [Ni*Nj, 1]), reshape(MeshData{3}, [Ni*Nj, 1])];
Ax = [Mx(1:Ni, 1), Mx(1:Ni, 2)];
M = size(Ax, 1);

X_base = reshape(Mx(:,1), Ni, Nj);
Y_base = reshape(Mx(:,2), Ni, Nj);

Nx_idx_base = [1]; % Seed with Trailing Edge

%% 2. DEFORMATIONS
% Deformation 1: Rotation
DMx_rot = Rotate(Mx, 30, 0.5);
DAx_rot = DMx_rot(1:Ni, :);
RAx_rot = DAx_rot - Ax;

% Deformation 2: Camber
DMx_cam = NACACamber(Mx, 0.15, 0.5); % Applying to full mesh for orthogonality
DAx_cam = DMx_cam(1:Ni, :);
RAx_cam = DAx_cam - Ax;

%% 3. ORTHOGONALITY CALCULATIONS
Q_base = calculateOrthogonality(X_base, Y_base);

X_rot = reshape(DMx_rot(:,1), Ni, Nj);
Y_rot = reshape(DMx_rot(:,2), Ni, Nj);
Q_rot = calculateOrthogonality(X_rot, Y_rot);
dQ_rot = Q_rot - Q_base;

X_cam = reshape(DMx_cam(:,1), Ni, Nj);
Y_cam = reshape(DMx_cam(:,2), Ni, Nj);
Q_cam = calculateOrthogonality(X_cam, Y_cam);
dQ_cam = Q_cam - Q_base;

% Output Table
fprintf('\n==================================================\n');
fprintf('ORTHOGONALITY ANALYSIS (Target Analytical Meshes)\n');
fprintf('==================================================\n');
fprintf('Base Mesh Average Orthogonality:    %.6f\n', mean(Q_base(:)));
fprintf('Rotation Average Orthogonality:     %.6f\n', mean(Q_rot(:)));
fprintf('Camber Average Orthogonality:       %.6f\n', mean(Q_cam(:)));
fprintf('==================================================\n\n');

%% 4. CALCULATE NORMALIZED CURVATURE (CENTRAL DIFFERENCES)
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

%% 5. PARAMETER SWEEP SETUP
SF_R = 3;                  
omega_vals = 0:0.01:1.0;    
num_omega = length(omega_vals);

N_vals = [25, 50, 100];
num_N_tests = length(N_vals);
err_threshold = 1e-6;

% Storage Arrays (Rot = Rotation, Cam = Camber)
err_rot = zeros(num_N_tests, num_omega);
rmse_rot = zeros(num_N_tests, num_omega);
nodes_rot = zeros(1, num_omega);

err_cam = zeros(num_N_tests, num_omega);
rmse_cam = zeros(num_N_tests, num_omega);
nodes_cam = zeros(1, num_omega);

vis_N = 25;
vis_omega = [0, 0.4, 0.75];
vis_idx_rot = cell(3, 1);
vis_idx_cam = cell(3, 1);

%% 6. DATA GENERATION LOOP
fprintf('Generating Data for Curvature Sweeps...\n');
for j = 1:num_omega
    omega = omega_vals(j);
    
    % SWEEP 1: Fixed Nodes
    for i = 1:num_N_tests
        N_G = N_vals(i) - length(Nx_idx_base); 
        
        [~, max_err_r, rmse_r] = GreedyCholesky_Weighted(Nx_idx_base, Ax, RAx_rot, SF_R, kappa_norm, omega, 'K', N_G);
        err_rot(i, j) = max_err_r(end);
        rmse_rot(i, j) = rmse_r(end);
        
        [~, max_err_c, rmse_c] = GreedyCholesky_Weighted(Nx_idx_base, Ax, RAx_cam, SF_R, kappa_norm, omega, 'K', N_G);
        err_cam(i, j) = max_err_c(end);
        rmse_cam(i, j) = rmse_c(end);
    end
    
    % SWEEP 2: Fixed Error Threshold
    [idx_out_r, ~, ~] = GreedyCholesky_Weighted(Nx_idx_base, Ax, RAx_rot, SF_R, kappa_norm, omega, 'err_tol', err_threshold, 'K', Inf);
    nodes_rot(j) = length(idx_out_r);
    
    [idx_out_c, ~, ~] = GreedyCholesky_Weighted(Nx_idx_base, Ax, RAx_cam, SF_R, kappa_norm, omega, 'err_tol', err_threshold, 'K', Inf);
    nodes_cam(j) = length(idx_out_c);
end

% SWEEP 3: Visualization Points
N_G_vis = vis_N - length(Nx_idx_base);
for k = 1:3
    [idx_r, ~, ~] = GreedyCholesky_Weighted(Nx_idx_base, Ax, RAx_rot, SF_R, kappa_norm, vis_omega(k), 'K', N_G_vis);
    vis_idx_rot{k} = idx_r;
    
    [idx_c, ~, ~] = GreedyCholesky_Weighted(Nx_idx_base, Ax, RAx_cam, SF_R, kappa_norm, vis_omega(k), 'K', N_G_vis);
    vis_idx_cam{k} = idx_c;
end

%% 7. PLOT 1: PERFORMANCE SWEEPS (1x3)
fig_sweeps = figure('Name', 'Deformation Sweeps', 'Units', 'centimeters', 'Position', [5, 15, 17, 6]);
tiledlayout(1, 3, "TileSpacing", "compact", "Padding", "compact");
colors = colororder;

% Tile 1: Max Error
ax1 = nexttile;
hold(ax1, 'on'); grid(ax1, 'on'); set(ax1, 'YScale', 'log', 'LineWidth', 1.2, 'FontSize', 10, 'FontWeight', 'bold');
xlabel(ax1, 'Curvature Weighting (\omega)'); ylabel(ax1, 'Max Error');
for i = 1:num_N_tests
    c_idx = mod(i-1, size(colors,1)) + 1;
    plot(ax1, omega_vals, err_rot(i, :), '-', 'LineWidth', 1.5, 'Color', colors(c_idx, :), 'DisplayName', sprintf('Rot: N=%d', N_vals(i)));
    plot(ax1, omega_vals, err_cam(i, :), '--', 'LineWidth', 1.5, 'Color', colors(c_idx, :), 'DisplayName', sprintf('Cam: N=%d', N_vals(i)));
end
legend(ax1, 'Location', 'northwest', 'FontSize', 8, 'FontWeight', 'normal');

% Tile 2: RMSE
ax2 = nexttile;
hold(ax2, 'on'); grid(ax2, 'on'); set(ax2, 'YScale', 'log', 'LineWidth', 1.2, 'FontSize', 10, 'FontWeight', 'bold');
xlabel(ax2, 'Curvature Weighting (\omega)'); ylabel(ax2, 'RMSE');
for i = 1:num_N_tests
    c_idx = mod(i-1, size(colors,1)) + 1;
    plot(ax2, omega_vals, rmse_rot(i, :), '-', 'LineWidth', 1.5, 'Color', colors(c_idx, :), 'DisplayName', sprintf('Rot: N=%d', N_vals(i)));
    plot(ax2, omega_vals, rmse_cam(i, :), '--', 'LineWidth', 1.5, 'Color', colors(c_idx, :), 'DisplayName', sprintf('Cam: N=%d', N_vals(i)));
end
legend(ax2, 'Location', 'northwest', 'FontSize', 8, 'FontWeight', 'normal');

% Tile 3: Total Nodes
ax3 = nexttile;
hold(ax3, 'on'); grid(ax3, 'on'); set(ax3, 'LineWidth', 1.2, 'FontSize', 10, 'FontWeight', 'bold');
xlabel(ax3, 'Curvature Weighting (\omega)'); ylabel(ax3, 'Total Nodes');
plot(ax3, omega_vals, nodes_rot, '-', 'LineWidth', 1.5, 'Color', 'k', 'DisplayName', 'Rotation');
plot(ax3, omega_vals, nodes_cam, '--', 'LineWidth', 1.5, 'Color', 'k', 'DisplayName', 'Camber');
legend(ax3, 'Location', 'northwest', 'FontSize', 8, 'FontWeight', 'normal');

%% 8. PLOT 2: ROTATION VISUALIZATIONS (1x3)
% Bulletproof padding for pcolor: make C exactly match X and Y sizes with NaNs at the edges
dQ_rot_pad = nan(size(X_rot)); 
dQ_rot_pad(1:size(dQ_rot,1), 1:size(dQ_rot,2)) = dQ_rot;

fig_rot = figure('Name', 'Rotation Vis', 'Units', 'centimeters', 'Position', [5, 8, 17, 6]);
tiledlayout(1, 3, "TileSpacing", "tight", "Padding", "tight");

for k = 1:3
    ax = nexttile;
    hold(ax, 'on'); set(ax, 'LineWidth', 1.2, 'FontSize', 10, 'FontWeight', 'bold');
    xlim([-2, 3]); ylim([-1.5, 1.5]);
    % Orthogonality Background
    pc = pcolor(ax, X_rot, Y_rot, dQ_rot_pad);
    shading(ax, 'flat'); colormap(ax, 'jet');
    set(pc, 'HandleVisibility', 'off'); 
    
    % Profile and Points
    plot(ax, DAx_rot(:,1), DAx_rot(:,2), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Deformed Aerofoil');
    c_points = DAx_rot(vis_idx_rot{k}, :);
    scatter(ax, c_points(:,1), c_points(:,2), 30, 'w', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Control Nodes');
    
    xlabel(ax, 'X'); ylabel(ax, 'Y');
  
    
    cb = colorbar(ax); 
    if k == 3
        ylabel(cb, '\Delta Orthogonality (q)', 'FontSize', 10, 'FontWeight', 'bold'); 
    else
        set(cb, 'Visible', 'off');
    end
end

%% 9. PLOT 3: CAMBER VISUALIZATIONS (1x3)
% Bulletproof padding for pcolor
dQ_cam_pad = nan(size(X_cam)); 
dQ_cam_pad(1:size(dQ_cam,1), 1:size(dQ_cam,2)) = dQ_cam;

fig_cam = figure('Name', 'Camber Vis', 'Units', 'centimeters', 'Position', [5, 1, 17, 6]);
tiledlayout(1, 3, "TileSpacing", "compact", "Padding", "compact");

for k = 1:3
    ax = nexttile;
    hold(ax, 'on'); set(ax, 'LineWidth', 1.2, 'FontSize', 10, 'FontWeight', 'bold');
    xlim([-2, 3]); ylim([-1.5, 1.5]);
    % Orthogonality Background
    pc = pcolor(ax, X_cam, Y_cam, dQ_cam_pad);
    shading(ax, 'flat'); colormap(ax, 'jet');
    set(pc, 'HandleVisibility', 'off'); 
    
    % Profile and Points
    plot(ax, DAx_cam(:,1), DAx_cam(:,2), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Deformed Aerofoil');
    c_points = DAx_cam(vis_idx_cam{k}, :);
    scatter(ax, c_points(:,1), c_points(:,2), 30, 'w', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Control Nodes');
    
    xlabel(ax, 'X'); ylabel(ax, 'Y');
    
   
    cb = colorbar(ax); 
    if k == 3
        ylabel(cb, '\Delta Orthogonality (q)', 'FontSize', 8, 'FontWeight', 'bold'); 
    else
        set(cb, 'Visible', 'off');
    end
end

%% 10. SAVE PLOTS
if ~exist('06_results', 'dir'), mkdir('06_results'); end
date_str = datestr(now, 'yyyymmdd');

% Save Sweeps
save_name1 = sprintf('Test20_Sweeps_%s', date_str);
savefig(fig_sweeps, fullfile('06_results', [save_name1, '.fig']));
exportgraphics(fig_sweeps, fullfile('06_results', [save_name1, '.pdf']), 'ContentType', 'vector');

% Save Rotation
save_name2 = sprintf('Test20_RotVis_%s', date_str);
savefig(fig_rot, fullfile('06_results', [save_name2, '.fig']));
exportgraphics(fig_rot, fullfile('06_results', [save_name2, '.pdf']), 'ContentType', 'vector');

% Save Camber
save_name3 = sprintf('Test20_CamVis_%s', date_str);
savefig(fig_cam, fullfile('06_results', [save_name3, '.fig']));
exportgraphics(fig_cam, fullfile('06_results', [save_name3, '.pdf']), 'ContentType', 'vector');

fprintf('All three figures successfully saved and exported as PDFs.\n');

%% HELPER FUNCTIONS
function Q = calculateOrthogonality(X, Y)
    [N, M] = size(X);
    Q = zeros(N-1, M-1);
    for i = 1:N-1
        for j = 1:M-1
            x1 = X(i, j);     y1 = Y(i, j);
            x2 = X(i+1, j);   y2 = Y(i+1, j);
            x3 = X(i+1, j+1); y3 = Y(i+1, j+1);
            x4 = X(i, j+1);   y4 = Y(i, j+1);
            
            v1 = [x2 - x1, y2 - y1];
            v2 = [x3 - x2, y3 - y2];
            v3 = [x4 - x3, y4 - y3];
            v4 = [x1 - x4, y1 - y4];
            
            v1_sq = dot(v1, v1); v2_sq = dot(v2, v2);
            v3_sq = dot(v3, v3); v4_sq = dot(v4, v4);
            
            if v1_sq == 0 || v2_sq == 0 || v3_sq == 0 || v4_sq == 0
                Q(i, j) = 1; continue;
            end
            
            Q(i, j) = 0.25 * ((dot(v1, v2)^2)/(v1_sq*v2_sq) + (dot(v2, v3)^2)/(v2_sq*v3_sq) + ...
                              (dot(v3, v4)^2)/(v3_sq*v4_sq) + (dot(v4, v1)^2)/(v4_sq*v1_sq));
        end
    end
end

function DPoints = Rotate(Points, angle_deg, x_origin)
    theta = deg2rad(angle_deg);
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    DPoints = Points;
    DPoints(:,1) = DPoints(:,1) - x_origin;
    DPoints = (R * DPoints')';
    DPoints(:,1) = DPoints(:,1) + x_origin;
end