%% VARIABLE NAMING CONVENTION                                              
%  Nx       Known (Control) Points                                         
%  Nx_idx   Nx indices in Ax                                               
%  Ax       Aerofoil Mesh Points                                           
%  RAx      Applied Deformation                                            
%  DAx      Actual Aerofoil Position After Rx                              

clear; clc; close all;

%% LOAD 3D MESH       
filename = '05_meshes\01_3d_meshes\surfacepoints140K.plt';
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
Ax = RawData(:, 1:3); % 3D Coordinates [X, Y, Z]
M = size(Ax, 1);

% Base Initial Point - Seed with the Root Trailing Edge (Index 1)                                       
Nx_idx_base = [1];
 
%% DEFORMATION
% Dummy 3D twist along the span for testing
DAx = BendTwist3D(Ax, [24.9,0,-3.28], [47, 35.4, 1.55], 20, -20);
RAx = DAx - Ax;

%% SINGLE RUN EXECUTION
SF_R = 100; 
max_err_threshold = 1e-4; 
N_IP = 100;

fprintf('Running FPS for %d Initial Points...\n', N_IP);
Nx_idx_fps = Nx_idx_base;
Nx_idx_fps = IP_Distance(Nx_idx_fps, Ax, N_IP);

fprintf('Running Greedy Cholesky to E_tol = %g...\n', max_err_threshold);
[Nx_idx_final, max_err_hist, rmse_hist] = GreedyCholesky(Nx_idx_fps, Ax, RAx, SF_R, "err_tol", max_err_threshold);

N_G = length(max_err_hist);
fprintf('Solver complete. Added %d Greedy Points (Total Nodes: %d).\n', N_G, length(Nx_idx_final));

%% DATA EXTRACTION FOR PLOTTING
% Extract spatial coordinates for the FPS points
fps_points = Ax(Nx_idx_fps, :);

% Find greedy points by taking the difference between the final list and the FPS list
greedy_idx = setdiff(Nx_idx_final, Nx_idx_fps, 'stable');
greedy_points = Ax(greedy_idx, :);

%% PLOTTING
fprintf('Generating visualization...\n');
fig_main = figure('Name', 'Single Run Visualization', 'Position', [100, 100, 1200, 500]);
tiledlayout(1, 2, "TileSpacing", "compact");

% --- TILE 1: 3D Deformation ---
ax1 = nexttile;
axes(ax1); % Make ax1 the active axis for the custom function
hold on; grid on; axis equal;
PlotDeformation3D(Ax, DAx);
title('3D Surface Deformation');
%view(3);
hold off;

% --- TILE 2: X-Z Plane Control Nodes ---
ax2 = nexttile;
hold(ax2, 'on'); grid(ax2, 'on'); axis(ax2, 'equal');

% Plot the background surface mesh in light grey (X and Z coordinates)
scatter(ax2, Ax(:,1), Ax(:,2), 2, [0.8 0.8 0.8], 'filled', 'DisplayName', 'Surface Mesh');

% Plot the Initial (FPS) Points in Blue
scatter(ax2, fps_points(:,1), fps_points(:,2), 25, 'b', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('Initial Points (N=%d)', length(Nx_idx_fps)));

% Plot the Greedy Selected Points in Red
if ~isempty(greedy_points)
    scatter(ax2, greedy_points(:,1), greedy_points(:,2), 25, 'r', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('Greedy Points (N=%d)', length(greedy_idx)));
end

xlabel(ax2, 'X Coordinate (Chordwise)');
ylabel(ax2, 'Z Coordinate (Vertical)');
title(ax2, 'Chosen Control Nodes (X-Y Projection)');
legend(ax2, 'Location', 'best');
hold(ax2, 'off');

%% SAVE PLOT
if ~exist('06_results', 'dir')
    mkdir('06_results');
end

date_str = datestr(now, 'yyyymmdd');
save_filename = sprintf('Test5_SingleRun_140K_%s', date_str);

savefig(fig_main, fullfile('06_results', [save_filename, '.fig']));
%exportgraphics(fig_main, fullfile('06_results', [save_filename, '.png']), 'Resolution', 300);

fprintf('Plot saved to /06_results/%s.fig\n', save_filename);