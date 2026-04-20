%% Performance Benchmarking: Standard Greedy vs Recurrence Cholesky
clear; clc; close all;

% --- CONFIGURATION ---
% Setup paths and files (Update these filenames to match yours)
mesh_files = {'\05_meshes\01_3d_meshes\surfacepoints140K.plt', ...
              '\05_meshes\01_3d_meshes\surfacepoints1M.plt'};% ...
              %'\05_meshes\01_3d_meshes\surfacepoints8M.plt'};
num_runs = 1;
num_target_points = 500; % Set your fixed target Nc here
SF_R = 100; % Support radius placeholder

% Initialize storage arrays
mesh_sizes = zeros(length(mesh_files), 1);
avg_time_std = zeros(length(mesh_files), 1);
avg_time_imp = zeros(length(mesh_files), 1);
ops_std = zeros(length(mesh_files), 1);
ops_imp = zeros(length(mesh_files), 1);


%% MAIN BENCHMARKING LOOP
for i = 1:length(mesh_files)
    fprintf('Loading %s...\n', mesh_files{i});
    fprintf('Loading numeric data...\n');
    RawData = readmatrix(mesh_files{i}, 'NumHeaderLines', 1, 'FileType', 'text');
    Ax = RawData(:, 1:3); % 3D Coordinates [X, Y, Z]
    M = size(Ax, 1);
    fprintf('Loaded %s...\n', mesh_files{i});
    
    
    % --- PLACEHOLDER DEFORMATION ---
    % Example: Simple bending in the Z-axis based on X position
    DAx = BendTwist3D(Ax, [25,0,-3.3], [47, 35.4, 1.6], 30, -20);
    RAx = DAx - Ax; 
    
    
    mesh_sizes(i) = M;
    % Select an initial point (e.g., the origin or tip)
    Nx_idx_in = 1; 
    
    %% TIMING LOOP: STANDARD GREEDY
    fprintf('  Running Standard Greedy (%d iterations)...\n', num_runs);
    times_std = zeros(num_runs, 1);
    for run = 1:num_runs
        tic;
        [Nx_out_std, ~, ~] = Greedy(Nx_idx_in, Ax, RAx, SF_R, 'err_tol', 0, 'K', num_target_points);
        times_std(run) = toc;
    end
    avg_time_std(i) = mean(times_std);
    
    %% TIMING LOOP: RECURRENCE CHOLESKY
    fprintf('  Running Recurrence Cholesky (%d iterations)...\n', num_runs);
    times_imp = zeros(num_runs, 1);
    for run = 1:num_runs
        tic;
        [Nx_out_imp, ~, ~] = GreedyCholesky(Nx_idx_in, Ax, RAx, SF_R, 'err_tol', 0, 'K', num_target_points);
        times_imp(run) = toc;
    end
    avg_time_imp(i) = mean(times_imp);
    
    %% CALCULATE THEORETICAL OPERATIONS
    % Using the integrated sum over all k steps (k = 1 to Nc)
    k = 1:num_target_points;
    
    % Standard: Full Factorization (k^3/3) + Solve (2k^2) + Eval (M-k)*k
    ops_std(i) = sum((k.^3)/3 + 2*(k.^2) + (M - k).*k);
    
    % Improved: Incremental Update (k^2) + Solve (2k^2) + Eval (M-k)*k
    ops_imp(i) = sum(k.^2 + 2*(k.^2) + (M - k).*k);
    
    fprintf('  Done. Standard: %.2fs | Improved: %.2fs\n', avg_time_std(i), avg_time_imp(i));
end

%% PLOT 1: VERIFICATION OF CONTROL POINTS (Final Mesh)
% Verifies that both algorithms yield the exact same spatial points
fig_verify = figure('Name', 'Control Point Verification');
scatter3(Ax(:,1), Ax(:,2), Ax(:,3), 1, [0.8 0.8 0.8], '.'); hold on;
scatter3(Ax(Nx_out_std, 1), Ax(Nx_out_std, 2), Ax(Nx_out_std, 3), 50, 'ro', 'LineWidth', 1.5);
scatter3(Ax(Nx_out_imp, 1), Ax(Nx_out_imp, 2), Ax(Nx_out_imp, 3), 50, 'b+', 'LineWidth', 1.5);
title(sprintf('Selected Points on %d-Node Mesh', mesh_sizes(end)));
legend('Surface Mesh', 'Standard Greedy', 'Recurrence Cholesky');
axis equal; grid on; view(3);

%% PLOT 2: MAIN PERFORMANCE RESULTS
fig_main = figure('Name', 'Performance Benchmarking', 'Position', [100, 100, 1000, 400]);

% Left Subplot: Wall-Clock Time
subplot(1, 2, 1);
loglog(mesh_sizes, avg_time_std, '-ro', 'LineWidth', 2, 'MarkerSize', 6); hold on;
loglog(mesh_sizes, avg_time_imp, '-b^', 'LineWidth', 2, 'MarkerSize', 6);
title('Empirical Wall-Clock Execution Time');
xlabel('Number of Surface Mesh Nodes (M)');
ylabel('Average Time (Seconds)');
legend('Standard Greedy', 'Recurrence Cholesky', 'Location', 'northwest');
grid on;

% Right Subplot: Theoretical Operations
subplot(1, 2, 2);
loglog(mesh_sizes, ops_std, '-ro', 'LineWidth', 2, 'MarkerSize', 6); hold on;
loglog(mesh_sizes, ops_imp, '-b^', 'LineWidth', 2, 'MarkerSize', 6);
title('Theoretical Floating-Point Operations');
xlabel('Number of Surface Mesh Nodes (M)');
ylabel('Total Operations (FLOPs)');
legend('Standard Greedy', 'Recurrence Cholesky', 'Location', 'northwest');
grid on;

%% SAVE DATA
date_str = datestr(now, 'yyyymmdd');
save_filename = sprintf('Test1ImprovedTiming_%s', date_str);

% Save the main plot as a .fig file for future editing
savefig(fig_main, fullfile('04_figures', [save_filename, '.fig']));

% Optional: Also save as a high-res PDF for Overleaf
%exportgraphics(fig_main, fullfile('plots', [save_filename, '.pdf']), 'ContentType', 'vector');

fprintf('\nBenchmarking complete. Figure saved to /plots/%s.fig\n', save_filename);